# Copyright 2022 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import datetime
import os
import subprocess

import yaml
from cached_property import cached_property
from ebi_eva_common_pyutils.command_utils import run_command_with_output
from ebi_eva_common_pyutils.config import cfg
from ebi_eva_common_pyutils.config_utils import get_contig_alias_db_creds_for_profile
from ebi_eva_common_pyutils.contig_alias.contig_alias import ContigAliasClient
from ebi_eva_common_pyutils.logger import AppLogger
from ebi_eva_common_pyutils.metadata_utils import get_metadata_connection_handle, insert_new_assembly_and_taxonomy
from ebi_eva_common_pyutils.pg_utils import execute_query, get_all_results_for_query
from ebi_eva_common_pyutils.spring_properties import SpringPropertiesGenerator
from ebi_eva_common_pyutils.taxonomy.taxonomy import get_scientific_name_from_taxonomy
from psycopg2.extras import execute_values

from eva_assembly_ingestion.parse_counts import count_variants_extracted, count_variants_remapped, \
    count_variants_ingested


def pretty_print(header, table):
    cell_widths = [len(h) for h in header]
    for row in table:
        for i, cell in enumerate(row):
            cell_widths[i] = max(cell_widths[i], len(str(cell)))
    format_string = ' | '.join('{%s:>%s}' % (i, w) for i, w in enumerate(cell_widths))
    print('| ' + format_string.format(*header) + ' |')
    for row in table:
        print('| ' + format_string.format(*row) + ' |')


class AssemblyIngestionJob(AppLogger):
    all_tasks = ['load_tracker', 'remap_cluster', 'update_dbs']
    tracking_table = 'eva_progress_tracker.remapping_tracker'

    def __init__(self, taxonomy, target_assembly, release_version):
        self.taxonomy = taxonomy
        self.target_assembly = target_assembly
        self.release_version = release_version
        self.private_settings_file = cfg['maven']['settings_file']
        self.maven_profile = cfg['maven']['environment']
        self.properties_generator = SpringPropertiesGenerator(self.maven_profile, self.private_settings_file)

    @cached_property
    def scientific_name(self):
        return get_scientific_name_from_taxonomy(self.taxonomy)

    def run_all(self, tasks, instance, source_of_assembly, resume):
        if 'load_tracker' in tasks:
            self.load_tracker()
        if 'remap_cluster' in tasks:
            self.run_remapping_and_clustering(instance, resume)
        if 'update_dbs' in tasks:
            self.update_dbs(source_of_assembly)

    def load_tracker(self):
        """Load the tracking table with the source assemblies for this taxonomy. Will not load anything if jobs in
        the tracker already exist for this taxonomy/target assembly pair."""
        header_to_print = (
            'Sources', 'Taxonomy', 'Scientific Name', 'Assembly', 'Target Assembly', 'Num Studies', 'Status')
        existing_jobs = self.get_job_information_from_tracker()
        if existing_jobs:
            self.info(f'Jobs already exist for taxonomy {self.taxonomy} and target assembly {self.target_assembly}, '
                      f'not loading anything new.')
            pretty_print(header_to_print, existing_jobs)
            return

        column_names = ('source', 'taxonomy', 'scientific_name', 'origin_assembly_accession', 'assembly_accession',
                        'remapping_version', 'release_version', 'num_studies', 'num_ss_ids', 'study_accessions',
                        'remapping_status')
        rows = []
        rows_to_print = []
        for source_assembly, projects in self.get_source_assemblies_and_projects():
            rows.append(('EVA', self.taxonomy, self.scientific_name, source_assembly, self.target_assembly,
                         1, self.release_version, len(projects), 1, None, 'Pending'))
            rows_to_print.append(('EVA', self.taxonomy, self.scientific_name, source_assembly, self.target_assembly,
                                  len(projects), 'Pending'))
        for source_assembly, num_studies in self.get_source_assemblies_and_num_studies_dbsnp():
            rows.append(('DBSNP', self.taxonomy, self.scientific_name, source_assembly, self.target_assembly,
                         1, self.release_version, num_studies, 1, None, 'Pending'))
            rows_to_print.append(('DBSNP', self.taxonomy, self.scientific_name, source_assembly, self.target_assembly,
                                  num_studies, 'Pending'))
        if len(rows) == 0:
            self.info(f'Nothing to process for taxonomy {self.taxonomy} and target assembly {self.target_assembly}')
            return
        with get_metadata_connection_handle(self.maven_profile, self.private_settings_file) as pg_conn:
            with pg_conn.cursor() as cursor:
                insert_query = f"INSERT INTO {self.tracking_table} ({','.join(column_names)}) VALUES %s "
                execute_values(cursor, insert_query, rows)
        pretty_print(header_to_print, rows_to_print)

    def get_job_information_from_tracker(self):
        """Gets jobs from tracker by target assembly, taxonomy, and release version"""
        query = (
            f"SELECT source, taxonomy, scientific_name, origin_assembly_accession, assembly_accession, "
            f"num_studies, remapping_status "
            f"FROM eva_progress_tracker.remapping_tracker "
            f"WHERE release_version={self.release_version} "
            f"AND assembly_accession='{self.target_assembly}' AND taxonomy={self.taxonomy}"
        )
        with get_metadata_connection_handle(self.maven_profile, self.private_settings_file) as pg_conn:
            return get_all_results_for_query(pg_conn, query)

    def get_source_assemblies_and_projects(self):
        """Query metadata for all public projects with this taxonomy, of these getting all reference accessions
        for all analyses."""
        with get_metadata_connection_handle(self.maven_profile, self.private_settings_file) as pg_conn:
            query = (
                f"SELECT DISTINCT vcf_reference_accession, ARRAY_AGG(project_accession) "
                f"FROM evapro.project "
                f"LEFT OUTER JOIN evapro.project_taxonomy USING (project_accession) "
                f"LEFT OUTER JOIN evapro.project_analysis USING (project_accession) "
                f"LEFT OUTER JOIN evapro.analysis USING (analysis_accession) "
                f"WHERE taxonomy_id={self.taxonomy} AND ena_status=4 AND hidden_in_eva=0 "
                f"GROUP BY vcf_reference_accession"
            )
            return get_all_results_for_query(pg_conn, query)

    def get_source_assemblies_and_num_studies_dbsnp(self):
        # Source assemblies for dbSNP are not in metadata, so instead we get them from release 3 in the tracker
        with get_metadata_connection_handle(self.maven_profile, self.private_settings_file) as pg_conn:
            query = (
                f"SELECT origin_assembly_accession, num_studies FROM {self.tracking_table} "
                f"WHERE release_version=3 AND taxonomy={self.taxonomy} AND source='DBSNP'"
            )
            return get_all_results_for_query(pg_conn, query)

    def run_remapping_and_clustering(self, instance, resume):
        """Run remapping and clustering for all source assemblies in the tracker marked as not Complete, resuming
        the nextflow process if specified. (Note that this will also resume or rerun anything marked as Failed.)"""
        source_assemblies = self.get_incomplete_assemblies()
        self.info(f'Running remapping and clustering for the following assemblies: {source_assemblies}')
        for source_assembly in source_assemblies:
            self.process_one_assembly(source_assembly, instance, resume)

    def get_incomplete_assemblies(self):
        incomplete_assemblies = []
        for row in self.get_job_information_from_tracker():
            source_assembly = row[3]
            status = row[6]
            if status != 'Completed':
                incomplete_assemblies.append(source_assembly)
        return incomplete_assemblies

    def process_one_assembly(self, source_assembly, instance, resume):
        self.set_status_start(source_assembly)
        base_directory = cfg['remapping']['base_directory']
        nextflow_pipeline = os.path.join(os.path.dirname(__file__), 'nextflow', 'remap_cluster.nf')
        assembly_directory = os.path.join(base_directory, self.taxonomy, source_assembly)
        work_dir = os.path.join(assembly_directory, 'work')

        extraction_properties_file = self.create_extraction_properties(
            output_file_path=os.path.join(assembly_directory, 'remapping_extraction.properties'),
            source_assembly=source_assembly
        )
        ingestion_properties_file = self.create_ingestion_properties(
            output_file_path=os.path.join(assembly_directory, 'remapping_ingestion.properties'),
            source_assembly=source_assembly
        )
        clustering_template_file = self.create_clustering_properties(
            output_file_path=os.path.join(assembly_directory, 'clustering_template.properties'),
            instance=instance,
            source_assembly=source_assembly
        )

        os.makedirs(work_dir, exist_ok=True)
        remapping_log = os.path.join(assembly_directory, 'remapping_process.log')
        remap_cluster_config_file = os.path.join(assembly_directory, 'remap_cluster_config.yaml')
        remap_cluster_config = {
            'taxonomy_id': self.taxonomy,
            'source_assembly_accession': source_assembly,
            'target_assembly_accession': self.target_assembly,
            'species_name': self.scientific_name,
            'output_dir': assembly_directory,
            'genome_assembly_dir': cfg['genome_downloader']['output_directory'],
            'extraction_properties': extraction_properties_file,
            'ingestion_properties': ingestion_properties_file,
            'clustering_properties': clustering_template_file,
            'clustering_instance': instance,
            'remapping_config': cfg.config_file,
            'remapping_required': self.check_remapping_required(source_assembly)
        }

        for part in ['executable', 'nextflow', 'jar']:
            remap_cluster_config[part] = cfg[part]
        with open(remap_cluster_config_file, 'w') as open_file:
            yaml.safe_dump(remap_cluster_config, open_file)
        try:
            command = [
                cfg['executable']['nextflow'],
                '-log', remapping_log,
                'run', nextflow_pipeline,
                '-params-file', remap_cluster_config_file,
                '-work-dir', work_dir
            ]
            if resume:
                command.append('-resume')
            curr_working_dir = os.getcwd()
            os.chdir(assembly_directory)
            run_command_with_output('Nextflow remapping process', ' '.join(command))
        except subprocess.CalledProcessError as e:
            self.error('Nextflow remapping pipeline failed')
            self.set_status_failed(source_assembly)
            raise e
        finally:
            os.chdir(curr_working_dir)
        self.set_status_end(source_assembly)
        self.count_variants_from_logs(assembly_directory, source_assembly)

    def check_remapping_required(self, source_assembly):
        return source_assembly != self.target_assembly

    def create_extraction_properties(self, output_file_path, source_assembly):
        properties = self.properties_generator.get_remapping_extraction_properties(
            taxonomy=self.taxonomy,
            source_assembly=source_assembly
        )
        with open(output_file_path, 'w') as open_file:
            open_file.write(properties)
        return output_file_path

    def create_ingestion_properties(self, output_file_path, source_assembly):
        properties = self.properties_generator.get_remapping_ingestion_properties(
            source_assembly=source_assembly,
            target_assembly=self.target_assembly
        )
        with open(output_file_path, 'w') as open_file:
            open_file.write(properties)
        return output_file_path

    def create_clustering_properties(self, output_file_path, instance, source_assembly):
        properties = self.properties_generator.get_clustering_properties(
            instance=instance,
            source_assembly=source_assembly,
            target_assembly=self.target_assembly,
            rs_report_path=f'{source_assembly}_to_{self.target_assembly}_rs_report.txt'
        )
        with open(output_file_path, 'w') as open_file:
            open_file.write(properties)
        return output_file_path

    def set_status(self, source_assembly, status):
        query = (
            f"UPDATE {self.tracking_table} "
            f"SET remapping_status='{status}', remapping_start = '{datetime.now().isoformat()}' "
            f"WHERE release_version={self.release_version} "
            f"AND origin_assembly_accession='{source_assembly}' AND taxonomy={self.taxonomy}"
        )
        with get_metadata_connection_handle(self.maven_profile, self.private_settings_file) as pg_conn:
            execute_query(pg_conn, query)

    def set_status_start(self, source_assembly):
        self.set_status(source_assembly, 'Started')

    def set_status_end(self, source_assembly):
        self.set_status(source_assembly, 'Completed')

    def set_status_failed(self, source_assembly):
        self.set_status(source_assembly, 'Failed')

    def set_counts(self, source_assembly, source, nb_variant_extracted=None, nb_variant_remapped=None,
                   nb_variant_ingested=None):
        set_statements = []
        query = (
            f"SELECT * FROM {self.tracking_table} "
            f"WHERE release_version={self.release_version} AND origin_assembly_accession='{source_assembly}' "
            f"AND taxonomy='{self.taxonomy}' AND source='{source}'"
        )
        with get_metadata_connection_handle(self.maven_profile, self.private_settings_file) as pg_conn:
            # Check that this row exists
            results = get_all_results_for_query(pg_conn, query)
        if results:
            if nb_variant_extracted is not None:
                set_statements.append(f"num_ss_extracted = {nb_variant_extracted}")
            if nb_variant_remapped is not None:
                set_statements.append(f"num_ss_remapped = {nb_variant_remapped}")
            if nb_variant_ingested is not None:
                set_statements.append(f"num_ss_ingested = {nb_variant_ingested}")

        if set_statements:
            query = (
                f"UPDATE {self.tracking_table} SET {', '.join(set_statements)} "
                f"WHERE release_version={self.release_version} AND origin_assembly_accession='{source_assembly}' "
                f"AND taxonomy='{self.taxonomy}' AND source='{source}'"
            )
            with get_metadata_connection_handle(cfg['maven']['environment'], cfg['maven']['settings_file']) as pg_conn:
                execute_query(pg_conn, query)

    def count_variants_from_logs(self, assembly_directory, source_assembly):
        vcf_extractor_log = os.path.join(assembly_directory, 'logs', source_assembly + '_vcf_extractor.log')
        eva_remapping_count = os.path.join(assembly_directory, 'eva', source_assembly + '_eva_remapped_counts.yml')
        dbsnp_remapping_count = os.path.join(assembly_directory, 'dbsnp', source_assembly + '_dbsnp_remapped_counts.yml')
        eva_ingestion_log = os.path.join(assembly_directory, 'logs', source_assembly + '_eva_remapped.vcf_ingestion.log')
        dbsnp_ingestion_log = os.path.join(assembly_directory, 'logs', source_assembly + '_dbsnp_remapped.vcf_ingestion.log')

        eva_total, eva_written, dbsnp_total, dbsnp_written = count_variants_extracted(vcf_extractor_log)
        eva_candidate, eva_remapped, eva_unmapped = count_variants_remapped(eva_remapping_count)
        dbsnp_candidate, dbsnp_remapped, dbsnp_unmapped = count_variants_remapped(dbsnp_remapping_count)
        # Use the number of variant read rather than the number of variant ingested to get the total number of variant
        # when some might have been written in previous execution.
        eva_ingestion_candidate, eva_ingested, eva_duplicates = count_variants_ingested(eva_ingestion_log)
        dbsnp_ingestion_candidate, dbsnp_ingested, dbsnp_duplicates = count_variants_ingested(dbsnp_ingestion_log)

        self.set_counts(
            source_assembly, 'EVA',
            nb_variant_extracted=eva_written,
            nb_variant_remapped=eva_remapped,
            nb_variant_ingested=eva_ingestion_candidate
        )
        self.set_counts(
            source_assembly, 'DBSNP',
            nb_variant_extracted=dbsnp_written,
            nb_variant_remapped=dbsnp_remapped,
            nb_variant_ingested=dbsnp_ingestion_candidate
        )

        self.info(f'For Taxonomy: {self.taxonomy} and Assembly: {source_assembly} Source: EVA ')
        self.info(f'Number of variant read:{eva_total}, written:{eva_written}, attempt remapping: {eva_candidate}, '
                  f'remapped: {eva_remapped}, failed remapped {eva_unmapped}')
        self.info(f'For Taxonomy: {self.taxonomy} and Assembly: {source_assembly} Source: DBSNP ')
        self.info(
            f'Number of variant read:{dbsnp_total}, written:{dbsnp_written}, attempt remapping: {dbsnp_candidate}, '
            f'remapped: {dbsnp_remapped}, failed remapped {dbsnp_unmapped}')

    def update_dbs(self, source_of_assembly):
        """Update all relevant databases to reflect the new assembly."""
        incomplete_assemblies = self.get_incomplete_assemblies()
        if incomplete_assemblies:
            self.info(f'Processing for the following source assemblies is not yet complete: {incomplete_assemblies}')
            self.info('Not updating databases.')
            return
        self.add_to_supported_assemblies(source_of_assembly)
        self.add_to_metadata()
        self.add_to_contig_alias()

    def add_to_supported_assemblies(self, source_of_assembly):
        today = datetime.date.today().strftime('%Y-%m-%d')
        with get_metadata_connection_handle(self.maven_profile, self.private_settings_file) as pg_conn:
            # First check if the current assembly is already target - if so don't do anything
            current_query = (
                f"SELECT assembly_id FROM evapro.supported_assembly_tracker "
                f"WHERE taxonomy_id={self.taxonomy} AND current=true;"
            )
            results = get_all_results_for_query(pg_conn, current_query)
            if len(results) > 0 and results[0][0] == self.target_assembly:
                self.info(f'Current assembly for taxonomy {self.taxonomy} is already {self.target_assembly}')
                return

            # Deprecate the last current assembly
            update_query = (
                f"UPDATE evapro.supported_assembly_tracker "
                f"SET current=false, end_date='{today}' "
                f"WHERE taxonomy_id={self.taxonomy} AND current=true;"
            )
            execute_query(pg_conn, update_query)

            # Then insert the new assembly
            insert_query = (
                f"INSERT INTO evapro.supported_assembly_tracker "
                f"(taxonomy_id, source, assembly_id, current, start_date) "
                f"VALUES({self.taxonomy}, '{source_of_assembly}', '{self.target_assembly}', true, '{today}');"
            )
            execute_query(pg_conn, insert_query)

    def add_to_metadata(self):
        with get_metadata_connection_handle(self.maven_profile, self.private_settings_file) as pg_conn:
            insert_new_assembly_and_taxonomy(pg_conn, self.target_assembly, self.taxonomy)

    def add_to_contig_alias(self):
        contig_alias_url, contig_alias_user, contig_alias_pass = get_contig_alias_db_creds_for_profile(
            self.maven_profile, self.private_settings_file)
        client = ContigAliasClient(contig_alias_url, contig_alias_user, contig_alias_pass)
        client.insert_assembly(self.target_assembly)
