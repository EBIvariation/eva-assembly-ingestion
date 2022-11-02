import datetime

from cached_property import cached_property
from ebi_eva_common_pyutils.config import cfg
from ebi_eva_common_pyutils.config_utils import get_contig_alias_db_creds_for_profile
from ebi_eva_common_pyutils.contig_alias.contig_alias import ContigAliasClient
from ebi_eva_common_pyutils.logger import AppLogger
from ebi_eva_common_pyutils.metadata_utils import get_metadata_connection_handle, insert_new_assembly_and_taxonomy
from ebi_eva_common_pyutils.pg_utils import execute_query, get_all_results_for_query
from ebi_eva_common_pyutils.taxonomy.taxonomy import get_scientific_name_from_taxonomy
from psycopg2.extras import execute_values


class TaxonomyRemappingJob(AppLogger):
    """
    load_tracker:
        - find source assemblies
        - get number of studies and variants (?)
        - update tracking table
    remap_cluster:
        - remap all source assemblies in tracker
        - cluster target assembly
        - update tracker
    update_dbs:
        - update supported assemblies
        - update metadata (used by webservices)
        - update contig alias
    """
    all_tasks = ['load_tracker', 'remap_cluster', 'update_dbs']
    tracking_table = 'eva_progress_tracker.remapping_tracker'

    def __init__(self, taxonomy, target_assembly, source_of_assembly):
        self.taxonomy = taxonomy
        self.target_assembly = target_assembly
        self.source_of_assembly = source_of_assembly
        self.private_settings_file = cfg['maven']['settings_file']
        self.maven_profile = cfg['maven']['environment']

    @cached_property
    def scientific_name(self):
        return get_scientific_name_from_taxonomy(self.taxonomy)

    def run_all(self, tasks=None, resume=False):
        if not tasks:
            tasks = self.all_tasks
        if 'load_tracker' in tasks:
            self.load_tracker()
        if 'remap_cluster' in tasks:
            self.run_remapping_and_clustering()
        if 'update_dbs' in tasks:
            self.update_dbs()

    def load_tracker(self):
        """Load the tracking table with the source assemblies for this taxonomy"""
        # TODO check whether there's anything in the tracker already for target_assembly
        # TODO which of these should be removed? which should be set?
        column_names = ('source', 'taxonomy', 'scientific_name', 'origin_assembly_accession', 'assembly_accession',
                        'remapping_version', 'release_version', 'num_studies', 'num_ss_ids', 'study_accessions')
        rows = []
        source_assemblies = []
        for source_assembly, projects in self.get_source_assemblies_and_projects():
            source_assemblies.append(source_assembly)
            rows.append(('EVA, DBSNP', self.taxonomy, self.scientific_name, source_assembly, self.target_assembly,
                         None, None, len(projects), 1, None))
        if len(rows) == 0:
            self.info(f'Nothing to remap for taxonomy {self.taxonomy} and target assembly {self.target_assembly}')
            return
        with get_metadata_connection_handle(self.maven_profile, self.private_settings_file) as pg_conn:
            with pg_conn.cursor() as cursor:
                insert_query = (
                    f"INSERT INTO {self.tracking_table} ({','.join(column_names)}) VALUES %s "
                    f"ON CONFLICT () "
                    f"DO UPDATE SET "
                )
                # TODO fill in the conflict clause
                execute_values(cursor, insert_query, rows)
        self.info(f'Added the following assemblies to {self.tracking_table}: {source_assemblies}')

    def get_source_assemblies_and_projects(self):
        """Query metadata for all public projects with this taxonomy, of these getting all reference accessions
        for all analyses."""
        # TODO get dbSNP assemblies as well
        # TODO what if there are multiple analyses with different reference assemblies
        with get_metadata_connection_handle(self.maven_profile, self.private_settings_file) as pg_conn:
            query = (f"SELECT DISTINCT vcf_reference_accession, ARRAY_AGG(project_accession) "
                     f"FROM evapro.project "
                     f"LEFT OUTER JOIN evapro.project_taxonomy USING (project_accession) "
                     f"LEFT OUTER JOIN evapro.project_analysis USING (project_accession) "
                     f"LEFT OUTER JOIN evapro.analysis USING (analysis_accession) "
                     f"WHERE taxonomy_id={self.taxonomy} AND ena_status=4 AND hidden_in_eva=0 "
                     f"GROUP BY vcf_reference_accession"
                     )
            return get_all_results_for_query(pg_conn, query)

    def run_remapping_and_clustering(self):
        # TODO get sources from tracking table and only run if not complete to ensure idempotent
        pass

    def update_dbs(self):
        """Update all relevant databases to reflect the new assembly."""
        # TODO check tracking table if status is complete before updating
        self.add_to_supported_assemblies()
        self.add_to_metadata()
        self.add_to_contig_alias()

    def add_to_supported_assemblies(self):
        # TODO not idempotent - if target_assembly is already the current don't do anything!
        today = datetime.date.today().strftime('%Y-%m-%d')
        with get_metadata_connection_handle(self.maven_profile, self.private_settings_file) as pg_conn:
            # First deprecate the last current assembly
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
                f"VALUES({self.taxonomy}, '{self.source_of_assembly}', '{self.target_assembly}', true, '{today}');"
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
