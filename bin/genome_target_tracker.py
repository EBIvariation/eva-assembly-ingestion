#!/usr/bin/env python

# Copyright 2023 EMBL - European Bioinformatics Institute
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
from argparse import ArgumentParser

import requests
from ebi_eva_common_pyutils.logger import logging_config
from ebi_eva_common_pyutils.metadata_utils import get_metadata_connection_handle, insert_new_assembly_and_taxonomy
from ebi_eva_common_pyutils.pg_utils import get_all_results_for_query
from ebi_eva_common_pyutils.taxonomy.taxonomy import get_scientific_name_from_ensembl

logger = logging_config.get_logger(__name__)
logging_config.add_stdout_handler()

remapping_genome_target_table = 'evapro.supported_assembly_tracker'
ensembl_url = 'http://rest.ensembl.org/info/assembly'


def get_all_taxonomies_from_eva(private_config_xml_file):
    taxonomy_list = []
    with get_metadata_connection_handle("production_processing", private_config_xml_file) as pg_conn:
        query = 'SELECT DISTINCT taxonomy_id FROM evapro.taxonomy'
        for taxonomy in get_all_results_for_query(pg_conn, query):
            taxonomy_list.append(taxonomy[0])

    return taxonomy_list


def get_tax_latest_asm_from_eva(private_config_xml_file):
    eva_tax_asm = {}
    with get_metadata_connection_handle("production_processing", private_config_xml_file) as pg_conn:
        query = f"""SELECT DISTINCT taxonomy_id, source, assembly_id FROM {remapping_genome_target_table} 
        WHERE current=TRUE"""
        for tax_id, source, assembly in get_all_results_for_query(pg_conn, query):
            eva_tax_asm[tax_id] = {'assembly': assembly, 'source': source}

    return eva_tax_asm


def add_assembly_to_accessioned_assemblies(private_config_xml_file, taxonomies_to_assemblies):
    with get_metadata_connection_handle("production_processing", private_config_xml_file) as pg_conn:
        for taxonomy in taxonomies_to_assemblies:
            assembly = taxonomies_to_assemblies[taxonomy]['assembly']
            insert_new_assembly_and_taxonomy(pg_conn, assembly, taxonomy)


def get_tax_asm_from_sources(eva_tax_asm_source):
    source_tax_asm = {}
    for tax_id in eva_tax_asm_source:
        # Check for each taxonomy which source is present in table and try to get supported assembly from that source.
        # Currently only Ensembl is supported
        source_in_eva = eva_tax_asm_source[tax_id]["source"]
        if source_in_eva == 'Ensembl':
            tax_asm_from_ensembl = get_tax_asm_from_ensembl(tax_id)
            if tax_asm_from_ensembl is not None:
                source_tax_asm[tax_id] = tax_asm_from_ensembl
        else:
            logger.error(
                f'No implementation present to check assembly supported by Source ({source_in_eva}) for taxonomy {tax_id}')

    return source_tax_asm


def get_tax_asm_from_ensembl(tax_id):
    try:
        logger.info(f'Query Ensembl for species name using taxonomy {tax_id}')
        sp_name = get_scientific_name_from_ensembl(tax_id)
    except Exception:
        logger.warning(f'Could not get species name for taxonomy {tax_id} in Ensembl')
        return None
    # Get assembly from Ensembl
    logger.info(f'Query Ensembl for supported assembly for taxonomy {tax_id} using species_name "{sp_name}"')
    assembly = get_supported_asm_from_ensembl(sp_name)
    if assembly != 'None':
        return {'assembly': assembly, 'source': 'Ensembl'}
    else:
        logger.warning(
            f'Could not find supported assembly for taxonomy_id {tax_id} using species_name "{sp_name}" in Ensembl')
        return None


def get_supported_asm_from_ensembl(scientific_name):
    url = ensembl_url + '/' + scientific_name.lower().replace(' ', '_')
    response = requests.get(url, params={'content-type': 'application/json'})
    data = response.json()
    assembly_accession = str(data.get('assembly_accession'))
    return assembly_accession


def check_supported_target_assembly(private_config_xml_file):
    taxonomy_list = get_all_taxonomies_from_eva(private_config_xml_file)
    eva_tax_asm_source = get_tax_latest_asm_from_eva(private_config_xml_file)
    source_tax_asm = get_tax_asm_from_sources(eva_tax_asm_source)

    taxonomy_with_mismatch_assembly = {}
    taxonomy_not_tracked_by_eva = []
    taxonomy_tracked_but_not_retrieved_from_source = []

    for tax_id in taxonomy_list:
        if tax_id in eva_tax_asm_source and tax_id in source_tax_asm:
            if eva_tax_asm_source[tax_id]["assembly"] != source_tax_asm[tax_id]["assembly"]:
                taxonomy_with_mismatch_assembly[tax_id] = {
                    "eva": eva_tax_asm_source[tax_id]["assembly"],
                    "source": source_tax_asm[tax_id]["assembly"]
                }
        else:
            if tax_id not in eva_tax_asm_source:
                taxonomy_not_tracked_by_eva.append(tax_id)
            elif tax_id in eva_tax_asm_source and tax_id not in source_tax_asm:
                taxonomy_tracked_but_not_retrieved_from_source.append(tax_id)

    print_summary(taxonomy_with_mismatch_assembly, taxonomy_not_tracked_by_eva,
                  taxonomy_tracked_but_not_retrieved_from_source)

    return taxonomy_with_mismatch_assembly, taxonomy_not_tracked_by_eva, taxonomy_tracked_but_not_retrieved_from_source


def print_summary(taxonomy_with_mismatch_assembly, taxonomy_not_tracked_by_eva,
                  taxonomy_tracked_but_not_retrieved_from_source):
    if taxonomy_not_tracked_by_eva:
        logger.warning(f'The following taxonomy are not tracked by EVA: {taxonomy_not_tracked_by_eva}')

    if taxonomy_tracked_but_not_retrieved_from_source:
        logger.error(
            f'The following taxonomy are tracked by EVA but their corresponding assemblies could not be retrieved from sources for checking: '
            f'{taxonomy_tracked_but_not_retrieved_from_source}')

    if taxonomy_with_mismatch_assembly:
        logger.error(f'Taxonomy with different supported assemblies in EVA and Source:')
        for tax_id, diff_asm in taxonomy_with_mismatch_assembly.items():
            logger.error(
                f"Taxonomy {tax_id} has different supported assemblies in EVA({diff_asm['eva']}) "
                f"and Source({diff_asm['source']})")


def main():
    argparse = ArgumentParser(description='Track the currently supported assembly by Ensembl for every species '
                                          'and report if its matches with what is supported by EVA')
    argparse.add_argument('--private_config_xml_file', required=True,
                          help='Path to the file containing the username/passwords to access '
                               'production and development databases')
    args = argparse.parse_args()

    check_supported_target_assembly(args.private_config_xml_file)


if __name__ == "__main__":
    main()
