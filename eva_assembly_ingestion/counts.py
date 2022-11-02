#!/usr/bin/env python
import re

import yaml
from ebi_eva_common_pyutils.command_utils import run_command_with_output


# TODO this is really only for updating tracking table
def count_variants_remapped(count_yml_file):
    with open(count_yml_file) as open_file:
        data = yaml.safe_load(open_file)
    candidate_variants = data.get('all')
    remapped_variants = data.get('Flank_50', {}).get('Remapped', 0) + \
                        data.get('Flank_2000', {}).get('Remapped', 0) + \
                        data.get('Flank_50000', {}).get('Remapped', 0)
    unmapped_variants = data.get('filtered') or 0 + \
                        (data.get('Flank_50000', {}).get('total', 0) - data.get('Flank_50000', {}).get('Remapped', 0))
    return candidate_variants, remapped_variants, unmapped_variants


def parse_log_line(line, regex_list=None):
    if not regex_list:
        regex_list = [r'Items read = (\d+)', r'items written = (\d+)']
    results = []
    for regex in regex_list:
        match = re.search(regex, line)
        if match:
            results.append(int(match.group(1)))
        else:
            results.append(None)
    return tuple(results)


def count_variants_extracted(extraction_log):
    command = f'grep "EXPORT_EVA_SUBMITTED_VARIANTS_STEP" {extraction_log} | tail -1'
    log_line = run_command_with_output('Get total number of eva variants written', command, return_process_output=True)
    eva_total, eva_written = parse_log_line(log_line)
    command = f'grep "EXPORT_DBSNP_SUBMITTED_VARIANTS_STEP" {extraction_log} | tail -1'
    log_line = run_command_with_output('Get total number of dbsnp variants written', command, return_process_output=True)
    dbsnp_total, dbnp_written = parse_log_line(log_line)
    return eva_total, eva_written, dbsnp_total, dbnp_written


def count_variants_ingested(ingestion_log):
    command = f'grep "INGEST_REMAPPED_VARIANTS_FROM_VCF_STEP" {ingestion_log} | tail -1'
    log_line = run_command_with_output('Get total number of variants written', command, return_process_output=True)
    regex_list = [r'Items \(remapped ss\) read = (\d+)', r'ss ingested = (\d+)', r'ss skipped \(duplicate\) = (\d+)']
    ss_read, ss_written, ss_duplicates = parse_log_line(log_line, regex_list)
    return ss_read, ss_written, ss_duplicates
