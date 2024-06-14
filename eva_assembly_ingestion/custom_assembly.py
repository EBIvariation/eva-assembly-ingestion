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
import os
import re
import shutil
import urllib
from copy import copy
from csv import DictReader, excel_tab, DictWriter
from typing import List, Dict

from cached_property import cached_property
from ebi_eva_common_pyutils.config import cfg
from ebi_eva_common_pyutils.logger import AppLogger
from ebi_eva_internal_pyutils.metadata_utils import get_metadata_connection_handle
from ebi_eva_internal_pyutils.pg_utils import get_all_results_for_query
from retry import retry


class CustomAssembly(AppLogger):
    """
    This class creates a custom assembly based on existing assembly fasta and report.
    It adds additional contigs provided by the required_contigs function after checking that they are not already in
    the assembly.
    It also renames all the sequence to INSDC accession.
    """
    def __init__(self, assembly_accession, assembly_fasta_path, assembly_report_path, no_rename=False,
                 eutils_api_key=None):
        self.assembly_accession = assembly_accession
        self.assembly_fasta_path = assembly_fasta_path
        self.assembly_report_path = assembly_report_path
        self.no_rename = no_rename
        self.eutils_api_key = eutils_api_key or cfg.get('eutils_api_key')
        self.assembly_report_headers = None

    @property
    def output_assembly_report_path(self):
        base, ext = os.path.splitext(self.assembly_report_path)
        return base + '_custom' + ext

    @property
    def output_assembly_fasta_path(self):
        base, ext = os.path.splitext(self.assembly_fasta_path)
        return base + '_custom' + ext

    @property
    def assembly_directory(self):
        return os.path.dirname(self.assembly_fasta_path)

    @cached_property
    def required_contigs(self) -> List[Dict]:
        """
        The contigs that are required to be added to the custom assembly. This provides a list of dict where each dict
        contains at least one key call 'genbank' and the value should be the INSDC accession of the sequence.
        The only other key supported is refseq.
        """
        raise NotImplementedError

    @staticmethod
    def _get_assembly_report(assembly_report):
        """Parse the assembly report and return each row as a dict."""
        headers = None
        with open(assembly_report) as open_file:
            # Parse the assembly report file to find the header then stop
            for line in open_file:
                if line.lower().startswith("# sequence-name") and "sequence-role" in line.lower():
                    headers = line.strip().split('\t')
                    break
            reader = DictReader(open_file, fieldnames=headers, dialect=excel_tab)
            return headers, [record for record in reader]

    @cached_property
    def assembly_report_rows(self):
        """Provides assembly report rows with each row as a dict."""
        headers, rows = self._get_assembly_report(self.assembly_report_path)
        self.assembly_report_headers = headers
        return rows

    @cached_property
    def extended_report_rows(self):
        """Provide the list of assembly report rows extended with additional ones if there are any."""
        if self.genbank_contig_to_add:
            extended_report_rows = copy(self.assembly_report_rows)
            for contig_dict in self.genbank_contig_to_add:
                row = {
                    "# Sequence-Name": contig_dict['genbank'],
                    "Sequence-Role": "scaffold",
                    "GenBank-Accn": contig_dict['genbank'],
                    'Relationship': '<>'
                }
                if 'refseq' in contig_dict:
                    row['RefSeq-Accn'] = contig_dict['refseq']
                    row['Relationship'] = '='
                extended_report_rows.append(row)
        else:
            extended_report_rows = self.assembly_report_rows
        return extended_report_rows

    @cached_property
    def genbank_contig_to_add(self):
        genbank_contigs = set(row['GenBank-Accn'] for row in self.assembly_report_rows)
        return [contig_dict for contig_dict in self.required_contigs if contig_dict['genbank'] not in genbank_contigs]

    @staticmethod
    def _get_contig_accessions_in_fasta(fasta_path):
        written_contigs = []
        match = re.compile(r'>(.*?)\s')
        if os.path.isfile(fasta_path):
            with open(fasta_path, 'r') as file:
                for line in file:
                    written_contigs.extend(match.findall(line))
        return written_contigs

    @staticmethod
    def rewrite_changing_names(input_fasta, output_fasta, contig_to_rename):
        with open(input_fasta) as open_input, open(output_fasta, 'w') as open_output:
            for line in open_input:
                if line.startswith('>'):
                    contig_name = line.split()[0][1:]
                    if contig_name in contig_to_rename:
                        line = '>' + contig_to_rename[contig_name] + '\n'
                open_output.write(line)

    @cached_property
    def contig_names_in_fasta(self):
        return self._get_contig_accessions_in_fasta(self.assembly_fasta_path)

    @cached_property
    def contig_to_rename(self):
        rename_map = {}
        genbank_contigs = set()
        map_to_genbank = {}
        for row in self.assembly_report_rows:
            if row['GenBank-Accn'] != 'na' and row['Relationship'] != '<>':
                genbank_contigs.add(row['GenBank-Accn'])
                map_to_genbank[row['RefSeq-Accn']] = row['GenBank-Accn']
                map_to_genbank[row['# Sequence-Name']] = row['GenBank-Accn']

        for name in self.contig_names_in_fasta:
            if name not in genbank_contigs:
                if name not in map_to_genbank:
                    self.error(f'Sequence {name} in fasta file does not match any INSDC sequence')
                else:
                    rename_map[name] = map_to_genbank[name]

        return rename_map

    @retry(tries=4, delay=2, backoff=1.2, jitter=(1, 3))
    def download_contig_from_ncbi(self, contig_accession):
        sequence_tmp_path = os.path.join(self.assembly_directory, contig_accession + '.fa')
        parameters = {
            'db': 'nuccore',
            'id': contig_accession,
            'rettype': 'fasta',
            'retmode': 'text',
            'tool': 'eva',
            'email': 'eva-dev@ebi.ac.uk'
        }
        if self.eutils_api_key:
            parameters['api_key'] = self.eutils_api_key
        url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' + urllib.parse.urlencode(parameters)
        self.info('Downloading ' + contig_accession)
        urllib.request.urlretrieve(url, sequence_tmp_path)
        return sequence_tmp_path

    def generate_assembly_report(self):
        if self.genbank_contig_to_add:
            self.info(f'Create custom assembly report for {self.assembly_accession}')
            with open(self.output_assembly_report_path, 'w') as open_output:
                writer = DictWriter(open_output, fieldnames=self.assembly_report_headers,  dialect=excel_tab, restval='na')
                writer.writeheader()
                for row in self.extended_report_rows:
                    writer.writerow(row)
        else:
            os.symlink(self.assembly_report_path, self.output_assembly_report_path)

    def generate_fasta(self):
        """
        Check if custom contig needs to be added to the assembly. If yes then copy the fasta file and append the new
        contigs otherwise create a symlink to the normal assembly.
        """
        # Find out what are the contigs that needs to be appended to the assembly
        contig_to_append = []
        for contig_dict in self.genbank_contig_to_add:
            if contig_dict['genbank'] not in self.contig_names_in_fasta:
                contig_to_append.append(self.download_contig_from_ncbi(contig_dict['genbank']))

        if contig_to_append or self.contig_to_rename:
            self.info(f'Create custom assembly fasta for {self.assembly_accession}')
            if self.contig_to_rename and not self.no_rename:
                self.rewrite_changing_names(self.assembly_fasta_path, self.output_assembly_fasta_path, self.contig_to_rename)
            else:
                shutil.copy(self.assembly_fasta_path, self.output_assembly_fasta_path, follow_symlinks=True)
            if contig_to_append:
                with open(self.output_assembly_fasta_path, 'a+') as fasta:
                    for contig_path in contig_to_append:
                        with open(contig_path) as sequence:
                            for line in sequence:
                                # Check that the line is not empty
                                if line.strip():
                                    fasta.write(line)
                        os.remove(contig_path)
        else:
            os.symlink(self.assembly_fasta_path, self.output_assembly_fasta_path)


class CustomAssemblyFromDatabase(CustomAssembly):

    @cached_property
    def required_contigs(self):
        """Return list of dict retrieve from the eva_tasks.eva2469_contig_analysis table."""
        self.info('Retrieve required contigs from database')
        with get_metadata_connection_handle(cfg['maven']['environment'], cfg['maven']['settings_file']) as pg_conn:
            query = ("select distinct contig_accession,refseq_contig_from_equiv_table from eva_tasks.eva2469_contig_analysis "
                     "where source_table in ('dbsnpSubmittedVariantEntity', 'submittedVariantEntity') "
                     f"and assembly_accession='{self.assembly_accession}'")
            return [dict([
                ('genbank', genbank_accession.strip() if genbank_accession else ''),
                ('refseq', refseq_accession.strip() if refseq_accession else '')
            ]) for genbank_accession, refseq_accession in get_all_results_for_query(pg_conn, query)]
