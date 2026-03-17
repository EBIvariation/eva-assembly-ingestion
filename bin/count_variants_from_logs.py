#!/usr/bin/env python

# Copyright 2026 EMBL - European Bioinformatics Institute
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

from ebi_eva_common_pyutils.logger import logging_config

from eva_assembly_ingestion.config import load_config
from eva_assembly_ingestion.assembly_ingestion_job import AssemblyIngestionJob


def main():
    argparse = ArgumentParser(description='Gather counts from logs')
    argparse.add_argument('--taxonomy', required=True, type=int, help='Taxonomy id')
    argparse.add_argument('--source_assembly', required=True, type=str, help='Source assembly accession')
    argparse.add_argument('--target_assembly', required=True, type=str, help='Target assembly accession')
    argparse.add_argument('--output_directory', required=True, type=str, help='Path to processing directory')
    argparse.add_argument('--release_version', required=True, type=int, help='Release version')
    args = argparse.parse_args()

    load_config()

    job = AssemblyIngestionJob(args.taxonomy, args.target_assembly, args.release_version)
    logging_config.add_stdout_handler()

    job.count_variants_from_logs(args.output_directory, args.source_assembly, [args.taxonomy])


if __name__ == "__main__":
    main()
