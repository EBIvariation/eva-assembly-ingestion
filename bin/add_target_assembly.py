#!/usr/bin/env python

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
from argparse import ArgumentParser

from ebi_eva_common_pyutils.logger import logging_config

from eva_assembly_ingestion.config import load_config
from eva_assembly_ingestion.assembly_ingestion_job import AssemblyIngestionJob


def main():
    argparse = ArgumentParser(description='Add a new target assembly for a given taxonomy')
    argparse.add_argument('--taxonomy', required=True, type=int, help='Taxonomy id to be processed')
    argparse.add_argument('--target_assembly', required=True, type=str, help='New target assembly accession')
    argparse.add_argument('--source_of_assembly', required=False, type=str, default='Ensembl',
                          help='Source of new target assembly (default Ensembl)')
    argparse.add_argument('--tasks', required=False, type=str, nargs='+',
                          default=AssemblyIngestionJob.all_tasks, choices=AssemblyIngestionJob.all_tasks,
                          help='Task or set of tasks to perform (defaults to all)')
    argparse.add_argument('--release_version', required=True, type=int,
                          help='Release version this assembly will be processed for')
    argparse.add_argument('--resume', help='If a process has been run already this will resume it.',
                          action='store_true', default=False)
    args = argparse.parse_args()

    load_config()

    job = AssemblyIngestionJob(args.taxonomy, args.target_assembly, args.release_version)
    logging_config.add_stdout_handler()

    job.run_all(
        tasks=args.tasks,
        source_of_assembly=args.source_of_assembly,
        resume=args.resume
    )


if __name__ == "__main__":
    main()
