#!/usr/bin/env python
from argparse import ArgumentParser, ArgumentError

from ebi_eva_common_pyutils.logger import logging_config

from eva_assembly_ingestion.config import load_config
from eva_assembly_ingestion.job import TaxonomyRemappingJob


def main():
    argparse = ArgumentParser(description='Add a new target assembly for a given taxonomy')
    argparse.add_argument('--taxonomy', help='Taxonomy id to be processed')
    argparse.add_argument('--target_assembly', help='New target assembly accession')
    argparse.add_argument('--source_of_assembly', default='Ensembl',
                          help='Source of new target assembly (default Ensembl)')
    argparse.add_argument('--instance', help="Accessioning instance id for clustering", required=False, default=6,
                          type=int, choices=range(1, 13))
    argparse.add_argument('--tasks', required=False, type=str, nargs='+',
                          default=TaxonomyRemappingJob.all_tasks, choices=TaxonomyRemappingJob.all_tasks,
                          help='Task or set of tasks to perform (defaults to all)')
    argparse.add_argument('--resume', help='If a process has been run already this will resume it.',
                          action='store_true', default=False)
    args = argparse.parse_args()

    load_config()

    if not args.taxonomy or not args.target_assembly:
        raise ArgumentError(None, 'Must provide both --taxonomy and --target_assembly')

    job = TaxonomyRemappingJob(args.taxonomy, args.target_assembly, args.source_of_assembly)
    logging_config.add_stdout_handler()

    job.run_all(args.tasks)


if __name__ == "__main__":
    main()
