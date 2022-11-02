#!/usr/bin/env python
import argparse

from ebi_eva_common_pyutils.logger import logging_config

from eva_assembly_ingestion.config import load_config
from eva_assembly_ingestion.custom_assembly import CustomAssemblyFromDatabase


def main():
    parser = argparse.ArgumentParser(description='Generate custom assembly report for a given assembly',
                                     add_help=False)
    parser.add_argument("-a", "--assembly-accession",
                        help="Assembly for which the process has to be run, e.g. GCA_000002315.3",
                        required=True)
    parser.add_argument("-f", "--fasta-file", help="Path to the fasta file containing the assembly", required=True)
    parser.add_argument("-r", "--report-file",
                        help="Path to the assembly report file containing the assembly", required=True)
    parser.add_argument("--no-rename", help="Disable renaming of contigs", default=False, action='store_true')
    parser.add_argument('--help', action='help', help='Show this help message and exit')

    args = parser.parse_args()

    load_config()
    logging_config.add_stdout_handler()

    assembly = CustomAssemblyFromDatabase(args.assembly_accession, args.fasta_file, args.report_file, args.no_rename)
    assembly.generate_assembly_report()
    assembly.generate_fasta()


if __name__ == "__main__":
    main()
