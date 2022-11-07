#!/usr/bin/env python
import argparse


def touch(f):
    open(f, 'w').close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--assembly-accession", required=True)
    parser.add_argument("-f", "--fasta-file", required=True)
    parser.add_argument("-r", "--report-file", required=True)
    parser.add_argument("--no-rename", action='store_true')
    args = parser.parse_args()
    touch(args.fasta_file.replace('.fa', '_custom.fa'))
    touch(args.report_file.replace('.txt', '_custom.txt'))
