#!/usr/bin/env python
import argparse
import os


def touch(f):
    open(f, 'w').close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--assembly-accession", required=True)
    parser.add_argument("-s", "--species", required=True)
    parser.add_argument("-o", "--output-directory")
    args = parser.parse_args()
    d = os.path.join(args.output_directory, args.species.lower().replace(' ', '_'), args.assembly_accession)
    os.makedirs(d, exist_ok=True)
    touch(os.path.join(d, args.assembly_accession + '.fa'))
    touch(os.path.join(d, args.assembly_accession + '_assembly_report.txt'))
