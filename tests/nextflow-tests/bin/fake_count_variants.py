#!/usr/bin/env python
import argparse


if __name__ == "__main__":
    argparse = argparse.ArgumentParser(description='Gather counts from logs')
    argparse.add_argument('--taxonomy', required=True)
    argparse.add_argument('--source_assembly', required=True)
    argparse.add_argument('--target_assembly', required=True)
    argparse.add_argument('--output_directory', required=True)
    argparse.add_argument('--release_version', required=True)
    args = argparse.parse_args()
