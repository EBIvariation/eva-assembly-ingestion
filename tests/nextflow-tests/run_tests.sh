#!/bin/bash

set -Eeuo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
SOURCE_DIR="$(dirname $(dirname $SCRIPT_DIR))"

cwd=${PWD}
cd ${SCRIPT_DIR}

mkdir -p ${SCRIPT_DIR}/genomes
PATH=${SCRIPT_DIR}/bin:$PATH

printf "\e[32m===== REMAPPING AND CLUSTERING PIPELINE =====\e[0m\n"
nextflow run ${SOURCE_DIR}/eva_assembly_ingestion/nextflow/remap_cluster.nf -params-file test_config.yaml \
	 --target_assembly_accession GCA_0000002 \
	 --species_name "Thingy thungus" \
	 --genome_assembly_dir ${SCRIPT_DIR}/genomes \
	 --extraction_properties ${SCRIPT_DIR}/template.properties \
	 --ingestion_properties ${SCRIPT_DIR}/template.properties \
	 --clustering_properties ${SCRIPT_DIR}/template.properties \
	 --output_dir ${SCRIPT_DIR}/output \
	 --remapping_config ${SCRIPT_DIR}/test_config.yaml \
	 --release_version 7 \
	 -resume

ls ${SCRIPT_DIR}/output/dbsnp/GCA_0000001.1_1233_dbsnp_remapped.vcf \
   ${SCRIPT_DIR}/output/dbsnp/GCA_0000001.1_1233_dbsnp_remapped_unmapped.vcf \
   ${SCRIPT_DIR}/output/dbsnp/GCA_0000001.1_1233_dbsnp_remapped_counts.yml \
   ${SCRIPT_DIR}/output/eva/GCA_0000001.1_1233_eva_remapped.vcf \
   ${SCRIPT_DIR}/output/eva/GCA_0000001.1_1233_eva_remapped_unmapped.vcf \
   ${SCRIPT_DIR}/output/eva/GCA_0000001.1_1233_eva_remapped_counts.yml \
   ${SCRIPT_DIR}/output/dbsnp/GCA_0000001.1_1234_dbsnp_remapped.vcf \
   ${SCRIPT_DIR}/output/dbsnp/GCA_0000001.1_1234_dbsnp_remapped_unmapped.vcf \
   ${SCRIPT_DIR}/output/dbsnp/GCA_0000001.1_1234_dbsnp_remapped_counts.yml \
   ${SCRIPT_DIR}/output/eva/GCA_0000001.1_1234_eva_remapped.vcf \
   ${SCRIPT_DIR}/output/eva/GCA_0000001.1_1234_eva_remapped_unmapped.vcf \
   ${SCRIPT_DIR}/output/eva/GCA_0000001.1_1234_eva_remapped_counts.yml \
   ${SCRIPT_DIR}/output/dbsnp/GCA_0000001.2_1234_dbsnp_remapped.vcf \
   ${SCRIPT_DIR}/output/dbsnp/GCA_0000001.2_1234_dbsnp_remapped_unmapped.vcf \
   ${SCRIPT_DIR}/output/dbsnp/GCA_0000001.2_1234_dbsnp_remapped_counts.yml \
   ${SCRIPT_DIR}/output/eva/GCA_0000001.2_1234_eva_remapped.vcf \
   ${SCRIPT_DIR}/output/eva/GCA_0000001.2_1234_eva_remapped_unmapped.vcf \
   ${SCRIPT_DIR}/output/eva/GCA_0000001.2_1234_eva_remapped_counts.yml

# Test we have 16 log files in the logs directory:
# 3 extraction, 6 ingestion, 1 process remapped, 1 clustering, 3 QC, 2 backpropagate
[[ $(find ${SCRIPT_DIR}/output/logs/ -type f -name "*.log" | wc -l) -eq 16 ]]

# Test we have 2 rs_reports in the logs directory
[[ $(find ${SCRIPT_DIR}/output/logs/ -type f -name "*.txt" | wc -l) -eq 2 ]]

# clean up
rm -rf work .nextflow* output genomes
cd ${cwd}
