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
   --taxonomy_id 1234 \
   --source_assembly_accession GCA_0000001 \
	 --target_assembly_accession GCA_0000002 \
	 --species_name "Thingy thungus" \
	 --genome_assembly_dir ${SCRIPT_DIR}/genomes \
	 --extraction_properties ${SCRIPT_DIR}/template.properties \
	 --ingestion_properties ${SCRIPT_DIR}/template.properties \
	 --clustering_properties ${SCRIPT_DIR}/template.properties \
	 --clustering_instance 1 \
	 --output_dir ${SCRIPT_DIR}/output \
	 --remapping_config ${SCRIPT_DIR}/test_config.yaml \
	 --remapping_required 1

ls ${SCRIPT_DIR}/output/dbsnp/GCA_0000001_dbsnp_remapped.vcf \
   ${SCRIPT_DIR}/output/dbsnp/GCA_0000001_dbsnp_remapped_unmapped.vcf \
   ${SCRIPT_DIR}/output/dbsnp/GCA_0000001_dbsnp_remapped_counts.yml \
   ${SCRIPT_DIR}/output/eva/GCA_0000001_eva_remapped.vcf \
   ${SCRIPT_DIR}/output/eva/GCA_0000001_eva_remapped_unmapped.vcf \
   ${SCRIPT_DIR}/output/eva/GCA_0000001_eva_remapped_counts.yml

# Test we have 5 log files in the logs directory (1 extraction, 2 ingestion, 2 clustering)
[[ $(find ${SCRIPT_DIR}/output/logs/ -type f -name "*.log" | wc -l) -eq 5 ]]

# clean up
rm -rf work .nextflow* output genomes
cd ${cwd}
