#!/bin/bash
# Test the case where all source assemblies are the same as the target assembly.

set -Eeuo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
SOURCE_DIR="$(dirname $(dirname $SCRIPT_DIR))"

cwd=${PWD}
cd ${SCRIPT_DIR}

mkdir -p ${SCRIPT_DIR}/genomes
PATH=${SCRIPT_DIR}/bin:$PATH

printf "\e[32m===== CLUSTERING ONLY PIPELINE (source == target assembly) =====\e[0m\n"
nextflow run ${SOURCE_DIR}/eva_assembly_ingestion/nextflow/remap_cluster.nf -params-file test_config_no_remapping.yaml \
     --target_assembly_accession GCA_0000002 \
     --species_name "Thingy thungus" \
     --genome_assembly_dir ${SCRIPT_DIR}/genomes \
     --extraction_properties ${SCRIPT_DIR}/template.properties \
     --ingestion_properties ${SCRIPT_DIR}/template.properties \
     --clustering_properties ${SCRIPT_DIR}/template.properties \
     --output_dir ${SCRIPT_DIR}/output \
     --remapping_config ${SCRIPT_DIR}/test_config_no_remapping.yaml \
     --release_version 7 \
     -resume


# Verify clustering RS report was produced
ls ${SCRIPT_DIR}/output/logs/GCA_0000002_new_rs_report.txt

# clean up
rm -rf work .nextflow* output genomes
cd ${cwd}
