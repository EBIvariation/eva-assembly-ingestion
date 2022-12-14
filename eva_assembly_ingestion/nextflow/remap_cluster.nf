#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    Remap one assembly version to another, cluster, and QC.

    Inputs:
            --taxonomy_id                   taxonomy id of submitted variants that needs to be remapped.
            --source_assembly_accession     assembly accession of the submitted variants are currently mapped to.
            --target_assembly_accession     assembly accession the submitted variants will be remapped to.
            --species_name                  scientific name to be used for the species.
            --genome_assembly_dir           path to the directory where the genome should be downloaded.
            --extraction_properties         path to extraction properties file
            --ingestion_properties          path to ingestion properties file
            --clustering_properties         path to clustering properties file
            --clustering_instance           instance id to use for clustering
            --output_dir                    path to the directory where the output file should be copied.
            --remapping_config              path to the remapping configuration file
            --remapping_required            flag that sets the remapping as required if true otherwise the remapping is skipped and only the clustering can be run
            --memory                        memory in GB to use for memory-hungry processes (e.g. Java), default 8GB
    """
}

params.source_assembly_accession = null
params.target_assembly_accession = null
params.species_name = null
params.memory = 8
// help
params.help = null

// Show help message
if (params.help) exit 0, helpMessage()

// Test input files
if (!params.taxonomy_id || !params.source_assembly_accession || !params.target_assembly_accession || !params.species_name || !params.genome_assembly_dir ) {
    if (!params.taxonomy_id) log.warn('Provide the taxonomy id of the source submitted variants using --taxonomy_id')
    if (!params.source_assembly_accession) log.warn('Provide the source assembly using --source_assembly_accession')
    if (!params.target_assembly_accession) log.warn('Provide the target assembly using --target_assembly_accession')
    if (!params.species_name) log.warn('Provide a species name using --species_name')
    if (!params.genome_assembly_dir) log.warn('Provide a path to where the assembly should be downloaded using --genome_assembly_dir')
    exit 1, helpMessage()
}

species_name = params.species_name.toLowerCase().replace(" ", "_")
source_to_target = "${params.source_assembly_accession}_to_${params.target_assembly_accession}"

// Create an channel that will either be empty if remapping will take place or contain a dummy value if not
// This will allow to trigger the clustering even if no remapping is required
// We're using params.genome_assembly_dir because the clustering process needs to receive a file object
empty_ch = params.remapping_required ? Channel.empty() : Channel.of(params.genome_assembly_dir)

process retrieve_source_genome {

    output:
    path "${params.source_assembly_accession}.fa" into source_fasta
    path "${params.source_assembly_accession}_assembly_report.txt" into source_report

    """
    $params.executable.genome_downloader --assembly-accession ${params.source_assembly_accession} --species ${species_name} --output-directory ${params.genome_assembly_dir}
    ln -s ${params.genome_assembly_dir}/${species_name}/${params.source_assembly_accession}/${params.source_assembly_accession}.fa
    ln -s ${params.genome_assembly_dir}/${species_name}/${params.source_assembly_accession}/${params.source_assembly_accession}_assembly_report.txt
    """
}


process retrieve_target_genome {

    output:
    path "${params.target_assembly_accession}.fa" into target_fasta
    path "${params.target_assembly_accession}_assembly_report.txt" into target_report

    """
    $params.executable.genome_downloader --assembly-accession ${params.target_assembly_accession} --species ${species_name} --output-directory ${params.genome_assembly_dir}
    ln -s ${params.genome_assembly_dir}/${species_name}/${params.target_assembly_accession}/${params.target_assembly_accession}.fa
    ln -s ${params.genome_assembly_dir}/${species_name}/${params.target_assembly_accession}/${params.target_assembly_accession}_assembly_report.txt
    """
}

process update_source_genome {

    input:
    path source_fasta from source_fasta
    path source_report from source_report
    env REMAPPINGCONFIG from params.remapping_config

    output:
    path  "${source_fasta.getBaseName()}_custom.fa" into updated_source_fasta
    path  "${source_report.getBaseName()}_custom.txt" into updated_source_report

    """
    ${params.executable.custom_assembly} --assembly-accession ${params.source_assembly_accession} --fasta-file ${source_fasta} --report-file ${source_report}
    """
}

process update_target_genome {

    input:
    path target_fasta from target_fasta
    path target_report from target_report
    env REMAPPINGCONFIG from params.remapping_config

    output:
    path "${target_fasta.getBaseName()}_custom.fa" into updated_target_fasta
    path "${target_report.getBaseName()}_custom.txt" into updated_target_report

    """
    ${params.executable.custom_assembly} --assembly-accession ${params.target_assembly_accession} --fasta-file ${target_fasta} --report-file ${target_report} --no-rename
    """
}


/*
 * Extract the submitted variants to remap from the accesioning warehouse and store them in a VCF file.
 */
process extract_vcf_from_mongo {
    memory "${params.memory}GB"
    clusterOptions "-g /accession"

    when:
    params.remapping_required

    input:
    path source_fasta from updated_source_fasta
    path source_report from updated_source_report

    output:
    // Store both vcfs (eva and dbsnp) into one channel
    path '*.vcf' into source_vcfs
    path "${params.source_assembly_accession}_vcf_extractor.log" into log_filename

    publishDir "$params.output_dir/logs", overwrite: true, mode: "copy", pattern: "*.log*"

    """
    java -Xmx8G -jar $params.jar.vcf_extractor \
        --spring.config.name=${params.extraction_properties} \
        --parameters.fasta=${source_fasta} \
        --parameters.assemblyReportUrl=file:${source_report} \
        > ${params.source_assembly_accession}_vcf_extractor.log
    """
}


/*
 * Variant remapping pipeline
 */
process remap_variants {
    memory "${params.memory}GB"

    input:
    path source_fasta from updated_source_fasta
    path target_fasta from updated_target_fasta
    path source_vcf from source_vcfs.flatten()

    output:
    path "${basename_source_vcf}_remapped.vcf" into remapped_vcfs
    path "${basename_source_vcf}_remapped_unmapped.vcf" into unmapped_vcfs
    path "${basename_source_vcf}_remapped_counts.yml" into remapped_ymls

    publishDir "$params.output_dir/eva", overwrite: true, mode: "copy", pattern: "*_eva_remapped*"
    publishDir "$params.output_dir/dbsnp", overwrite: true, mode: "copy", pattern: "*_dbsnp_remapped*"

    script:
    basename_source_vcf = source_vcf.getBaseName()
    """
    # Setup the PATH so that the variant remapping pipeline can access its dependencies
    mkdir bin
    for P in $params.executable.bcftools $params.executable.samtools $params.executable.bedtools $params.executable.minimap2 $params.executable.bgzip $params.executable.tabix
      do ln -s \$P bin/
    done
    PATH=`pwd`/bin:\$PATH
    # Nextflow needs the full path to the input parameters hence the pwd
    $params.executable.nextflow run $params.nextflow.remapping -resume \
      --oldgenome `pwd`/${source_fasta} \
      --newgenome `pwd`/${target_fasta} \
      --vcffile `pwd`/${source_vcf} \
      --outfile `pwd`/${basename_source_vcf}_remapped.vcf
    """
}


/*
 * Ingest the remapped submitted variants from a VCF file into the accessioning warehouse.
 */
process ingest_vcf_into_mongo {
    memory "${params.memory}GB"
    clusterOptions "-g /accession"

    input:
    path remapped_vcf from remapped_vcfs.flatten()
    path target_report from updated_target_report

    output:
    path "${remapped_vcf}_ingestion.log" into ingestion_log_filename

    publishDir "$params.output_dir/logs", overwrite: true, mode: "copy", pattern: "*.log*"

    script:
    """
    # Check the file name to know which database to load the variants into
    if [[ $remapped_vcf == *_eva_remapped.vcf ]]
    then
        loadTo=EVA
    else
        loadTo=DBSNP
    fi

    java -Xmx8G -jar $params.jar.vcf_ingestion \
        --spring.config.name=${params.ingestion_properties} \
        --parameters.vcf=${remapped_vcf} \
        --parameters.assemblyReportUrl=file:${target_report} \
        --parameters.loadTo=\${loadTo} \
        > ${remapped_vcf}_ingestion.log
    """
}

process process_remapped_variants {
    memory "${params.memory}GB"
    clusterOptions "-g /accession"

    input:
    path ingestion_log from empty_ch.mix(ingestion_log_filename.collect())

    output:
    path "${source_to_target}_process_remapped.log" into process_remapped_log_filename
    // TODO this also generates a rs report, for "newly remapped" rs - should we QC this separately?
    // path "${source_to_target}_rs_report.txt" optional true into rs_report_filename

    publishDir "$params.output_dir/logs", overwrite: true, mode: "copy", pattern: "*.log*"

    """
    java -Xmx8G -jar $params.jar.clustering \
        --spring.config.name=${params.clustering_properties} \
        --spring.batch.job.names=PROCESS_REMAPPED_VARIANTS_WITH_RS_JOB \
        > ${source_to_target}_process_remapped.log
    """
}

process cluster_unclustered_variants {
    memory "${params.memory}GB"
    clusterOptions "-g /accession/instance-${params.clustering_instance}"

    input:
    path process_remapped_log from process_remapped_log_filename

    output:
    path "${source_to_target}_clustering.log" into clustering_log_filename
    path "${source_to_target}_rs_report.txt" optional true into rs_report_filename

    publishDir "$params.output_dir/logs", overwrite: true, mode: "copy", pattern: "*.log*"

    """
    java -Xmx8G -jar $params.jar.clustering \
        --spring.config.name=${params.clustering_properties} \
        --spring.batch.job.names=CLUSTER_UNCLUSTERED_VARIANTS_JOB \
        > ${source_to_target}_clustering.log
    """
}

/*
 * Run clustering QC job
 */
process qc_clustering {
    memory "${params.memory}GB"
    clusterOptions "-g /accession"

    input:
    path rs_report from rs_report_filename

    output:
    path "${source_to_target}_clustering_qc.log" into clustering_qc_log_filename

    publishDir "$params.output_dir/logs", overwrite: true, mode: "copy", pattern: "*.log*"

    """
    java -Xmx8G -jar $params.jar.clustering \
        --spring.config.name=${params.clustering_properties} \
        --spring.batch.job.names=NEW_CLUSTERED_VARIANTS_QC_JOB \
        > ${source_to_target}_clustering_qc.log
    """
}
