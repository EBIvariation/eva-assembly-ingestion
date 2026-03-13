#!/usr/bin/env nextflow


nextflow.enable.dsl=2

def helpMessage() {
    log.info"""
    Remap one assembly version to another, cluster, and QC.

    Inputs:
            --release_version                   release version
            --source_assemblies_and_taxonomies  source assemblies and taxonomies needing remapping
            --target_assembly_accession         assembly accession the submitted variants will be remapped to.
            --species_name                      scientific name to be used for the species.
            --genome_assembly_dir               path to the directory where the genome should be downloaded.
            --extraction_properties             path to extraction properties file
            --ingestion_properties              path to ingestion properties file
            --clustering_properties             path to clustering properties file
            --output_dir                        path to the directory where the output file should be copied.
            --remapping_config                  path to the remapping configuration file
    """
}

params.release_version = null
params.source_assemblies_and_taxonomies = null
params.target_assembly_accession = null
params.species_name = null
// help
params.help = null

// Show help message
if (params.help) exit 0, helpMessage()

// Test input files
if (!params.source_assemblies_and_taxonomies || !params.target_assembly_accession || !params.species_name || !params.genome_assembly_dir ) {
    if (!params.source_assemblies_and_taxonomies) log.warn('Provide the source assemblies and taxonomy ids of the source submitted variants using --source_assemblies_and_taxonomies')
    if (!params.target_assembly_accession) log.warn('Provide the target assembly using --target_assembly_accession')
    if (!params.species_name) log.warn('Provide a species name using --species_name')
    if (!params.genome_assembly_dir) log.warn('Provide a path to where the assembly should be downloaded using --genome_assembly_dir')
    exit 1, helpMessage()
}


workflow {
    species_name = params.species_name.toLowerCase().replace(" ", "_")
    remapping_required = params.source_assemblies_and_taxonomies.any { it[0] != params.target_assembly_accession }

    if (remapping_required) {
        // Get only source assemblies and tax IDs that require remapping
        source_asm_and_tax_ids = Channel.fromList(params.source_assemblies_and_taxonomies)
            .filter { (source_asm, tax_ids) -> source_asm != params.target_assembly_accession }

        // Process source genomes
        retrieve_source_genome(source_asm_and_tax_ids, species_name)
        update_source_genome(
            retrieve_source_genome.out.fasta_and_report,
            params.remapping_config)

        // Process target genome
        retrieve_target_genome(params.target_assembly_accession, species_name)
        update_target_genome(
            retrieve_target_genome.out.target_fasta,
            retrieve_target_genome.out.target_report,
            params.remapping_config)

        // Remap required source assemblies
        asm_tax_fasta_report = source_asm_and_tax_ids.combine(update_source_genome.out.updated_fasta_and_report, by: 0)
            .transpose()
        extract_vcf_from_mongo(asm_tax_fasta_report)
        remap_variants(
            extract_vcf_from_mongo.out.source_vcfs.transpose(),
            update_target_genome.out.updated_target_fasta)
        ingest_vcf_into_mongo(
            remap_variants.out.remapped_vcfs, 
            update_target_genome.out.updated_target_report)

        gather_counts(ingest_vcf_into_mongo.out.ingestion_log_filename)

        // Cluster target assembly
        process_remapped_variants(ingest_vcf_into_mongo.out.ingestion_log_filename.collect())
        qc_process_remapped(process_remapped_variants.out.rs_report_filename)
        cluster_unclustered_variants(qc_process_remapped.out.qc_log_filename)
        qc_clustering(cluster_unclustered_variants.out.rs_report_filename)
        qc_clustering_duplicate_rs_acc(cluster_unclustered_variants.out.rs_report_filename)
        
        // Backpropagate to source assemblies
        backpropagate_clusters(
            source_asm_and_tax_ids,
            qc_clustering.out.qc_log_filename,
            qc_clustering_duplicate_rs_acc.out.qc_log_filename)
    } else {
        // Only perform clustering on target assembly
        // We're using params.genome_assembly_dir because cluster_unclustered_variants needs to receive a file object
        cluster_unclustered_variants(params.genome_assembly_dir)
        qc_clustering(cluster_unclustered_variants.out.rs_report_filename)
        qc_clustering_duplicate_rs_acc(cluster_unclustered_variants.out.rs_report_filename)
    }
}


process retrieve_source_genome {
    label 'short_time', 'med_mem'

    input:
    tuple val(source_assembly_accession), val(taxonomy_list)
    val species_name

    output:
    tuple val(source_assembly_accession), path("${source_assembly_accession}.fa"), path("${source_assembly_accession}_assembly_report.txt"), emit: fasta_and_report

    script:
    """
    $params.executable.genome_downloader --assembly-accession ${source_assembly_accession} --species ${species_name} --output-directory ${params.genome_assembly_dir}
    ln -s ${params.genome_assembly_dir}/${species_name}/${source_assembly_accession}/${source_assembly_accession}.fa
    ln -s ${params.genome_assembly_dir}/${species_name}/${source_assembly_accession}/${source_assembly_accession}_assembly_report.txt
    """
}


process retrieve_target_genome {
    label 'short_time', 'med_mem'

    input:
    val target_assembly_accession
    val species_name

    output:
    path "${target_assembly_accession}.fa", emit: target_fasta
    path "${target_assembly_accession}_assembly_report.txt", emit: target_report

    script:
    """
    $params.executable.genome_downloader --assembly-accession ${params.target_assembly_accession} --species ${species_name} --output-directory ${params.genome_assembly_dir}
    ln -s ${params.genome_assembly_dir}/${species_name}/${target_assembly_accession}/${target_assembly_accession}.fa
    ln -s ${params.genome_assembly_dir}/${species_name}/${target_assembly_accession}/${target_assembly_accession}_assembly_report.txt
    """
}

process update_source_genome {
    label 'short_time', 'med_mem'

    input:
    tuple val(source_assembly_accession), path(source_fasta), path(source_report)
    env REMAPPINGCONFIG

    output:
    tuple val(source_assembly_accession), path("${source_fasta.getBaseName()}_custom.fa"), path("${source_report.getBaseName()}_custom.txt"), emit: updated_fasta_and_report

    script:
    """
    ${params.executable.custom_assembly} --assembly-accession ${source_assembly_accession} --fasta-file ${source_fasta} --report-file ${source_report}
    """
}

process update_target_genome {
    label 'short_time', 'med_mem'

    input:
    path target_fasta
    path target_report
    env REMAPPINGCONFIG

    output:
    path "${target_fasta.getBaseName()}_custom.fa", emit: updated_target_fasta
    path "${target_report.getBaseName()}_custom.txt", emit: updated_target_report

    script:
    """
    ${params.executable.custom_assembly} --assembly-accession ${params.target_assembly_accession} --fasta-file ${target_fasta} --report-file ${target_report} --no-rename
    """
}


/*
 * Extract the submitted variants to remap from the accessioning warehouse and store them in a VCF file.
 */
process extract_vcf_from_mongo {
    label 'long_time', 'med_mem'

    input:
    tuple val(source_assembly_accession), val(taxonomy), path(source_fasta), path(source_report)

    output:
    // Store both vcfs (eva and dbsnp), emit: one channel
    tuple val(source_assembly_accession), val(taxonomy), path(source_fasta), path(source_report), path('*.vcf'), emit: source_vcfs
    path "${source_assembly_accession}_${taxonomy}_vcf_extractor.log", emit: log_filename

    publishDir "$params.output_dir/logs", overwrite: true, mode: "copy", pattern: "*.log*"

    script:
    """
    java -Xmx${task.memory.toGiga()-1}G -jar $params.jar.vcf_extractor \
        --spring.config.location=file:${params.extraction_properties} \
        --parameters.assemblyAccession=${source_assembly_accession} \
        --parameters.fasta=${source_fasta} \
        --parameters.assemblyReportUrl=file:${source_report} \
        --parameters.taxonomy=${taxonomy} \
        > ${source_assembly_accession}_${taxonomy}_vcf_extractor.log
    """
}


/*
 * Variant remapping pipeline
 */
process remap_variants {
    label 'long_time', 'med_mem'

    input:
    tuple val(source_assembly_accession), val(taxonomy), path(source_fasta), path(source_report), path(source_vcf)
    path target_fasta

    output:
    tuple val(source_assembly_accession), val(taxonomy), path("${basename_source_vcf}_remapped.vcf"), emit: remapped_vcfs
    path "${basename_source_vcf}_remapped_unmapped.vcf", emit: unmapped_vcfs
    path "${basename_source_vcf}_remapped_counts.yml", emit: remapped_ymls

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
    source $params.executable.python_activate

    # Set nextflow config if needed - has to be passed via commandline arg rather than env var
    if [[ -z "\${ASSEMBLY_NEXTFLOW_CONFIG:-}" ]]
    then
        nextflow_config_flag=""
    else
        nextflow_config_flag="-c \${ASSEMBLY_NEXTFLOW_CONFIG}"
    fi

    # Nextflow needs the full path to the input parameters hence the pwd
    $params.executable.nextflow run $params.nextflow.remapping -resume \
      --oldgenome `pwd`/${source_fasta} \
      --newgenome `pwd`/${target_fasta} \
      --vcffile `pwd`/${source_vcf} \
      --outfile `pwd`/${basename_source_vcf}_remapped.vcf \
      \${nextflow_config_flag}
    """
}


/*
 * Ingest the remapped submitted variants from a VCF file into the accessioning warehouse.
 */
process ingest_vcf_into_mongo {
    label 'long_time', 'med_mem'

    // Run ingestions in serial to avoid race conditions when writing variants to Mongo.
    // Note this applies across source assemblies as well as across EVA/dbSNP from the same assembly.
    maxForks 1

    input:
    tuple val(source_assembly_accession), val(taxonomy), path(remapped_vcf)
    path target_report

    output:
    tuple val(source_assembly_accession), val(taxonomy), path("${remapped_vcf}_ingestion.log"), emit: ingestion_log_filename

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

    java -Xmx${task.memory.toGiga()-1}G -jar $params.jar.vcf_ingestion \
        --spring.config.location=file:${params.ingestion_properties} \
        --parameters.remappedFrom=${source_assembly_accession} \
        --parameters.vcf=${remapped_vcf} \
        --parameters.assemblyReportUrl=file:${target_report} \
        --parameters.loadTo=\${loadTo} \
        > ${remapped_vcf}_ingestion.log
    """
}

/*
 * Gather counts from remapping processing (extraction, remapping and ingestion)
 */
process gather_counts {
    label 'default_time', 'default_mem'

    input:
    tuple val(source_assembly_accession), val(taxonomy), path(ingestion_log)

    script:
    """
    ${params.executable.count_variants_from_logs} \
        --source_assembly ${source_assembly_accession} \
        --taxonomy ${taxonomy} \
        --target_assembly ${params.target_assembly_accession} \
        --output_directory ${params.output_dir} \
        --release_version ${params.release_version}
    """
}


process process_remapped_variants {
    label 'long_time', 'med_mem'

    input:
    val ingestion_output

    output:
    path "${params.target_assembly_accession}_process_remapped.log", emit: process_remapped_log_filename
    path "${params.target_assembly_accession}_remapped_rs_report.txt", emit: rs_report_filename

    publishDir "$params.output_dir/logs", overwrite: true, mode: "copy"

    script:
    """
    java -Xmx${task.memory.toGiga()-1}G -jar $params.jar.clustering \
        --spring.config.location=file:${params.clustering_properties} \
        --spring.batch.job.names=PROCESS_REMAPPED_VARIANTS_WITH_RS_JOB \
        > ${params.target_assembly_accession}_process_remapped.log

    # Ensure we always create an RS report, to trigger downstream processing
    touch ${params.target_assembly_accession}_rs_report.txt
    mv ${params.target_assembly_accession}_rs_report.txt ${params.target_assembly_accession}_remapped_rs_report.txt
    """
}

process cluster_unclustered_variants {
    label 'long_time', 'med_mem'

    input:
    path process_remapped_log

    output:
    path "${params.target_assembly_accession}_clustering.log", emit: clustering_log_filename
    path "${params.target_assembly_accession}_new_rs_report.txt", emit: rs_report_filename

    publishDir "$params.output_dir/logs", overwrite: true, mode: "copy"

    script:
    """
    java -Xmx${task.memory.toGiga()-1}G -jar $params.jar.clustering \
        --spring.config.location=file:${params.clustering_properties} \
        --spring.batch.job.names=CLUSTER_UNCLUSTERED_VARIANTS_JOB \
        > ${params.target_assembly_accession}_clustering.log

    # Ensure we always create an RS report, to trigger downstream processing
    touch ${params.target_assembly_accession}_rs_report.txt
    mv ${params.target_assembly_accession}_rs_report.txt ${params.target_assembly_accession}_new_rs_report.txt
    """
}

/*
 * Run QC for remapped RS IDs
 */
process qc_process_remapped {
    label 'long_time', 'med_mem'

    input:
    path rs_report

    output:
    path "${params.target_assembly_accession}_process_remapped_qc.log", emit: qc_log_filename

    publishDir "$params.output_dir/logs", overwrite: true, mode: "copy", pattern: "*.log*"

    script:
    """
    java -Xmx${task.memory.toGiga()-1}G -jar $params.jar.clustering \
        --spring.config.location=file:${params.clustering_properties} \
        --spring.batch.job.names=NEW_CLUSTERED_VARIANTS_QC_JOB \
        > ${params.target_assembly_accession}_process_remapped_qc.log
    """
}

/*
 * Run QC for new RS IDs
 */
process qc_clustering {
    label 'long_time', 'med_mem'

    input:
    path rs_report

    output:
    path "${params.target_assembly_accession}_clustering_qc.log", emit: qc_log_filename

    publishDir "$params.output_dir/logs", overwrite: true, mode: "copy", pattern: "*.log*"

    script:
    """
    java -Xmx${task.memory.toGiga()-1}G -jar $params.jar.clustering \
        --spring.config.location=file:${params.clustering_properties} \
        --spring.batch.job.names=NEW_CLUSTERED_VARIANTS_QC_JOB \
        > ${params.target_assembly_accession}_clustering_qc.log
    """
}

/*
 * Run QC job for detecting duplicate RS accessions
 */
process qc_clustering_duplicate_rs_acc {
    label 'long_time', 'med_mem'

    input:
    path rs_report

    output:
    path "${params.target_assembly_accession}_clustering_qc_duplicate_rs_acc.log", emit: qc_log_filename

    publishDir "$params.output_dir/logs", overwrite: true, mode: "copy", pattern: "*.log*"

    script:
    """
    java -Xmx${task.memory.toGiga()-1}G -jar $params.jar.clustering \
         --spring.config.location=file:${params.clustering_properties} \
         --spring.batch.job.names=DUPLICATE_RS_ACC_QC_JOB \
         --parameters.duplicateRSAccFile=${params.target_assembly_accession}_duplicate_rs_accessions.txt \
         > ${params.target_assembly_accession}_clustering_qc_duplicate_rs_acc.log

    # Fail if the file is not empty
    if [ -s ${params.target_assembly_accession}_duplicate_rs_accessions.txt ]; then
        echo "Duplicate RS accessions detected! Failing the process."
        exit 1
    fi
    """
}


/*
 * Run Back propagation of new clustered RS
 */
process backpropagate_clusters {
    label 'long_time', 'med_mem'

    input:
    tuple val(source_assembly_accession), val(taxonomy)
    path clustering_qc_log
    path clustering_duplicate_qc_log

    output:
    path "${params.target_assembly_accession}_backpropagate_to_${source_assembly_accession}.log", emit: backpropagate_log_filename

    publishDir "$params.output_dir/logs", overwrite: true, mode: "copy", pattern: "*.log*"

    script:
    """
    java -Xmx${task.memory.toGiga()-1}G -jar $params.jar.clustering \
        --spring.config.location=file:${params.clustering_properties} \
        --parameters.remappedFrom=${source_assembly_accession} \
        --spring.batch.job.names=BACK_PROPAGATE_SPLIT_OR_MERGED_RS_JOB \
        > ${params.target_assembly_accession}_backpropagate_to_${source_assembly_accession}.log
    """
}
