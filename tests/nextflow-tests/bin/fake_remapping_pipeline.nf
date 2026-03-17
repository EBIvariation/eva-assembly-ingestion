#!/usr/bin/env nextflow

nextflow.enable.dsl=2

outfile_basename = file(params.outfile).getName()
outfile_no_extension = file(params.outfile).getBaseName()

workflow {
    remap_vcf()
}

process remap_vcf {

    publishDir workflow.launchDir, overwrite: true, mode: "copy", pattern: "*"

    output:
    path "${outfile_basename}", emit: remapped_vcf
    path "${outfile_no_extension}_counts.yml", emit: remapped_yml
    path "${outfile_no_extension}_unmapped.vcf", emit: unmapped_vcf

    script:
    """
    touch ${outfile_basename}
    touch ${outfile_no_extension}_counts.yml
    touch ${outfile_no_extension}_unmapped.vcf
    """
}
