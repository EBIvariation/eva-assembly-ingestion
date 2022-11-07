#!/usr/bin/env nextflow


outfile_basename = file(params.outfile).getName()
outfile_no_extension = file(params.outfile).getBaseName()

process remap_vcf {

    publishDir workflow.launchDir

    output:
    path "${outfile_basename}" into remapped_vcf
    path "${outfile_no_extension}_counts.yml" into remapped_yml
    path "${outfile_no_extension}_unmapped.vcf" into unmapped_vcf

    script:
    """
    touch ${outfile_basename}
    touch ${outfile_no_extension}_counts.yml
    touch ${outfile_no_extension}_unmapped.vcf
    """
}
