# eva-assembly-ingestion
Automation for ingesting a new assembly for a species (taxonomy) in the EVA databases.

## Scripts

### Adding a target assembly
Primary job for ingesting a new assembly, including finding assemblies for the taxonomy that need to be remapped, doing
the remapping and clustering, and updating all related metadata and other databases.

Supports the following tasks:

* `load_tracker`: Retrieves source assemblies and number of studies and loads into the tracker.
  Will not load if any jobs exist for this taxonomy / target assembly pair.
* `remap_cluster`: Remaps all source assemblies in the tracker and clusters on the target assembly.
  Will only start or resume jobs not marked as `complete`.
* `update_dbs`: Updates the following with the new assembly: supported assembly table, metadata, and contig alias.
  Will not do any updates if any incomplete jobs are present in the tracker.

Example usage:
```bash
# Run everything
add_target_assembly.py --taxonomy 9031 --target_assembly GCA_016699485.1 --release_version 5

# Run remapping and clustering only, resume
add_target_assembly.py --taxonomy 9031 --target_assembly GCA_016699485.1 --release_version 5 --tasks remap_cluster --resume
```

### Custom assembly generation
Executable to generate custom assemblies and assembly reports.
This is called in the main target assembly job and can be used for other remapping jobs as well.
```bash
# Standard run
get_custom_assembly.py --assembly-accession GCA_016699485.1 --fasta-file /path/to/fasta --report-file /path/to/report

# Disable contig renaming
get_custom_assembly.py --assembly-accession GCA_016699485.1 --fasta-file /path/to/fasta --report-file /path/to/report --no-rename
```

### Genome target tracker

For every species in EVA metadata, check which assembly is currently supported by Ensembl and report if it matches with what is supported by EVA.
```bash
genome_target_tracker.py --private_config_xml_file /path/to/config.xml
```

## Configuration

The scripts require a configuration YAML file proving the locations of executables and other parameters.
The default config is `.assembly_config.yml` in the user's home, or a path can be provided via the environment variable `ASSEMBLYCONFIG`.

A complete config looks like the following:

```yaml
maven:
  environment: development
  settings_file: /path/to/settings.xml

remapping:
  base_directory: /path/to/remapping_dir

eutils_api_key: 12345

genome_downloader:
  output_directory: /path/to/genomes_dir

executable:
  python_activate: /path/to/remapping_env
  nextflow: /path/to/nextflow
  bcftools: /path/to/bcftools
  samtools: /path/to/samtools
  bedtools: /path/to/bedtools
  minimap2: /path/to/minimap
  bgzip: /path/to/bgzip
  tabix: /path/to/tabix
  genome_downloader: /path/to/genome_downloader
  custom_assembly: /path/to/custom_assembly

jar:
  vcf_extractor: /path/to/extraction.jar
  vcf_ingestion: /path/to/ingestion.jar
  clustering: /path/to/clustering.jar

nextflow:
  remapping: /path/to/remapping.nf
```
