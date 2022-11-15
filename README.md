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

# Run remapping and clustering only, resume and run on a specific instance
add_target_assembly.py --taxonomy 9031 --target_assembly GCA_016699485.1 --release_version 5 --tasks remap_cluster --instance 3 --resume
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
