import os

from eva_assembly_ingestion.parse_counts import count_variants_extracted, count_variants_ingested


def test_count_variants_remapped():
    pass


def test_count_variants_extracted():
    log_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'resources', 'vcf_extractor.log'))
    assert count_variants_extracted(log_file) == (25434599, 25434599, 72446279, 72446277)


def test_count_variants_ingested():
    log_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'resources', 'vcf_ingestion.log'))
    assert count_variants_ingested(log_file) == (1, 1, 0)
