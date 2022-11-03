import os

from eva_assembly_ingestion.parse_counts import count_variants_extracted, count_variants_ingested, count_variants_remapped


def test_count_variants_remapped():
    log_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'resources', 'remapped_counts.yml'))
    assert count_variants_remapped(log_file) == (7147, 7002, 170)


def test_count_variants_extracted():
    log_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'resources', 'vcf_extractor.log'))
    assert count_variants_extracted(log_file) == (7147, 7147, 0, 0)


def test_count_variants_ingested():
    log_file = os.path.abspath(os.path.join(os.path.dirname(__file__), 'resources', 'vcf_ingestion.log'))
    assert count_variants_ingested(log_file) == (7002, 7002, 0)
