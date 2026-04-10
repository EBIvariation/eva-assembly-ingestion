import datetime
import os
import unittest
from unittest.mock import patch, MagicMock

from eva_assembly_ingestion.config import load_config
from eva_assembly_ingestion.assembly_ingestion_job import AssemblyIngestionJob


class TestAddToClusteredVariantUpdate(unittest.TestCase):
    resources_folder = os.path.join(os.path.dirname(__file__), 'resources')

    def setUp(self):
        config_file = os.path.join(self.resources_folder, 'remapping_config.yml')
        load_config(config_file)
        self.remapping_job = AssemblyIngestionJob(taxonomy=9913, target_assembly='GCA_000003055.3', release_version=5)

    def test_inserts_one_row_per_taxonomy_per_source_assembly(self):
        # Two tracker rows with two origin assemblies, each with two taxonomies → 4 inserts
        tracker_rows = [
            ('EVA', '9913,9940', 'Cattle', 'GCA_000000001.1', 'GCA_000003055.3', 10, 'Completed'),
            ('EVA', '9940', 'Cattle', 'GCA_000000002.1', 'GCA_000003055.3', 5, 'Completed'),
        ]
        mocked_now = datetime.datetime.now()
        with patch('eva_assembly_ingestion.assembly_ingestion_job.get_metadata_connection_handle'), \
                patch('eva_assembly_ingestion.assembly_ingestion_job.execute_query') as mock_execute_query, \
                patch.object(self.remapping_job, 'get_job_information_from_tracker', return_value=tracker_rows), \
                patch('eva_assembly_ingestion.assembly_ingestion_job.datetime') as mock_datetime:
            mock_datetime.datetime.now.return_value = mocked_now
            self.remapping_job.add_to_clustered_variant_update()
        captured_queries = [c[1][1] for c in mock_execute_query.mock_calls]
        expected_queries = [
            f"INSERT INTO evapro.clustered_variant_update (taxonomy_id, assembly_accession, source, ingestion_time) VALUES (9913, 'GCA_000003055.3', 'GCA_000000001.1', '{mocked_now.strftime('%Y-%m-%d %H:%M:%S.%f')}')",
            f"INSERT INTO evapro.clustered_variant_update (taxonomy_id, assembly_accession, source, ingestion_time) VALUES (9940, 'GCA_000003055.3', 'GCA_000000001.1', '{mocked_now.strftime('%Y-%m-%d %H:%M:%S.%f')}')",
            f"INSERT INTO evapro.clustered_variant_update (taxonomy_id, assembly_accession, source, ingestion_time) VALUES (9940, 'GCA_000003055.3', 'GCA_000000002.1', '{mocked_now.strftime('%Y-%m-%d %H:%M:%S.%f')}')"
        ]
        assert captured_queries == expected_queries
