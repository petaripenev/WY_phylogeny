
import unittest
from unittest.mock import patch, mock_open
from pandas.testing import assert_frame_equal
from build_ecosys_trees_and_annotate import *

class TestReadS3Clusters(unittest.TestCase):
    def setUp(self):
        self.mock_data = '''S	29	1311	*	*	*	*	*	WaterYear_AB7_D5_3235493_9	*
S	30	1311	*	*	*	*	*	WaterYear_AB9_D3_16863982_1	*
H	29	1311	99.0	+	0	0	1311M	WaterYear_HB5_D5_16973911_15	WaterYear_AB7_D5_3235493_9
H	29	1311	99.0	+	0	0	1311M	WaterYear_HB9_D5_9402319_5	WaterYear_AB7_D5_3235493_9
H	29	603	99.3	+	0	0	603M708I	WaterYear_HA9_D5_16727361_1	WaterYear_AB7_D5_3235493_9'''

    def test_read_S3_clusters(self):
        with patch(
                'builtins.open',
                mock_open(read_data=self.mock_data),
        ) as file_mock:
            result = read_S3_clusters(file_mock)
            res_keys = list(result.keys())
            self.assertEqual(len(result), 2)
            self.assertEqual(len(result[res_keys[0]]), 4)
            self.assertIn(res_keys[0], result[res_keys[0]])
            self.assertIn('WaterYear_AB9_D3_16863982_1', res_keys)
            self.assertIn('WaterYear_AB7_D5_3235493_9', res_keys)

class TestReadNRAnnotations(unittest.TestCase):
    def setUp(self):
        self.mock_data = [
            ['gene1', 'x', 'x', 'lineage1;lineage2', 'rank1;rank2'],
            ['gene2', 'x', 'x', 'lineage3;lineage4', 'rank3;rank4']
        ]
        self.expected_df = pd.DataFrame({
            'rank1': {'gene1': 'lineage1', 'gene2': None},
            'rank2': {'gene1': 'lineage2', 'gene2': None},
            'rank3': {'gene1': 'Unclassified lineage2', 'gene2': 'lineage3'},
            'rank4': {'gene1': 'Unclassified lineage2', 'gene2': 'lineage4'}
        })

    @patch('csv.reader')
    @patch('builtins.open', new_callable=mock_open)
    def test_read_NR_annotations(self, mock_file, mock_csv):
        mock_csv.return_value = self.mock_data
        result_df = read_NR_annotations('path_to_NR', ['rank1', 'rank2', 'rank3', 'rank4'])
        assert_frame_equal(result_df, self.expected_df)

class TestSplitClusters(unittest.TestCase):
    def setUp(self):
        self.s3_clusters_dict = {
            'WaterYear_AB9_D3_16863982_1': ['WaterYear_AB9_D3_16863982_1', 'WaterYear_HB5_D5_16973911_15', 'WaterYear_HB9_D5_9402319_5', 'WaterYear_HA9_D5_16727361_1'],
            'WaterYear_HB9_D5_9402320_5': ['WaterYear_HB9_D5_9402320_5']
        }
        self.expected_eco_dict = {
            'A': ['WaterYear_AB9_D3_16863982_1'],
            'H': ['WaterYear_AB9_D3_16863982_1', 'WaterYear_HB9_D5_9402320_5']
        }

    def test_split_clusters(self):
        result_eco_dict = split_clusters(self.s3_clusters_dict)
        self.assertEqual(result_eco_dict, self.expected_eco_dict)
        result_all_dict = split_clusters(self.s3_clusters_dict, do_not_split=True)
        self.assertEqual(result_all_dict, {'ALL': ['WaterYear_AB9_D3_16863982_1', 'WaterYear_HB9_D5_9402320_5']})

class TestParseEafPerEcosystem(unittest.TestCase):
    def setUp(self):
        self.eaf_dict = {
            'eco1': [
                {'cluster_gene': 'gene1', 'core': 'core1', 'mean_resampled_EAF': 0.1},
                {'cluster_gene': 'gene1', 'core': 'core2', 'mean_resampled_EAF': 0.2},
                {'cluster_gene': 'gene2', 'core': 'core1', 'mean_resampled_EAF': 0.3}
            ],
            'eco2': [
                {'cluster_gene': 'gene1', 'core': 'core1', 'mean_resampled_EAF': 0.4}
            ]
        }
        self.expected_eco_gene_eaf = {
            'eco1': {
                'gene1': {'core1': 0.1, 'core2': 0.2},
                'gene2': {'core1': 0.3}
            },
            'eco2': {
                'gene1': {'core1': 0.4}
            }
        }
        self.expected_eco_cores = {
            'eco1': ['core1', 'core2'],
            'eco2': ['core1']
        }

    def test_parse_eaf_per_ecosystem(self):
        from build_ecosys_trees_and_annotate import parse_eaf_per_ecosystem
        result_eco_gene_eaf, result_eco_cores = parse_eaf_per_ecosystem(self.eaf_dict)
        self.assertEqual(result_eco_gene_eaf, self.expected_eco_gene_eaf)
        self.assertEqual(result_eco_cores, self.expected_eco_cores)

if __name__ == '__main__':
    unittest.main()