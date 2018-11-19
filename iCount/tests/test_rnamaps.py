# pylint: disable=missing-docstring, protected-access
import os
import unittest
import warnings

from iCount.analysis import rnamaps
from iCount.files import remove_extension
from iCount.tests.utils import get_temp_file_name, make_file_from_list, make_list_from_file


class TestMakeLandmarksFile(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", (ResourceWarning, ImportWarning))

    def test_basic(self):
        regions = make_file_from_list([
            ['chr1', '.', 'CDS', '150', '200', '.', '+', '.', ''],
            ['chr1', '.', 'intron', '201', '351', '.', '+', '.', ''],
        ])
        fn = rnamaps.make_landmarks_file(regions, 'exon-intron')
        self.assertEqual(make_list_from_file(fn), [
            ['chr1', '200', '201', '.', '.', '+'],
        ])

        regions = make_file_from_list([
            ['chr1', '.', 'CDS', '150', '200', '.', '-', '.', ''],
            ['chr1', '.', 'intron', '201', '351', '.', '-', '.', ''],
        ])
        fn = rnamaps.make_landmarks_file(regions, 'intron-exon')
        self.assertEqual(make_list_from_file(fn), [
            ['chr1', '199', '200', '.', '.', '-'],
        ])

    def test_limits_upstream(self):
        """Landmarks with too short upstream segment should not be used."""
        regions = make_file_from_list([
            ['chr1', '.', 'CDS', '151', '200', '.', '+', '.', ''],
            ['chr1', '.', 'intron', '201', '400', '.', '+', '.', ''],
        ])
        fn = rnamaps.make_landmarks_file(regions, 'exon-intron')
        self.assertEqual(make_list_from_file(fn), [])

        regions = make_file_from_list([
            ['chr1', '.', 'CDS', '150', '200', '.', '-', '.', ''],
            ['chr1', '.', 'intron', '201', '350', '.', '-', '.', ''],
        ])
        fn = rnamaps.make_landmarks_file(regions, 'intron-exon')
        self.assertEqual(make_list_from_file(fn), [])

    def test_limits_downstream(self):
        """Landmarks with too short upstream segment should not be used."""
        regions = make_file_from_list([
            ['chr1', '.', 'CDS', '150', '200', '.', '+', '.', ''],
            ['chr1', '.', 'intron', '201', '350', '.', '+', '.', ''],
        ])
        fn = rnamaps.make_landmarks_file(regions, 'exon-intron')
        self.assertEqual(make_list_from_file(fn), [])

        regions = make_file_from_list([
            ['chr1', '.', 'CDS', '151', '200', '.', '-', '.', ''],
            ['chr1', '.', 'intron', '201', '351', '.', '-', '.', ''],
        ])
        fn = rnamaps.make_landmarks_file(regions, 'intron-exon')
        self.assertEqual(make_list_from_file(fn), [])


class TestComputeDistances(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)

        self.landmarks = make_file_from_list([
            ['chr1', '10', '11', '.', '.', '+'],
            ['chr1', '20', '21', '.', '.', '+'],
            ['chr1', '20', '21', '.', '.', '-'],
            ['chr2', '10', '11', '.', '.', '+'],
        ])

    def test_basic(self):
        xlinks = make_file_from_list([
            ['chr1', '12', '13', '.', '3', '+'],
        ])
        distances, total_cdna = rnamaps.compute_distances(self.landmarks, xlinks, max_width=10)
        self.assertEqual(total_cdna, 3)
        self.assertEqual(distances, {
            'chr1_+_10': {2: 3},
        })

    def test_scores_sum(self):
        xlinks = make_file_from_list([
            ['chr1', '12', '13', '.', '3', '+'],
            ['chr1', '12', '13', '.', '1', '+'],
        ])
        distances, total_cdna = rnamaps.compute_distances(self.landmarks, xlinks, max_width=10)
        self.assertEqual(total_cdna, 4)
        self.assertEqual(distances, {
            'chr1_+_10': {2: 4},
        })

    def test_strands_not_mixed(self):
        xlinks = make_file_from_list([
            ['chr1', '12', '13', '.', '3', '+'],
            ['chr1', '12', '13', '.', '1', '-'],
        ])
        distances, total_cdna = rnamaps.compute_distances(self.landmarks, xlinks, max_width=10)
        self.assertEqual(total_cdna, 4)
        self.assertEqual(distances, {
            'chr1_+_10': {2: 3},
            'chr1_-_20': {8: 1},
        })

    def test_chroms_not_mixed(self):
        xlinks = make_file_from_list([
            ['chr1', '12', '13', '.', '3', '+'],
            ['chr2', '12', '13', '.', '1', '+'],
        ])
        distances, total_cdna = rnamaps.compute_distances(self.landmarks, xlinks, max_width=10)
        self.assertEqual(total_cdna, 4)
        self.assertEqual(distances, {
            'chr1_+_10': {2: 3},
            'chr2_+_10': {2: 1},
        })

    def test_max_width(self):
        xlinks = make_file_from_list([
            ['chr1', '22', '23', '.', '3', '+'],
            ['chr2', '32', '33', '.', '1', '+'],
        ])
        distances, total_cdna = rnamaps.compute_distances(self.landmarks, xlinks, max_width=10)
        self.assertEqual(total_cdna, 4)
        self.assertEqual(distances, {
            'chr1_+_20': {2: 3},
        })

    def test_no_landmark(self):
        """Landmark is missing on this chromosome / stramd."""
        xlinks = make_file_from_list([
            ['chrX', '22', '23', '.', '3', '+'],
        ])
        distances, total_cdna = rnamaps.compute_distances(self.landmarks, xlinks, max_width=10)
        self.assertEqual(total_cdna, 3)
        self.assertEqual(distances, {})


class TestDistancesToFile(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)
        self.fname = get_temp_file_name()

    def test_distances_to_file(self):
        distances = {
            'chr1_+_10': {2: 3},
            'chr2_+_20': {-3: 3, 1: 5},
        }
        rnamaps.distances_to_file(distances, self.fname, total_cdna=11)
        self.assertEqual(make_list_from_file(self.fname, fields_separator='\t'), [
            ['total_cdna:11'],
            ['chr1_+_10', '2:3'],
            ['chr2_+_20', '-3:3', '1:5'],
        ])


class TestRun(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", ResourceWarning)
        self.outdir = get_temp_file_name()

    def test_run(self):
        regions = make_file_from_list(sort=True, data=[
            ['chr1', '.', 'intergenic', '1', '210', '.', '+', '.', ' '],
            ['chr1', '.', 'UTR5', '211', '270', '.', '+', '.', ' '],
            ['chr1', '.', 'CDS', '271', '330', '.', '+', '.', ' '],
            ['chr1', '.', 'intron', '331', '490', '.', '+', '.', ' '],
            ['chr1', '.', 'CDS', '491', '550', '.', '+', '.', ' '],
            ['chr1', '.', 'UTR3', '551', '760', '.', '+', '.', ' '],
            ['chr1', '.', 'intergenic', '761', '1100', '.', '+', '.', ' '],

            ['chr1', '.', 'intergenic', '1', '300', '.', '-', '.', ' '],
            ['chr1', '.', 'ncRNA', '301', '500', '.', '-', '.', ' '],
            ['chr1', '.', 'intron', '501', '600', '.', '-', '.', ' '],
            ['chr1', '.', 'ncRNA', '601', '750', '.', '-', '.', ' '],
            ['chr1', '.', 'intergenic', '751', '1000', '.', '-', '.', ' '],
        ])

        sites = make_file_from_list([
            ['chr1', '120', '121', '.', '1', '+'],
            ['chr1', '350', '351', '.', '1', '+'],
            ['chr1', '550', '551', '.', '1', '+'],
            ['chr1', '750', '751', '.', '1', '-'],
        ])

        rnamaps.run(sites, regions, outdir=self.outdir)

        self.assertTrue(os.path.isdir(self.outdir))

        sites_name = remove_extension(sites, ['.bed', '.bed.gz'])
        for maptype in rnamaps.MAP_TYPES:
            fname = os.path.join(self.outdir, '{}_{}.tsv'.format(sites_name, maptype))
            self.assertTrue(os.path.isfile(fname))
            self.assertGreater(os.path.getsize(fname), 1)


if __name__ == '__main__':
    unittest.main()
