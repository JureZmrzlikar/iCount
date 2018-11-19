# pylint: disable=missing-docstring, protected-access
import os
import unittest
import warnings

from iCount.plotting import plot_rnamap, plot_rnaheatmap, plot_combined
from iCount.tests.utils import get_temp_file_name, make_file_from_list


class TestPlotRnaMap(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", (ResourceWarning, ImportWarning))

    def test_normalize_cpm(self):

        normalized = plot_rnamap.normalize_cpm(1, 10**6)
        self.assertEqual(normalized, 1)

        normalized = plot_rnamap.normalize_cpm(2, 10**6 * 2)
        self.assertEqual(normalized, 1)

    def test_parse(self):
        fname = make_file_from_list(bedtool=False, data=[
            ['total_cdna:1000000'],
            ['chr1_+_100', '1:3', '2:5'],
            ['chr1_+_200', '-1:3', '2:1'],
            ['chr1_+_300', '-51:3'],
            ['chr2_+_100', '100:8'],
        ])
        data, landmark_count = plot_rnamap.parse_results(fname, up_limit=-50, down_limit=50)
        self.assertEqual(landmark_count, 2)
        self.assertEqual(data, {
            -1: 3,
            1: 3,
            2: 6,
        })

    def test_plot(self):
        outfile = get_temp_file_name(extension='png')
        fname = make_file_from_list(bedtool=False, data=[
            ['total_cdna:1000000'],
            ['chr1_+_100', '1:3', '2:5'],
            ['chr1_+_200', '-1:3', '2:1'],
        ])
        plot_rnamap.plot_rnamap(fname, outfile)
        self.assertTrue(os.path.isfile(outfile))


class TestSmooth(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", (ResourceWarning, ImportWarning))

    def test_smoothe(self):
        # pylint: disable=missing-docstring, protected-access
        values = [0, 2, 0]
        self.assertEqual(plot_rnamap.smooth(values, 1), [1, 2 / 3, 1])

        values = [0, 2, 4, 1, 0]
        self.assertEqual(plot_rnamap.smooth(values, 1), [1, 2, 7 / 3, 5 / 3, 0.5])
        self.assertEqual(plot_rnamap.smooth(values, 2), [2, 7 / 4, 7 / 5, 7 / 4, 5 / 3])


class TestPlotRnaHeatMap(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", (ResourceWarning, ImportWarning))

    def test_make_position_to_bin1(self):
        pos2bins = plot_rnaheatmap.make_position_to_bin([0, 2, 4])
        self.assertEqual(pos2bins, {0: 0, 1: 0, 2: 2, 3: 2, 4: 2})

    def test_parse(self):
        fname = make_file_from_list(bedtool=False, data=[
            ['total_cdna:1000000'],
            ['chr1_+_100', '1:3', '12:5'],
            ['chr1_+_200', '-1:3', '2:1'],
            ['chr1_+_300', '-51:3'],
            ['chr2_+_100', '100:8'],
        ])
        data = plot_rnaheatmap.parse_results(fname, up_limit=-50, down_limit=50, top_n=2, binsize=10)
        self.assertEqual(list(data.index), ['chr1_+_100', 'chr1_+_200'])
        self.assertEqual(list(data.columns), list(range(-50, 50, 10)))
        self.assertEqual(data.at['chr1_+_100', 0], 3.0)
        self.assertEqual(data.at['chr1_+_100', 10], 5.0)
        self.assertEqual(data.at['chr1_+_200', -10], 3.0)
        self.assertEqual(data.at['chr1_+_200', 0], 1.0)

    def test_plot(self):
        outfile = get_temp_file_name(extension='png')
        fname = make_file_from_list(bedtool=False, data=[
            ['total_cdna:1000000'],
            ['chr1_+_100', '1:3', '2:5'],
            ['chr1_+_200', '-1:3', '2:1'],
        ])
        plot_rnaheatmap.plot_rnaheatmap(fname, outfile, top_n=2, binsize=10)
        self.assertTrue(os.path.isfile(outfile))


class TestPlotRnaCombined(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore", (ResourceWarning, ImportWarning))

    def test_plot(self):
        outfile = get_temp_file_name(extension='png')
        fname = make_file_from_list(bedtool=False, data=[
            ['total_cdna:1000000'],
            ['chr1_+_100', '1:3', '2:5'],
            ['chr1_+_200', '-1:3', '2:1'],
        ])
        plot_combined.plot_combined(fname, outfile, top_n=2, nbins=50)
        self.assertTrue(os.path.isfile(outfile))
