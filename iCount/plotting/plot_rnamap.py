""".. Line to protect from pydocstyle D205, D400.

Plot distribution RNA-map
-------------------------

Plot distribution of crosslinks relative to landmark of specific type.
"""
import csv
import os

import iCount

# pylint: disable=wrong-import-order
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
from matplotlib import pyplot as plt  # pylint: disable=wrong-import-position
# pylint: enable=wrong-import-order


def normalize_cpm(value, total):
    """Normalize by CPM."""
    return value / total * 10**6


def parse_results(fname, up_limit, down_limit):
    """Parse RNA-maps results file."""
    data = {}
    landmark_count = 0
    with open(fname, 'rt') as handle:
        reader = csv.reader(handle, delimiter='\t')
        total_cdna_line = next(reader)
        total_cdna = int(total_cdna_line[0].strip().split(':')[1])
        for row in reader:
            is_contributing = False
            # First item is landmark name
            for item in row[1:]:
                pos, score = item.split(':')
                pos, score = int(pos), int(score)
                if up_limit <= pos <= down_limit:
                    is_contributing = True
                    data[pos] = data.get(pos, 0) + score

            if is_contributing:
                landmark_count += 1

    # Perform CPM normalization
    normalized_data = {pos: normalize_cpm(score, total_cdna) for pos, score in data.items()}

    return normalized_data, landmark_count


def smooth(list_, half_window=1):
    """Use box averaging smoothing."""
    new_list = []
    for i, _ in enumerate(list_):
        jjs = [j for j in range(i - half_window, i + half_window + 1) if 0 <= j < len(list_)]
        new_list.append(sum([list_[k] for k in jjs]) / len(jjs))

    return new_list


def make_outfile_name(fname, imgfmt):
    """Make output filename."""
    basename = iCount.files.remove_extension(fname, ['.tsv'])
    dirname = os.path.dirname(fname)
    return os.path.join(dirname, '{}_distro.{}'.format(basename, imgfmt))


def guess_maptype(fname):
    """Try to get RNA maptype from filename."""
    for mtype in iCount.analysis.rnamaps.MAP_TYPES:
        if fname.endswith('{}.tsv'.format(mtype)):
            return mtype


def plot_rnamap(fnames,
                outfile=None,
                up_limit=100,
                down_limit=100,
                ax=None,
                ylim=None,
                imgfmt='png',
                smoothing=1,
                ):
    """
    Plot distribution RNA-map.

    Parameters
    ----------
    fnames : list_str
        List of rnamaps result files to plot.
    outfile : str
        Output file.
    up_limit : int
        Upstream plot limit.
    down_limit : int
        Sownstream plot limit.
    ax : str
        An ``matplotlib.axes.Axes`` instance onto which this plot can
        be drawn. This is useful if you would like to use this function
        to plot this image as a subplot of a more complex figure.
    ylim : int
        Limit of the y-axis.
    imgfmt : str
        Image format. Note that image format is automatically
        determined from outfile. This parameters only applies if
        outfile is None.
    smoothing : int
        Smoothing half-window. Average smoothing is used.

    Returns
    -------
    None

    """
    # Make sure llimits have correct signs.
    up_limit = -abs(int(up_limit))
    down_limit = abs(int(down_limit))

    if not isinstance(fnames, list):
        fnames = [fnames]

    if not outfile:
        outfile = make_outfile_name(fnames[0], imgfmt)

    # User can provide axes instance (subplot) into which to plot this heatmap.
    # If not given a figure and axes instances are created.
    if ax:
        is_independent = False
    else:
        is_independent = True
        fig = plt.figure()
        ax = plt.subplot(1, 1, 1)

    for fname in fnames:
        data, landmark_count = parse_results(fname, up_limit, down_limit)
        if not data:
            continue

        positions, scores = zip(*sorted(data.items()))
        label = '{} ({} landmarks)'.format(
            guess_maptype(os.path.basename(fname)),
            landmark_count,
        )
        ax.plot(positions, smooth(scores, smoothing), label=label)

    ax.set_xlim((up_limit, down_limit))
    ax.set_xlabel('Position')
    if ylim:
        ax.set_ylim((0, ylim))
    ax.set_ylabel('Score [CPM]')
    ax.set_title('RNA-map')
    ax.grid(b=True, which='major', axis='both')
    ax.legend()

    if is_independent:
        fig.savefig(outfile)
