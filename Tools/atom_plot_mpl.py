#!/usr/bin/env python3
"""Publication plots with matplotlib for the atom/pseudopotential codes.

The Fortran plotting routines (lib_plot) write, next to the gnuplot files,
a plain text description of each figure, fig_<name>.dat, see
lib_plot/atom_plot_spec.f90.  This script draws those figures.

usage:
    python3 atom_plot_mpl.py                  all fig_*.dat in the directory
    python3 atom_plot_mpl.py wfct bas0        only those figures
    python3 atom_plot_mpl.py --format svg     SVG instead of PDF
    python3 atom_plot_mpl.py --dir some/dir   figures in another directory

Output is <name>_mpl.pdf (or .svg), so the gnuplot <name>.pdf is kept.
Text is kept as editable text in PDF and SVG, so labels can be changed
afterwards (e.g. LaTeX labels added with Inkscape).

Requires numpy and matplotlib.  Edit STYLE below to change the look.
"""

import argparse
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# Look of the figures.  Change here, not in the Fortran code.

STYLE = {
    'figsize': (6.0, 4.2),          # inches; (3.4, 2.6) for one journal column
    'font_size': 12,
    'line_width': 1.5,
    'show_title': True,             # journals usually prefer no title
    'legend_frame': False,
    'zero_line_color': '0.6',
    'zero_line_width': 0.6,
}

# Colors indexed by the gnuplot line style number of atom_plot_command
# (index 0 is used for auxiliary curves).  Colorblind-friendly cycle.

COLORS = ['0.3', '#0072B2', '#D55E00', '#009E73', '#CC79A7',
          '#E69F00', '#56B4E9', '#F0E442', '#000000', '#999999']

LINESTYLES = {'solid': '-', 'dashed': '--', 'dotted': ':', 'dashdot': '-.'}


def apply_style():
    """Sets matplotlib parameters from STYLE."""
    plt.rcParams.update({
        'font.size': STYLE['font_size'],
        'lines.linewidth': STYLE['line_width'],
        'axes.titlesize': STYLE['font_size'],
        'legend.fontsize': STYLE['font_size'] - 2,
        'legend.frameon': STYLE['legend_frame'],
        'pdf.fonttype': 42,         # editable text in PDF
        'svg.fonttype': 'none',     # editable text in SVG (Inkscape)
    })


def color(index):
    return COLORS[index % len(COLORS)]


def read_spec(filename):
    """Reads fig_<name>.dat into a dictionary."""
    spec = {'curves': [], 'functions': [], 'hlines': []}
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or '=' not in line:
                continue
            key, value = [s.strip() for s in line.split('=', 1)]
            if key in ('title', 'xlabel', 'ylabel'):
                spec[key] = value
            elif key in ('xrange', 'yrange'):
                spec[key] = [float(v) for v in value.split()]
            elif key == 'curve':
                fname, xcol, ycol, label, style, icolor = \
                    [s.strip() for s in value.split('|')]
                spec['curves'].append({'file': fname, 'x': int(xcol),
                                       'y': int(ycol), 'label': label,
                                       'style': style, 'color': int(icolor)})
            elif key in ('coulomb', 'centrifugal'):
                c, label, style, icolor = [s.strip() for s in value.split('|')]
                spec['functions'].append({'kind': key, 'c': float(c),
                                          'label': label, 'style': style,
                                          'color': int(icolor)})
            elif key == 'hline':
                spec['hlines'].append(float(value))
            else:
                print(f'  {filename}: unknown key {key}, ignored')
    return spec


_data_cache = {}


def read_data(dirname, fname):
    """Reads a gnuplot data file (columns, # comments), cached."""
    path = os.path.join(dirname, fname)
    if path not in _data_cache:
        _data_cache[path] = np.atleast_2d(np.loadtxt(path, comments='#'))
    return _data_cache[path]


def plot_figure(name, dirname, fmt):
    spec = read_spec(os.path.join(dirname, f'fig_{name}.dat'))

    fig, ax = plt.subplots(figsize=STYLE['figsize'])

    ymin, ymax = spec.get('yrange', (-np.inf, np.inf))
    span = ymax - ymin

    for c in spec['curves']:
        data = read_data(dirname, c['file'])
        y = data[:, c['y'] - 1].copy()
        # points far outside the y range are dropped, which breaks the
        # line at poles (log derivatives) instead of drawing a vertical line
        if np.isfinite(span):
            y[(y < ymin - span) | (y > ymax + span)] = np.nan
        ax.plot(data[:, c['x'] - 1], y, label=c['label'],
                color=color(c['color']), linestyle=LINESTYLES[c['style']])

    xmin, xmax = spec.get('xrange', ax.get_xlim())
    for f in spec['functions']:
        x = np.linspace(max(xmin, 1.0e-3 * (xmax - xmin)), xmax, 500)
        y = -f['c'] / x if f['kind'] == 'coulomb' else f['c'] / x**2
        ax.plot(x, y, label=f['label'], color=color(f['color']),
                linestyle=LINESTYLES[f['style']])

    for y in spec['hlines']:
        ax.axhline(y, color=STYLE['zero_line_color'],
                   linewidth=STYLE['zero_line_width'], zorder=0)

    if 'xrange' in spec:
        ax.set_xlim(*spec['xrange'])
    if 'yrange' in spec:
        ax.set_ylim(*spec['yrange'])
    ax.set_xlabel(spec.get('xlabel', ''))
    ax.set_ylabel(spec.get('ylabel', ''))
    if STYLE['show_title'] and 'title' in spec:
        ax.set_title(spec['title'])
    if spec['curves'] or spec['functions']:
        ax.legend(loc='best')

    fig.tight_layout()
    out = os.path.join(dirname, f'{name}_mpl.{fmt}')
    fig.savefig(out)
    plt.close(fig)
    print(f'  {out}')


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('names', nargs='*',
                        help='figures to plot (default: all fig_*.dat)')
    parser.add_argument('--format', default='pdf', choices=['pdf', 'svg'],
                        help='output format (default pdf)')
    parser.add_argument('--dir', default='.', help='directory with the files')
    args = parser.parse_args()

    names = args.names or sorted(
        os.path.basename(f)[4:-4]
        for f in glob.glob(os.path.join(args.dir, 'fig_*.dat')))
    if not names:
        sys.exit(f'  no fig_*.dat files in {args.dir}')

    apply_style()
    for name in names:
        try:
            plot_figure(name, args.dir, args.format)
        except (OSError, ValueError) as err:
            print(f'  {name}: not plotted ({err})')


if __name__ == '__main__':
    main()
