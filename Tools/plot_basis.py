#!/usr/bin/env python3
"""Plots the atomic basis functions written by kb_conv (atom_kb_orbitals)
to basis_orbitals.dat, with the style of atom_plot_mpl.py.

usage:
    python3 plot_basis.py                     all sets and l of basis_orbitals.dat
    python3 plot_basis.py --set 1             only basis set 1
    python3 plot_basis.py --format svg        SVG instead of PDF
    python3 plot_basis.py --dir some/dir      file in another directory
    python3 plot_basis.py --file other.dat    other file name

Output is basis_s<set>_l<l>_mpl.pdf (or .svg).  Confined orbitals (psi)
are solid lines, relaxed energy derivatives (phi) are dashed, the color
follows the zeta index.  Edit STYLE in atom_plot_mpl.py to change the look.
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from atom_plot_mpl import STYLE, LINESTYLES, apply_style, color, plt


def read_basis(filename):
    """Reads basis_orbitals.dat: blocks starting with a '# set' header."""
    blocks = []
    cur = None
    with open(filename) as f:
        for line in f:
            if line.startswith('# set'):
                t = line.split()
                cur = {'set': int(t[2]), 'l': int(t[4]), 'zeta': int(t[6]),
                       'kind': t[8], 'nodes': int(t[10]), 'rc': float(t[12]),
                       'E': float(t[15]), 'data': []}
                blocks.append(cur)
            elif line.strip() and not line.startswith('#') and cur is not None:
                cur['data'].append([float(x) for x in line.split()[:2]])
    for b in blocks:
        b['data'] = np.array(b['data'])
    return blocks


def plot_set_l(blocks, iset, l, dirname, fmt):
    fig, ax = plt.subplots(figsize=STYLE['figsize'])
    rmax = 0.0
    for b in blocks:
        r, rpsi = b['data'][:, 0], b['data'][:, 1]
        style = 'solid' if b['kind'] == 'psi' else 'dashed'
        ax.plot(r, rpsi, color=color(b['zeta']), linestyle=LINESTYLES[style],
                label=f"zeta {b['zeta']} {b['kind']}  r_c={b['rc']:.2f}  "
                      f"E={b['E']:.3f} Ry")
        rmax = max(rmax, r[-1])
    ax.axhline(0.0, color=STYLE['zero_line_color'],
               linewidth=STYLE['zero_line_width'], zorder=0)
    ax.set_xlim(0.0, 1.05 * rmax)
    ax.set_xlabel('r [a.u.]')
    ax.set_ylabel('rR(r)')
    if STYLE['show_title']:
        ax.set_title(f'Basis Functions, set {iset}, l = {l}')
    ax.legend(loc='best')
    fig.tight_layout()
    out = os.path.join(dirname, f'basis_s{iset}_l{l}_mpl.{fmt}')
    fig.savefig(out)
    plt.close(fig)
    print(f'  {out}')


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--set', type=int, default=None,
                        help='basis set to plot (default: all)')
    parser.add_argument('--format', default='pdf', choices=['pdf', 'svg'],
                        help='output format (default pdf)')
    parser.add_argument('--dir', default='.', help='directory with the file')
    parser.add_argument('--file', default='basis_orbitals.dat',
                        help='file name (default basis_orbitals.dat)')
    args = parser.parse_args()

    path = os.path.join(args.dir, args.file)
    if not os.path.isfile(path):
        sys.exit(f'  no file {path}')
    blocks = read_basis(path)

    apply_style()
    sets = sorted({b['set'] for b in blocks})
    if args.set is not None:
        sets = [s for s in sets if s == args.set]
    for iset in sets:
        for l in sorted({b['l'] for b in blocks if b['set'] == iset}):
            plot_set_l([b for b in blocks if b['set'] == iset and b['l'] == l],
                       iset, l, args.dir, args.format)


if __name__ == '__main__':
    main()
