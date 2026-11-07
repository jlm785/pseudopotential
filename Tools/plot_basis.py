#!/usr/bin/env python3
"""Plots the atomic basis functions written by kb_conv (or atom_all) to
plot_kb.dat, with the style of atom_plot_mpl.py.  In plot_kb.dat each basis
function is followed by a marker line with its set, zeta, kind, nodes, r_c and
energy.  The debugging file basis_orbitals.dat (full radial grid, written if the
!dbg lines of atom_kb_orbitals are uncommented) is also accepted.

usage:
    python3 plot_basis.py                     all sets and l of plot_kb.dat
    python3 plot_basis.py --set 1             only basis set 1
    python3 plot_basis.py --format svg        SVG instead of PDF
    python3 plot_basis.py --dir some/dir      file in another directory
    python3 plot_basis.py --file other.dat    other file name (for example basis_orbitals.dat)

Output is basis_s<set>_l<l>_mpl.pdf (or .svg).  Confined orbitals (psi) and
perturbative polarization orbitals (pol) are solid lines, relaxed energy
derivatives (phi) are dashed, the color follows the zeta index.  Edit STYLE in
atom_plot_mpl.py to change the look.
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from atom_plot_mpl import STYLE, LINESTYLES, apply_style, color, plt


LETTERS = 'spdfghi'


def read_plot_kb(lines):
    """Reads plot_kb.dat: the data of each curve precede its marker line.
    Basis functions have markers 'marker bas<l letter> set S zeta Z kind K
    nodes N r_c R E E'."""
    blocks = []
    data = []
    for line in lines:
        t = line.split()
        if not t:
            continue
        if t[0] == 'marker':
            if t[1].startswith('bas') and 'set' in t:
                kv = dict(zip(t[2::2], t[3::2]))
                blocks.append({'set': int(kv['set']), 'l': LETTERS.index(t[1][3]),
                               'zeta': int(kv['zeta']), 'kind': kv['kind'],
                               'nodes': int(kv['nodes']), 'rc': float(kv['r_c']),
                               'E': float(kv['E']), 'data': np.array(data)})
            data = []
        else:
            data.append([float(x) for x in t[:2]])
    return blocks


def read_basis_orbitals(lines):
    """Reads basis_orbitals.dat: blocks starting with a '# set' header."""
    blocks = []
    cur = None
    for line in lines:
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


def read_basis(filename):
    """Basis functions from plot_kb.dat or basis_orbitals.dat."""
    with open(filename) as f:
        lines = f.read().split('\n')
    if any(line.startswith('# set') for line in lines):
        return read_basis_orbitals(lines)
    return read_plot_kb(lines)


def plot_set_l(blocks, iset, l, dirname, fmt):
    fig, ax = plt.subplots(figsize=STYLE['figsize'])
    rmax = 0.0
    for b in blocks:
        r, rpsi = b['data'][:, 0], b['data'][:, 1]
        style = 'dashed' if b['kind'] == 'phi' else 'solid'
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
    parser.add_argument('--file', default='plot_kb.dat',
                        help='file name (default plot_kb.dat)')
    args = parser.parse_args()

    path = os.path.join(args.dir, args.file)
    if not os.path.isfile(path):
        sys.exit(f'  no file {path}')
    blocks = read_basis(path)
    if not blocks:
        sys.exit(f'  no labelled basis functions in {path} (file from an older version?)')

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
