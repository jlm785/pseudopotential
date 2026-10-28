"""Reads a psml file written by the atom code and checks it.
With a psp8 file of the same pseudopotential, compares the two.

usage: python3 check_psml.py file.psml [file.psp8]

Checks:
  - normalization of projectors (int r^2 chi^2 dr = 1) and pseudo
    wave-functions (int r^2 R^2 dr = 1), valence charge vs total-charge
  - relations between the sets: scalar_relativistic and spin_orbit
    semilocal potentials and projector operators from the lj set
  - with a psp8 file: KB energies, projectors r*chi, local potential,
    valence and core charges
"""
import sys
import xml.etree.ElementTree as ET
import numpy as np

sys.path.insert(0, __file__.rsplit('/', 1)[0] if '/' in __file__ else '.')
from check_psp8 import read_psp8

LCHAR = 'spdfghi'


def read_psml(fname):
    root = ET.parse(fname).getroot()
    ns = root.tag[1:].split('}')[0]
    q = lambda tag: f'{{{ns}}}{tag}'
    grid = np.array(root.find(f'{q("grid")}/{q("grid-data")}').text.split(), float)

    def radfunc(el):
        d = np.array(el.find(f'{q("radfunc")}/{q("data")}').text.split(), float)
        f = np.zeros(len(grid))
        f[:len(d)] = d
        return f

    p = {'grid': grid, 'root': root.attrib}
    spec = root.find(q('pseudo-atom-spec'))
    p['spec'] = spec.attrib
    p['shells'] = [s.attrib for s in spec.iter(q('shell'))]
    vc = root.find(q('valence-charge'))
    p['valence'] = (float(vc.get('total-charge')), radfunc(vc))
    cc = root.find(q('pseudocore-charge'))
    p['core'] = radfunc(cc) if cc is not None else None
    p['vloc'] = radfunc(root.find(q('local-potential')))
    p['slps'] = {}
    for s in root.iter(q('semilocal-potentials')):
        for v in s.iter(q('slps')):
            p['slps'][(s.get('set'), LCHAR.index(v.get('l')), v.get('j'))] = radfunc(v)
    p['proj'] = {}
    for s in root.iter(q('nonlocal-projectors')):
        for v in s.iter(q('proj')):
            key = (s.get('set'), LCHAR.index(v.get('l')), v.get('j'), int(v.get('seq')))
            p['proj'][key] = (float(v.get('ekb')), radfunc(v))
    p['pswf'] = {}
    for s in root.iter(q('pseudo-wave-functions')):
        for v in s.iter(q('pswf')):
            p['pswf'][(s.get('set'), LCHAR.index(v.get('l')), v.get('j'))] = radfunc(v)
    return p


def integ(f, r):
    return np.trapezoid(f, r) if hasattr(np, 'trapezoid') else np.trapz(f, r)


def operator(projs, r, rmax):
    """matrix of sum_i |r chi_i> ekb_i <r chi_i| on a coarse grid"""
    x = np.linspace(0, rmax, 400)
    m = np.zeros((len(x), len(x)))
    for ekb, chi in projs:
        u = np.interp(x, r, r*chi)
        m += ekb*np.outer(u, u)
    return m


def check(fname, fpsp8=None):
    p = read_psml(fname)
    r = p['grid']
    sp = p['spec']
    print(f'{fname}: {sp.get("atomic-label")} Z={sp.get("atomic-number")} z-pseudo={float(sp.get("z-pseudo"))}'
          f' relativity={sp.get("relativity")} core={sp.get("core-corrections")} grid npts={len(r)} r(1)={r[0]}')
    print('   shells: ' + ', '.join(f'{s["n"]}{s["l"]} {float(s["occupation"]):g}' for s in p['shells']))
    qt, qv = p['valence']
    print(f'   valence charge: int r^2 q dr = {integ(r*r*qv, r):.6f}  total-charge = {qt:.6f}')

    worst = 0.0
    for key, (ekb, chi) in sorted(p['proj'].items()):
        worst = max(worst, abs(integ(r*r*chi*chi, r) - 1))
    print(f'   projectors: {len(p["proj"])}, max |int r^2 chi^2 dr - 1| = {worst:.1e}')
    worst = 0.0
    for key, R in p['pswf'].items():
        worst = max(worst, abs(integ(r*r*R*R, r) - 1))
    print(f'   pseudo wave-functions: {len(p["pswf"])}, max |int r^2 R^2 dr - 1| = {worst:.1e}')

    if sp.get('relativity') == 'dirac':
        # semilocal: SR = ((l+1) V+ + l V-)/(2l+1), SO = 2 (V+ - V-)/(2l+1)
        e_sr = e_so = 0.0
        for (st, l, j), v in p['slps'].items():
            if st != 'scalar_relativistic' or l == 0:
                continue
            vp = next(f for (s2, l2, j2), f in p['slps'].items() if s2 == 'lj' and l2 == l and float(j2) > l)
            vm = next(f for (s2, l2, j2), f in p['slps'].items() if s2 == 'lj' and l2 == l and float(j2) < l)
            e_sr = max(e_sr, np.abs(v - ((l+1)*vp + l*vm)/(2*l+1)).max())
            so = p['slps'][('spin_orbit', l, None)]
            e_so = max(e_so, np.abs(so - 2*(vp - vm)/(2*l+1)).max())
        print(f'   semilocal sets from lj: max error SR {e_sr:.1e}  SO {e_so:.1e} Ha')
        # projector operators: SR + l/2 SO = V(j=l+1/2), SR - (l+1)/2 SO = V(j=l-1/2)
        lmax = max(l for (_, l, _, _) in p['proj'])
        for l in range(1, lmax+1):
            get = lambda st, cond=lambda j: True: [(e, c) for (s2, l2, j2, _), (e, c) in p['proj'].items()
                                                    if s2 == st and l2 == l and (j2 is None or cond(float(j2)))]
            sr, so = get('scalar_relativistic'), get('spin_orbit')
            ljp, ljm = get('lj', lambda j: j > l), get('lj', lambda j: j < l)
            if not (sr and so and ljp and ljm):
                continue
            rmax = 10.0
            osr, oso = operator(sr, r, rmax), operator(so, r, rmax)
            op, om = operator(ljp, r, rmax), operator(ljm, r, rmax)
            ep = np.abs(osr + l/2*oso - op).max()/np.abs(op).max()
            em = np.abs(osr - (l+1)/2*oso - om).max()/np.abs(om).max()
            print(f'   l={l}: projector operators SR+SO vs lj, rel. errors {ep:.1e} (j=l+1/2) {em:.1e} (j=l-1/2)')

    if fpsp8:
        q = read_psp8(fpsp8)
        x = q['vloc'][:, 0]
        n = len(x)
        st = 'scalar_relativistic' if sp.get('relativity') == 'dirac' else 'non_relativistic'
        print(f'   comparison with {fpsp8}:')
        for (nn, l), (ekb8, _, pr8) in sorted(q['blocks'].items()):
            sset = st if nn == 1 else 'spin_orbit'
            projs = sorted([(seq, e, c) for (s2, l2, j2, seq), (e, c) in p['proj'].items()
                            if s2 == sset and l2 == l])
            ek = np.array([e for _, e, _ in projs])
            dp = max(np.abs(np.interp(x, r, r*c) - pr8[:, k]).max() for k, (_, _, c) in enumerate(projs))
            print(f'      {"SR" if nn == 1 else "SO"} l={l}: ekb psml {np.round(ek, 6)}  psp8 {np.round(ekb8, 6)}'
                  f'  max|diff| {np.abs(ek-ekb8).max():.1e};  max|r chi - p| {dp:.1e}')
        rl = x[(x > 0.01) & (x < 0.8*x[-1])]
        dv = np.abs(np.interp(rl, r, p['vloc']) - np.interp(rl, x, q['vloc'][:, 1])).max()
        print(f'      local potential: max diff {dv:.1e} Ha (0.01 < r < {0.8*x[-1]:.1f})')
        if q['vale'] is not None:
            dq = np.abs(np.interp(x, r, qv) - q['vale'][:, 2]).max()
            print(f'      valence charge: max diff {dq:.1e} (max {np.abs(qv).max():.2e})')
        if q['core'] is not None and p['core'] is not None:
            dc = np.abs(np.interp(x, r, p['core']) - q['core'][:, 2]).max()
            print(f'      core charge: max diff {dc:.1e} (max {np.abs(p["core"]).max():.2e})')


if __name__ == '__main__':
    check(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
