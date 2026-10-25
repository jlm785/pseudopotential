"""Reads a psp8 file following abinit's psp8in (src/64_psp/m_psp8.F90) and checks it.
usage: python3 check_psp8.py file.psp8 [old_file.psp8]"""
import sys
import numpy as np

def tokens(line):
    return line.replace('D', 'E').replace('d', 'e').split()

def read_psp8(fname):
    lines = open(fname).read().splitlines()
    pos = 0
    def rec():
        nonlocal pos
        line = lines[pos]; pos += 1
        return tokens(line)
    p = {}
    p['title'] = lines[0]; pos = 1
    t = rec(); p['zatom'], p['zion'] = float(t[0]), float(t[1])
    t = rec(); pspcod, pspxc, lmax, lloc, mmax = map(int, t[:5])
    assert pspcod == 8, 'pspcod'
    p.update(pspxc=pspxc, lmax=lmax, lloc=lloc, mmax=mmax)
    t = rec(); rchrg, fchrg = float(t[0]), float(t[1])
    t = rec(); nproj = [int(x) for x in t[:lmax+1]]
    t = rec(); ext = int(t[0])
    nprojso = [0]*lmax
    if ext in (2, 3):
        t = rec(); nprojso = [int(x) for x in t[:lmax]]
    elif ext not in (0, 1):
        raise ValueError('invalid extension_switch')
    nproj_tmp = nproj + nprojso
    if lloc < 4 and lloc + 1 <= len(nproj_tmp) and nproj_tmp[lloc] != 0:
        raise ValueError(f'nproj={nproj_tmp[lloc]} for the local l={lloc} (abinit stops)')
    nso = 2 if ext in (2, 3) else 1
    blocks = {}
    vloc = None
    for nn in range(1, nso+1):
        for ll in (range(0, lmax+1) if nn == 1 else range(1, lmax+1)):
            npr = nproj[ll] if nn == 1 else nprojso[ll-1]
            if npr > 0:
                t = rec(); assert int(t[0]) == ll, f'expected l={ll}, got {t[0]}'
                ekb = [float(x) for x in t[1:1+npr]]
                d = np.array([[float(x) for x in rec()[:2+npr]] for _ in range(mmax)])
                blocks[(nn, ll)] = (np.array(ekb), d[:, 1], d[:, 2:])
            elif ll == lloc and nn == 1:
                t = rec(); assert int(t[0]) == ll
                d = np.array([[float(x) for x in rec()[:3]] for _ in range(mmax)])
                vloc = d[:, 1:]
        if nn == 1 and lloc > lmax:
            t = rec(); assert int(t[0]) == lloc, f'general local potential label {t[0]} != lloc {lloc}'
            d = np.array([[float(x) for x in rec()[:3]] for _ in range(mmax)])
            vloc = d[:, 1:]
    core = None
    if fchrg > 1e-15:
        core = np.array([[float(x) for x in rec()[:7]] for _ in range(mmax)])
    vale = None
    if ext in (1, 3):
        vale = np.array([[float(x) for x in rec()[:3]] for _ in range(mmax)])
    p.update(rchrg=rchrg, fchrg=fchrg, nproj=nproj, nprojso=nprojso, ext=ext,
             blocks=blocks, vloc=vloc, core=core, vale=vale, leftover=len(lines)-pos)
    return p

def simpson(f, h):
    return h/3*(f[0] + f[-1] + 4*f[1:-1:2].sum() + 2*f[2:-1:2].sum())

def check(fname, old=None):
    p = read_psp8(fname)
    print(f'{fname}: lmax={p["lmax"]} lloc={p["lloc"]} mmax={p["mmax"]} ext={p["ext"]} '
          f'nproj={p["nproj"]} nprojso={p["nprojso"]} fchrg={p["fchrg"]} leftover lines={p["leftover"]}')
    r = p['vloc'][:, 0]; h = r[1]-r[0]
    assert r[0] == 0 and np.abs(np.diff(r) - h).max() < 1e-8, 'grid not linear from 0'
    for (nn, ll), (ekb, rr, pr) in sorted(p['blocks'].items()):
        norms = [simpson(pr[:, k]**2, h) for k in range(pr.shape[1])]
        print(f'   {"SR" if nn == 1 else "SO"} l={ll}: ekb={np.round(ekb, 6)}  int p^2 dr={np.round(norms, 8)}')
    if p['vale'] is not None:
        rv = p['vale'][:, 2]
        print(f'   valence: int 4 pi rho r^2 dr = {simpson(rv*r**2, h):.6f}   (zion = {p["zion"]})')
    if p['core'] is not None:
        print(f'   core: 4 pi rho_c(0) = {p["core"][0, 2]:.6f}, min = {p["core"][:, 2].min():.2e}')
    # spin-orbit algebra: SR + <L.S> SO must give pure j operators
    if p['ext'] in (2, 3):
        for ll in range(1, p['lmax']+1):
            if (1, ll) not in p['blocks'] or (2, ll) not in p['blocks']:
                continue
            esr, _, psr = p['blocks'][(1, ll)]; eso, _, pso = p['blocks'][(2, ll)]
            same = np.abs(psr - pso).max()
            op = lambda e, P: (P * e) @ P.T * h
            vp = op(esr, psr) + ll/2 * op(eso, pso)              # j = l+1/2
            vm = op(esr, psr) - (ll+1)/2 * op(eso, pso)          # j = l-1/2
            ep = esr[0] + ll/2*eso[0]; em = esr[1] - (ll+1)/2*eso[1]
            errp = np.abs(vp - ep*np.outer(psr[:, 0], psr[:, 0])*h).max() / abs(ep*h)
            errm = np.abs(vm - em*np.outer(psr[:, 1], psr[:, 1])*h).max() / abs(em*h)
            print(f'   l={ll}: projectors SR==SO {same:.1e};  E(j=l+1/2)={ep:.6f}  E(j=l-1/2)={em:.6f} Ha;'
                  f'  rank-1 errors {errp:.1e} {errm:.1e}')
    if old:
        q = read_psp8(old)
        for ll in range(p['lmax']+1):
            if (1, ll) in p['blocks'] and (1, ll) in q['blocks']:
                e1, _, p1 = p['blocks'][(1, ll)]; e2, _, p2 = q['blocks'][(1, ll)]
                o1 = (p1*e1) @ p1.T; o2 = (p2*e2) @ p2.T
                print(f'   l={ll}: operator vs old file, max rel. diff {np.abs(o1-o2).max()/np.abs(o2).max():.1e}')
        print(f'   local potential vs old file: max diff {np.abs(p["vloc"][:, 1]-q["vloc"][:, 1]).max():.1e}')
        if p['core'] is not None:
            print(f'   core vs old file: max diff {np.abs(p["core"][:, 2:]-q["core"][:, 2:]).max():.1e}')

if __name__ == '__main__':
    check(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
