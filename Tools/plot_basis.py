"""Plots the basis functions of basis_orbitals.dat (written by kb_conv).
usage: python3 plot_basis.py basis_orbitals.dat out.png [set]"""
import sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
fn, out = sys.argv[1], sys.argv[2]
iset = int(sys.argv[3]) if len(sys.argv) > 3 else 1
B=[];cur=None
for line in open(fn):
    if line.startswith('# set'):
        t=line.split(); cur={'set':int(t[2]),'l':int(t[4]),'z':int(t[6]),'kind':t[8],'rc':float(t[12]),'E':float(t[15]),'d':[]}; B.append(cur)
    elif line.strip() and not line.startswith('#'): cur['d'].append([float(x) for x in line.split()])
B=[b for b in B if b['set']==iset]
ls=sorted({b['l'] for b in B})
fig,axs=plt.subplots(1,len(ls),figsize=(5.5*len(ls),4.5),squeeze=False)
for ax,l in zip(axs[0],ls):
    for b in [b for b in B if b['l']==l]:
        d=np.array(b['d'])
        ax.plot(d[:,0],d[:,1],ls='-' if b['kind']=='psi' else '--',
                label=f"zeta {b['z']} {b['kind']}  rc={b['rc']:.2f}  E={b['E']:.3f} Ry")
    ax.axhline(0,color='gray',lw=0.5); ax.set_title(f"set {iset}, l={l}"); ax.set_xlabel('r (bohr)'); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(out,dpi=100)
