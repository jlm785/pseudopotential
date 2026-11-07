# Files and Fortran units ("tapes") of the code

Source: `WORK/Src` (commit f337cf3 plus the changes of 7 November 2026: consistent units >= 10, TAUC files only for debugging). *Xy* is the chemical symbol.
Most routines receive the unit and the file name as arguments; the values are set in
the main programs (`atom.f90`, `psd_gen.f90`, `kb_conv.f90`, `atom_test.f90`, `kb_test.f90`,
`atom_all.f90`, `plot_show.f90`, `plot_kb_show.f90`, `plot_ln_show.f90`).

## Unit numbers (the same in every program)

| unit | variable | file | format | written by | read by |
|---|---|---|---|---|---|
| 5, 6 | | standard input and output | | | |
| 10 | `iowrite` (`atom_all`) | `atom.out` (the other programs write to unit 6) | text | `atom_all` | |
| 11 | `ioread` | `atom.dat` | text | `atom`, `atom_all` (if absent) | all except the plot programs |
| 12 | `ioae`, `iodata` | `datafile.dat` | unformatted | `atom`, `atom_all` | `psd_gen`, `plot_ln_show` |
| 13 | `iopsd` | `pseudo.dat` | unformatted | `psd_gen`, `atom_all` | `atom_test`, `plot_ln_show` |
| 14 | `ioreal` | `psd.pot`, *Xy*`_POTRE.DAT` (PARSEC), *Xy*`.psf` (SIESTA) | text | `psd_gen`, `atom_all` | `kb_conv` (`psd.pot`) |
| 15 | `ioplot` | `plot.dat` | text | `psd_gen`, `atom_all` | `plot_show` |
| 16 | `iopsdkb` (`iopsd` in `kb_test`, `iokb` in `plot_ln_show`) | `pseudokb.dat` | unformatted | `kb_conv`, `atom_all` | `atom_test`, `kb_test`, `plot_ln_show` |
| 17 | `ioplotkb` (`ioplot` in `plot_kb_show`) | `plot_kb.dat` (also the basis functions, each marker line with set, zeta, kind, nodes, r_c, E) | text | `kb_conv`, `atom_all` | `plot_kb_show`, `Tools/plot_basis.py` |
| 18 | `iokb` | *Xy*`_POTKB_F.DAT` (cpw2000) | text | `kb_conv`, `atom_all` | |
| 19 | `ioupf` | *Xy*`_LDA_TM.UPF` (QE), *Xy*`.psp8` (ABINIT), *Xy*`.psml` (SIESTA), one after the other | text | `kb_conv`, `atom_all` | |
| 20 | `iobas` | `basis_orbitals.dat` (debugging only, `!dbg` lines in `lib_kb/atom_kb_orbitals.f90`) | text | `kb_conv`, `atom_all` | `Tools/plot_basis.py --file basis_orbitals.dat` |
| 21 | `iocomm` | `command_`*name*`.gp`, `fig_`*name*`.dat` | text | plotting | gnuplot, `Tools/atom_plot_mpl.py` |
| 23 | `iotmp` | `tmp_mmga_x`, `tmp_mmga_t`, `tmp_mmga_c` (temporary) | text | plotting | |
| 30 | `iost` (`ioplot` in `plot_ln_show`) | `log_true.gp`, `log_ppot.gp`, `log_ppkb.gp`, and `command_log*.gp`, `fig_log*.dat` | text | `atom_all`, `plot_ln_show` | gnuplot |
| 31-34 | `iost+1`, ..., `iost+4` | `wfct.gp`, `fwft.gp`, `ppot.gp`, `fppt.gp`; `kbop.gp`, `kbtr.gp`, `kbek.gp`, `kbvk.gp`, `bas`*N*`.gp` (31) | text | `atom_all`, `plot_show`, `plot_kb_show` | gnuplot |

*name*`.pdf` files are written by gnuplot (no Fortran unit).

## Where the names are defined

| file | defined in |
|---|---|
| `atom.dat` | `atom.f90:39`, `psd_gen.f90:66`, `atom_test.f90:34`, `kb_test.f90:21`, `atom_all.f90:148` |
| `atom.out` | `atom_all.f90:183` (opened only if `iowrite /= 6`; same code in `atom.f90:46`, `atom_test.f90:40`, `kb_test.f90:27`) |
| `datafile.dat` | `atom.f90:42`, `psd_gen.f90:72`, `atom_all.f90:151`, `plot_ln_show.f90:65` |
| `pseudo.dat` | `atom.f90:36`, `psd_gen.f90:69`, `atom_test.f90:31`, `atom_all.f90:145`, `plot_ln_show.f90:68` |
| `psd.pot` | `psd_gen.f90:75`, `kb_conv.f90:74`, `atom_all.f90:154` |
| *Xy*`_POTRE.DAT`, *Xy*`.psf` | suffixes `psd_gen.f90:76-77`, `atom_all.f90:155-156`; names `lib_psd/atom_psd_sub.f90:381,389` |
| `plot.dat` | `psd_gen.f90:80`, `atom_all.f90:159`, `plot_show.f90:38` |
| `pseudokb.dat` | `kb_conv.f90:84`, `atom_test.f90:37`, `kb_test.f90:24`, `atom_all.f90:169`, `plot_ln_show.f90:70` |
| `plot_kb.dat` | `kb_conv.f90:87`, `atom_all.f90:172`, `plot_kb_show.f90:38` |
| `basis_orbitals.dat` (debugging) | `kb_conv.f90:90`, `atom_all.f90:175` (`iobas`, `filebas`, passed to `atom_kb_sub` and `atom_kb_orbitals`) |
| *Xy*`_POTKB_F.DAT`, `.UPF`, `.psp8`, `.psml` | suffixes `kb_conv.f90:77-81`, `atom_all.f90:162-166`; names `lib_kb/atom_kb_sub.f90` (`filekb = ...`, functional tag `scorr`) |
| plot files | `lib_plot/atom_plot_outdat.f90`, `atom_plot_kb_outdat.f90`, `atom_plot_ln_outpdf.f90`, `atom_plot_command.f90:49`, `atom_plot_spec.f90:78`, `atom_plot_gnuplot.f90`, `atom_plot_gnuplot_ok.f90:36` |

## Debugging output (commented with `!dbg`)

- `TAUC.dat` (core kinetic energy density, also the last record of `datafile.dat`) and
  `TAUC_orbitals.dat` (per orbital) in `lib_atom/atom_atm_datout.f90`; uncomment the `!dbg`
  lines to write them (unit 12, after `datafile.dat` is closed).
- `basis_orbitals.dat` (basis functions on the full radial grid) in `lib_kb/atom_kb_orbitals.f90`,
  unit `iobas` = 20.  The same functions, every 0.05 bohr, are in `plot_kb.dat`.
