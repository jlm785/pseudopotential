# pseudopotential
Modernized version of a pseudopotential generation code

This is a Fortran 90 version of the venerable atomic configuration and pseudopotential generation code, with extensive documentation and some new capabilities.

The initial code was written by Sverre Froyen at the University of California at Berkeley in the early 1980s.  Norm Troullier modified it to generate the Troullier-Martins pseudopotential at the University of Minnesota at Minneapolis.  It has been maintained by José Luís Martins and collaborators.  It has been on the web since early 1998 (https://web.archive.org/web/19980111014249/http://bohr.inesc.pt/~jlm/pseudo.html)

The 6.x.y version is a conversion from the old Fortran 77 to modern Fortran 90 with minimal changes to the core algorithms.

## What the code does

- All-electron atomic calculations: non-relativistic, spin-polarized, or relativistic (Dirac equation), with the LDA (Ceperley-Alder, Perdew-Wang) or the PBE GGA.
- Norm-conserving Troullier-Martins pseudopotentials, with optional partial core corrections for exchange and correlation and for the core kinetic energy density.
- Kleinman-Bylander form of the pseudopotential, with a ghost-state test.
- Tests of the pseudopotentials: ionization and excitation energies compared with all-electron calculations, for the semi-local and Kleinman-Bylander forms, and logarithmic derivatives.
- Atomic orbital basis sets (SZ, DZ, DZP or bespoke) for the cpw2000 code.
- Pseudopotential files for several plane-wave and real-space codes:

| file | code |
|------|------|
| *Xy*`_POTKB_F.DAT` | cpw2000 (with the atomic basis sets) |
| *Xy*`_LDA_TM.UPF` (or `_GGA_`) | Quantum ESPRESSO (UPF, with spin-orbit) |
| *Xy*`.psp8` | ABINIT (with spin-orbit) |
| *Xy*`.psml` | SIESTA (PSML 1.1) |
| *Xy*`.psf` | SIESTA |
| *Xy*`_POTRE.DAT` | PARSEC |

where *Xy* is the chemical symbol.

## Quick start

The code is self-contained, it only needs a Fortran compiler.  `gnuplot` is used for the plots if it is installed.

    cd Src
    make                      # choose the compiler in Src/make.inc
    cd ..
    mkdir WORK; cd WORK
    ../Src/atom_all_gfortran.exe Si        # the suffix is the compiler of make.inc, in this case gfortran

`atom_all` writes a default input file `atom.dat` for the element, generates the pseudopotential, tests it and writes all the files above.  If the element has tested pseudopotentials (for now only carbon, `HARD` and `SOFT`) the default is the first one, otherwise the parameters are an educated guess from the old tables, and the program says which.  Optional further arguments, in any order, choose a tested pseudopotential (`HARD`, `SOFT`, or `NONE` for the old tables), the partial core correction (yes/no) and the basis set (SZ/DZ/DZP).  The same steps can be run one at a time (`atom`, `psd_gen`, `kb_conv`, `atom_test`, `kb_test` and the plotting programs) to change the parameters interactively.

## Documentation

- `Doc/UserGuide/user_guide.pdf`: installation, running the code, the `atom.dat` input file and computational details.
- `Doc/Other_docs/files_and_units.md`: the files used by the programs.
- Source code documentation can be generated with doxygen (`Doc/doxy/Doxyfile`).

## Validation and tools

- `Validation/run.sh` runs `atom_all` for all the elements of `Validation/list_elements` and collects the test results in `ptb.out`.
- `Tools/`: python scripts to plot the results with matplotlib (`atom_plot_mpl.py`, `plot_basis.py`) and to check the PSML and psp8 files (`check_psml.py`, `check_psp8.py`).

## What's new in 6.2.0 (since 6.0.8)

- New output formats: UPF for Quantum ESPRESSO with spin-orbit, psp8 for ABINIT (with spin-orbit) and PSML for SIESTA; the PARSEC and SIESTA `.psf` outputs were repaired.
- Tested pseudopotentials: an element can have several pseudopotentials that were compared with all-electron calculations (for now `HARD` and `SOFT`), chosen by a label.  The programs say whether the pseudopotential generated is one of them or an educated guess, and the label goes to the titles of the output files (`var=`).  This will be updated slowly, started with carbon.
- Core kinetic energy density and its partial core correction (E. Gareis, I. Schelter and S. Kümmel, Phys. Rev. B 114, 165112 (2026)).
- Atomic orbital basis sets for cpw2000: several basis sets per element, confined orbitals with cutoff radii derived from the SIESTA energy-shift radius, the relaxed energy derivative as second zeta, and a single perturbative polarization orbital.  Bespoke basis sets in the interactive mode of `kb_conv`.
- OpenMP parallelization of the Bessel transforms.
- Plots with matplotlib, and the gnuplot data files are written even if gnuplot is not installed.
- Numerical improvements: removal of the noise of GGA potentials near the origin, values at the origin from symmetric points, full precision in the intermediate file `psd.pot`, consistent file units.
- Test configurations without a self-consistent solution are skipped instead of stopping the program.
- The core radii line of `atom.dat` is back to the original column convention.  In `atom.dat` files written by
  versions 6.0.x the partial core numbers must be moved ten columns to the left, otherwise they are
  read in the wrong fields without a warning (see the user guide).
- An improved DRAFT user guide, and a validation script for all the elements.
- Several bug fixes.

## References

Early references of the code are:

Nonlinear ionic pseudopotentials in spin-density-functional calculations, Steven G. Louie, Sverre Froyen, and Marvin L. Cohen, Phys. Rev. B 26, 1738 (1982)

Efficient pseudopotentials for plane-wave calculations, N. Troullier and José Luís Martins, Phys. Rev. B 43, 1993 (1991).

## License

GNU General Public License v2 (see `LICENSE`).
