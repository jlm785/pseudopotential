
include make.inc

#  The main programs are compiled without OpenMP, the executables are linked
#  with FCOMP (OpenMP runtime used by lib_num/ft_radial.f90 and ft_precise.f90)

all: atom.exe  atom_all.exe  psd_gen.exe  kb_conv.exe  plot_show.exe  plot_kb_show.exe  plot_ln_show.exe atom_test.exe kb_test.exe

atom.exe:  libatom  atom.f90
	$(F77) -c $(FCOPTS) atom.f90
	$(F77) -o atom$(SUFFIX).exe $(FCOPTS) $(FCOMP) atom.o libatom$(SUFFIX).a

atom_all.exe:  libatom  atom_all.f90
	$(F77) -c $(FCOPTS) atom_all.f90
	$(F77) -o atom_all$(SUFFIX).exe $(FCOPTS) $(FCOMP) atom_all.o libatom$(SUFFIX).a

psd_gen.exe:  libatom  psd_gen.f90
	$(F77) -c $(FCOPTS) psd_gen.f90
	$(F77) -o psd_gen$(SUFFIX).exe $(FCOPTS) $(FCOMP) psd_gen.o libatom$(SUFFIX).a

kb_conv.exe:  libatom  kb_conv.f90
	$(F77) -c $(FCOPTS) kb_conv.f90
	$(F77) -o kb_conv$(SUFFIX).exe $(FCOPTS) $(FCOMP) kb_conv.o libatom$(SUFFIX).a

plot_show.exe:  libatom  plot_show.f90
	$(F77) -c $(FCOPTS) plot_show.f90
	$(F77) -o plot_show$(SUFFIX).exe $(FCOPTS) $(FCOMP) plot_show.o libatom$(SUFFIX).a

plot_kb_show.exe:  libatom  plot_kb_show.f90
	$(F77) -c $(FCOPTS) plot_kb_show.f90
	$(F77) -o plot_kb_show$(SUFFIX).exe $(FCOPTS) $(FCOMP) plot_kb_show.o libatom$(SUFFIX).a

plot_ln_show.exe:  libatom  plot_ln_show.f90
	$(F77) -c $(FCOPTS) plot_ln_show.f90
	$(F77) -o plot_ln_show$(SUFFIX).exe $(FCOPTS) $(FCOMP) plot_ln_show.o libatom$(SUFFIX).a

atom_test.exe:  libatom  atom_test.f90
	$(F77) -c $(FCOPTS) atom_test.f90
	$(F77) -o atom_test$(SUFFIX).exe $(FCOPTS) $(FCOMP) atom_test.o libatom$(SUFFIX).a

kb_test.exe:  libatom  kb_test.f90
	$(F77) -c $(FCOPTS) kb_test.f90
	$(F77) -o kb_test$(SUFFIX).exe $(FCOPTS) $(FCOMP) kb_test.o libatom$(SUFFIX).a






libatom:
	(cd lib_atom          &&  make libadd )
	(cd lib_num           &&  make libadd )
	(cd lib_blas_pack     &&  make libadd )
	(cd lib_psd           &&  make libadd )
	(cd lib_kb            &&  make libadd )
	(cd lib_plot          &&  make libadd )



cleanobj:
	(cd lib_atom          &&  make clean )
	(cd lib_num           &&  make clean )
	(cd lib_blas_pack     &&  make clean )
	(cd lib_psd           &&  make clean )
	(cd lib_kb            &&  make clean )
	(cd lib_plot          &&  make clean )
	(/bin/rm -f *.mod *.o *~)

clean:  cleanobj
	(/bin/rm -f *.a *.exe)

.PHONY:   clean cleanobj libatom
