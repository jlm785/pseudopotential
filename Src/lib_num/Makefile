ROOT = ..

include $(ROOT)/make.inc

src90 = $(wildcard *.f90)



obj =  $(src90:.f90=.o)

libadd:    $(obj)
	ar -rv $(ROOT)/libatom$(SUFFIX).a *.o

%.o: %.f90
	$(F77) -c $(FCOPTS) $<

#  OpenMP only for the files with OpenMP directives and the function
#  called in their parallel regions (FCOMP is defined in make.inc)

ft_radial.o: ft_radial.f90
	$(F77) -c $(FCOPTS) $(FCOMP) $<

ft_precise.o: ft_precise.f90
	$(F77) -c $(FCOPTS) $(FCOMP) $<

sbessj.o: sbessj.f90
	$(F77) -c $(FCOPTS) $(FCOMP) $<

clean:
	rm -f *.o *.mod *genmod*

.PHONY:   clean libadd
