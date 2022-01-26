ROOT = ..

include $(ROOT)/make.inc

src77 = $(wildcard *.f)
src90 = $(wildcard *.f90)



obj =  $(src77:.f=.o) $(src90:.f90=.o)

libadd:    $(obj)
	ar -rv $(ROOT)/libatom$(SUFFIX).a *.o

%.o: %.f
	$(F77) -c $(FCOPTS) $(FCEXTRA) $<

%.o: %.f90
	$(F77) -c $(FCOPTS) $(FCEXTRA) $<

clean:
	rm -f *.o *.mod *genmod*

.PHONY:   clean libadd
