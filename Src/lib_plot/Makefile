ROOT = ..

include $(ROOT)/make.inc

src90 = $(wildcard *.f90)



obj =  $(src90:.f90=.o)

libadd:    $(obj)
	ar -rv $(ROOT)/libatom$(SUFFIX).a *.o

%.o: %.f90
	$(F77) -c $(FCOPTS) $<

clean:
	rm -f *.o *.mod *genmod*

.PHONY:   clean libadd
