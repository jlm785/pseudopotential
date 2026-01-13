#
#   definition of the variables for the Makefile that makes the pw library.
#

#   If you are lucky the defaults for each type of CPU may work and give you a
#   fast functionning code.  But as that depends on the compilers and libraries
#   installed on your machine, nothing is assured.
#
#   You will find below alternatives that may work by switching commented
#   with uncommented lines.

#   Indicate what type of cpu you are using (uncomment one).  It will only suggest the best
#   compiler for the job.  "cat /proc/cpuinfo" will tell you what CPU is on your machine.

CPU=Intel
# CPU=AMD
# CPU=ARM


#   Suggestions for compilers.  Feel free to override the suggestions.

ifeq ($(CPU),Intel)
   F77=ifx
   # F77=ifort
   # F77=gfortran
   # F77=pgfortran
else
   ifeq ($(CPU),AMD)
      F77=gfortran
      # F77=ifort
      # F77=pgfortran
   else
      F77=gfortran
   endif
endif


#   below this point you should have a good reason (fine tune debugging, curiosity) to change
#   the defaults.  Each compiler has its own section.


#   debugging option (true, false).  Use if the code throws a segmentation fault
#   Uncomment also optional changes for debugging below for fine tuning...

DEBUG=false
# DEBUG=true



#   suffix for executables if you want to compare compilers

SUFFIX:=_$(F77)
# SUFFIX=

ifeq ($(DEBUG),true)
   SUFFIX:=$(SUFFIX)_debug
endif



#  For Intel ifx compiler:

ifeq ($(F77),ifx)

   #   profiling option (-g,empty).  Use it for tunning the code.

   PROF=
   # PROF=-g

   ifeq ($(DEBUG),false)
      FCOPTS=-xHost -O3 $(PROF)
   else
      FCOPTS=-check all -warn all,nodec,interfaces,noexternal -gen-interfaces -traceback -fpe0 -g -O0
   endif

endif



#  For Intel ifort compiler:

ifeq ($(F77),ifort)

   #   profiling option (-g,empty).  Use it for tunning the code.

   PROF=
   # PROF=-g

   ifeq ($(DEBUG),false)
      FCOPTS=-xHost -O3 $(PROF) -diag-disable=10448
   else
      FCOPTS=-check all -warn all,nodec,interfaces,noexternal -gen-interfaces -traceback -fpe0 -fp-stack-check -g -O0 -diag-disable=10448
   endif

endif




#  For gcc gfortran compiler:


ifeq ($(F77),gfortran)

   #   profiling option (-pg,empty).  Use it for tunning the code.

   PROF=
   # PROF=-pg

   #   Compiling options and libraries.  Must have the chosen library installed.

   ifeq ($(DEBUG),false)

      FCOPTS=-fomit-frame-pointer -O3 -mtune=native $(PROF)

   else

      FCOPTS=-g -Wall -fbounds-check -pedantic-errors -ffpe-trap=zero,invalid,overflow

   endif

endif




#  For portland group pgfortran compiler, it is much less used than the previous, so beware.


ifeq ($(F77),pgfortran)

   #   Never tried profiling or debugging

   PROF=


   #   Compiling options and libraries.  Must have the chosen library installed.

   FCOPTS=-fast

endif


