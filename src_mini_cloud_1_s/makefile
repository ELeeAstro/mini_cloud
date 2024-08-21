# for gfortran Compiler
#===================

F90 = gfortran
F90LINKER =  gfortran

#Debugging and development flags
#FFLAGS	= -Og -g -pipe -Wall -Wextra -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow
#FFLAGS	= -Og -pipe -Wall -Wextra -g -fbacktrace

#Serial flags
FFLAGS  = -O3 -pipe

# for intel Compiler
#====================

#F90 = ifx
#F90LINKER = ifx

#FFLAGS   = -O0 -g -traceback -fpp -prec-div -fp-model source -fpe0 -ipo
#FFLAGS   = -O3 -fpp -fp-model source -ipo

# for nvfortran compiler
#====================

#F90 = nvfortran
#F90LINKER = nvfortran

#FFLAGS    = -O1

#====================

DEFS      =
INCLUDES  =
LFLAGS    = $(FFLAGS)

OBJECTS = \
dopri5.o \
mini_cloud_1_s_mod.o \
test_main.o \

# executable statement
EXECS  = ../mini_cloud_1_s

.SUFFIXES : .o .f90 .f .F90

default: mini_chem_1_s

mini_chem_1_s: $(OBJECTS)
	$(F90LINKER) $(LFLAGS) $(OBJECTS) -o $(EXECS) $(LAPACK) $(BLAS)

clean:
	rm -f *.o *.mod *~ *__genmod.f90 $(EXECS)

.f90.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<

.F90.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<

.f.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<
