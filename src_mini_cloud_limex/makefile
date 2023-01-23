# for gfortran Compiler
#===================

F90 = gfortran
F90LINKER =  gfortran

#Debugging and development flags
#FFLAGS	= -Og -g -pipe -Wall -Wextra -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow
#FFLAGS	= -Og -pipe -Wall -Wextra -g -fbacktrace

#Serial flags
FFLAGS  = -O2 -pipe

#Parallel flags
#FFLAGS	= -O2 -pipe -fopenmp

# for ifort Compiler
#====================

#F90 = ifort
#F90LINKER = ifort

#FFLAGS   = -O0 -g -traceback -fpp -prec-div -fp-model source -fpe0 -ipo
#FFLAGS    = -O0 -g -traceback -xHost -fpp -fp-model source -qopenmp -ipo
#FFLAGS    = -O3 -xHost -qopenmp -fpp -fp-model source -ipo

#====================

DEFS      =
INCLUDES  =
LFLAGS    = $(FFLAGS)


OBJECTS = \
limex.o \
limex_aux.o \
limex_dense.o \
mini_cloud_precision.o \
mini_cloud_class.o \
mini_cloud_saturation.o \
mini_cloud_nucleation.o \
mini_cloud_chem.o \
mini_cloud_rhs.o \
mini_cloud_i_limex.o \
mini_cloud_vf.o \
main.o

# executable statement
EXECS  = ../mini_cloud

.SUFFIXES : .o .f90 .f

default: mini_cloud

mini_cloud: $(OBJECTS)
	$(F90LINKER) $(LFLAGS) $(OBJECTS) -o $(EXECS)

clean:
	rm -f *.o *.mod *~ *__genmod.f90 $(EXECS)

.f90.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<

.f.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<
