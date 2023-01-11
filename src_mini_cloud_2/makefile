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
#FFLAGS	= -O3 -pipe -fopenmp

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
decsol.o \
dc_decsol.o \
seulex.o \
rodas.o \
radau5.o \
ros4.o \
mini_dihrt_precision.o \
mini_dihrt_class.o \
mini_dihrt_saturation.o \
mini_dihrt_nucleation.o \
mini_dihrt_chem.o \
mini_dihrt_rhs.o \
mini_dihrt_i_seulex.o \
mini_dihrt_i_rodas.o \
mini_dihrt_i_radau5.o \
mini_dihrt_i_ros4.o \
main.o

# executable statement
EXECS  = ../mini_dihrt

.SUFFIXES : .o .f90 .f

default: mini_dihrt

mini_dihrt: $(OBJECTS)
	$(F90LINKER) $(LFLAGS) $(OBJECTS) -o $(EXECS)

clean:
	rm -f *.o *.mod *~ *__genmod.f90 $(EXECS)

.f90.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<

.f.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<
