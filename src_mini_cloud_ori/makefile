# Compiler settings
#===================

F90 = gfortran
F90LINKER = $(F90)

# Default build type is release
FFLAGS = -O3 -pipe  # Default release flags

# Debugging flags
DEBUG_FLAGS = -Og -g -pipe -Wall -Wextra -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow

# Source files (mix of .f and .f90 files)
SRCS_F90 = \
    lxmie_mod.f90 \
    mini_cloud_precision.f90 \
    mini_cloud_class.f90 \
    mini_cloud_saturation.f90 \
    mini_cloud_nucleation.f90 \
    mini_cloud_chem.f90 \
    mini_cloud_rhs.f90 \
    mini_cloud_i_dlsode.f90 \
    mini_cloud_vf.f90 \
    mini_cloud_opac_mie_mod.f90 \
    main.f90

SRCS_F = dlsode.f

# Object files
OBJECTS = $(SRCS_F90:.f90=.o) $(SRCS_F:.f=.o)

# Executable
EXECS = ../mini_cloud_ori

#====================
# Phony targets
.PHONY: all clean debug release

# Default target (release by default)
all: $(EXECS)

# Debug build (overrides FFLAGS with debug flags)
debug: FFLAGS = $(DEBUG_FLAGS)
debug: clean all

# Release build (uses default release flags)
release: clean all

# Linking the executable
$(EXECS): $(OBJECTS)
	$(F90LINKER) $(FFLAGS) $(OBJECTS) -o $(EXECS)

# Pattern rule for compiling Fortran 90 files
%.o: %.f90
	$(F90) $(FFLAGS) -c $<

# Pattern rule for compiling Fortran 77 (.f) files
%.o: %.f
	$(F90) $(FFLAGS) -c $<

# Clean target
clean:
	rm -f *.o *.mod *~ *__genmod.f90 $(EXECS)
