# Compiler settings
#===================

F90 = gfortran
F90LINKER = $(F90)

# Default build type is release
FFLAGS = -O3 -pipe #Default release flags

# Debugging flags
DEBUG_FLAGS = -Og -g -pipe -Wall -Wextra -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow

# Source files (mix of .f and .f90 files)
SRCS_F90 = \
    vert_diff_exp_mod.f90 \
    vert_adv_exp_McCormack_mod.f90 \
    lxmie_mod.f90 \
    gamma_mod.f90 \
    mini_cloud_opac_mie_mod.f90 \
    mini_cloud_vf_mod.f90 \
    mini_cloud_3_gamma_mod.f90 \
    test_main.f90

SRCS_F = dlsode.f

# Object files
OBJECTS = $(SRCS_F90:.f90=.o) $(SRCS_F:.f=.o)

# Executable
EXECS = ../mini_cloud_3_gamma

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
