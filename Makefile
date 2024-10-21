# Fortran Compiler
FC = gfortran

# Executable
MAIN = main

# Use nf-config to get the necessary flags for NetCDF-Fortran
FFLAGS_NC = $(shell nf-config --fflags)
FLIBS_NC = $(shell nf-config --flibs)
# FFLAGS_NC =
# FLIBS_NC =
FSTDLIB_PATH = /home/lorenzolds/learn-fortran/stdlib/build/src
# Compiler and Linker flags
CFLAGS = -O3 -Wall -Wextra -std=f2008
CFLAGS += -fPIC -fcheck=all -fbacktrace -fmax-errors=3
CFLAGS += -g $(FFLAGS_NC) -I$(FSTDLIB_PATH)/mod_files #  flags for compiling
LFLAGS = $(FLIBS_NC) -L$(FSTDLIB_PATH)
LFLAGS += -lfortran_stdlib

# Source and object files
#TODO compile all .f90 files in cwd
SRC = fortran202x_split config_module utils_module netcdf_utils main
OBJ = $(SRC:=.o)
# Default target
default: $(MAIN)

# Linking the executable
$(MAIN): $(OBJ)
	$(FC) $(CFLAGS) -o $(MAIN) $(OBJ) $(LFLAGS)

# Pattern rule to compile .f90 to .o
%.o: %.f90
	$(FC) $(CFLAGS) -c $< -o $@

# Clean target to remove compiled files
clean:
	rm -f $(OBJ) *.mod

clean_main:
	rm -f $(OBJ) $(MAIN) *.mod
