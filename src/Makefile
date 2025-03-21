# Program name and version
PROGRAM = settling
VERSION = 1.0

# Directories
INST_DIR = $(HOME)/bin
SRC_DIR = .
BIN_DIR = ../bin
BUILD_DIR = ../build

# Compiler and flags
FC = gfortran


# Debug version    
#OPT = -g -fbounds-check -Wall -fcheck=all
#OPDE= d
# Optimized version
OPT = -O3
OPDE= o 


FCOPTS = -ffixed-line-length-132 -fcray-pointer -cpp -I../inc
FLAGS = -J$(BUILD_DIR) # Store .mod files in build directory

# Source files
SMODS = mod_common.f90 mod_rotation.f90 mod_nn.f90 mod_dragForce.f90

SFIX =

SF90 = main.f90 lib_sim.f90

# Object files
MODOBJ = $(patsubst %.f90,$(BUILD_DIR)/%.o,$(SMODS))
FIXOBJ = $(patsubst %.f,$(BUILD_DIR)/%.o,$(SFIX))
F90OBJ = $(patsubst %.f90,$(BUILD_DIR)/%.o,$(SF90))

# Libraries
#LDOPTO = -lblas -llis
LDOPTO = 


# Default target
all: $(BUILD_DIR) program

# Rule to create the build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Linking the final executable
program: $(MODOBJ) $(F90OBJ) $(FIXOBJ)
	$(FC) $(FLAGS) $(FCOPTS) $(OPT) -o $(BIN_DIR)/$(PROGRAM) $^ $(LDOPTO)

# Pattern rule to compile Fortran source files into object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FLAGS) $(FCOPTS) $(OPT) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f
	$(FC) $(FLAGS) $(FCOPTS) $(OPT) -c $< -o $@

# Clean up generated files
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)/$(PROGRAM)


install:
	cp $(BIN_DIR)/$(PROGRAM) $(INST_DIR)/$(PROGRAM)-$(VERSION)-$(OPDE)	