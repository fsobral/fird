# Paths

FKSS_HOME = $(CURDIR)
SRC = $(FKSS_HOME)/src
LIB = $(FKSS_HOME)/lib
TES = $(FKSS_HOME)/tests
SOL = $(FKSS_HOME)/solvers
BIN = $(FKSS_HOME)/bin
OBJ = $(FKSS_HOME)/objects

ifndef PROBLEM
   PROBLEM = tests/examples/hs14.f90
endif

# Compiler options

CC  = gcc
FC  = gfortran
FCC = 

# Solver configuration parameters

SOLVERLIB = /opt/tango/algencan-3.0.0/lib $(FKSS_HOME)/../trdf/lib
MODULES=$(FKSS_HOME)/../trdf/src
RESTORATION_INTERFACE = restoration
OPTIMIZATION_INTERFACE = trdf_solver_interface

# The following order do matters!
SLOPTS = -ldfoirfilter -ltrdf -lalgencan 

# Linking options

LOPTS = $(SLOPTS)

export

all: lib
	mkdir -p $(BIN)

base:
	$(MAKE) -C $(SRC) base

# Generate the main FKSS library
lib: base
	mkdir -p $(LIB)
	$(MAKE) -C $(SOL) install
	$(MAKE) -C $(SRC) all install

# Generate the solver interface object
#solver:
#	$(MAKE) -C $(SOL) install

# User-defined executable
fkss: all
	$(FC) $(foreach i,$(SOLVERLIB) $(LIB),-L$(i) ) \
	$(FCC) $(PROBLEM) $(SOL)/*.o $(LOPTS) -o $(BIN)/$@

# User-defined C executable
c_trdf: all
	$(CC) -L$(SOLVERLIB) -L$(LIB) $(OBJ)/solver.o \
	$(PROBLEM) $(LOPTS) -lgfortran -lm -o $(BIN)/$@

# Hock-Schittkowski test set executable
hstests: all
	$(MAKE) -C $(TES) hs

	$(FC) -L$(SOLVERLIB) -L$(LIB) $(OBJ)/solver.o \
	$(FCC) -I$(SRC) tests/hs/hstests.f -lhs $(LOPTS) \
	-o $(BIN)/$@ 

clean:
	rm -vf *~ $(LIB)/* $(BIN)/* $(OBJ)/*
	$(foreach i,$(SRC) $(SOL),$(MAKE) -C $(i) clean;)

.PHONY: lib all clean
