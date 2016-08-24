# Paths

FIRD_HOME = $(CURDIR)
SRC = $(FIRD_HOME)/src
LIB = $(FIRD_HOME)/lib
TES = $(FIRD_HOME)/tests
SOL = $(FIRD_HOME)/solvers
BIN = $(FIRD_HOME)/bin
OBJ = $(FIRD_HOME)/objects
DOC = $(FIRD_HOME)/doc

ifndef PROBLEM
   PROBLEM = tests/examples/hs14.f90
endif

ifndef PPROBLEM
   PPROBLEM = toyprob.f90
endif

# Compiler options

CC  = gcc
FC  = gfortran
FCC = 

# Solver configuration parameters

SOLVERLIB = /opt/tango/algencan-3.0.0/lib
RESTORATION_INTERFACE = restoration
OPTIMIZATION_INTERFACE = optimization

# The following order do matters!
SLOPTS = -ldfoirfilter -lalgencan 

# Linking options

LOPTS = $(SLOPTS)

export

all: lib
	mkdir -p $(BIN)

base:
	mkdir -p $(LIB)
	$(MAKE) -C $(SRC) base

# Generate the main FIRD library
lib: base
	$(MAKE) -C $(SOL) install
	$(MAKE) -C $(SRC) all install

# User-defined executable
fird: all
	$(FC) $(foreach i,$(SOLVERLIB) $(LIB),-L$(i) ) \
	$(FCC) $(PROBLEM) $(SOL)/*.o $(LOPTS) -o $(BIN)/$@

# Hock-Schittkowski test set executable
hstests: all
	$(MAKE) -C $(TES) hs

	$(FC) $(foreach i,$(SOLVERLIB) $(LIB),-L$(i) ) \
	$(FCC) tests/hs/hstests.f -lhs $(SOL)/*.o $(LOPTS) -o $(BIN)/$@

# Probability test
# To solve other probability problems using Genz's MVNDST subroutine,
# please add the new problem in the tests/probability directory and
# set the PPROBLEM variable accordingly
probability: all
	$(MAKE) -C $(TES) probability

	$(FC) $(foreach i,$(SOLVERLIB) $(LIB),-L$(i) ) \
	$(FCC) tests/probability/$(PPROBLEM) -lprobability   \
	$(SOL)/*.o $(LOPTS) -o $(BIN)/$@

# Documentation
doc:
	$(MAKE) -C $(DOC) all

clean:
	rm -vfr *~ $(LIB)/* $(BIN)/* $(OBJ)/*
	$(foreach i,$(SRC) $(SOL) $(TES) $(DOC),$(MAKE) -C $(i) clean;)

.PHONY: lib all clean doc
