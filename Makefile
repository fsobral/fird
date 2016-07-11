# Paths

FKSS_HOME = $(CURDIR)
SRC = $(FKSS_HOME)/src
LIB = $(FKSS_HOME)/lib
TES = $(FKSS_HOME)/tests
SOL = $(FKSS_HOME)/solvers
BIN = $(FKSS_HOME)/bin
OBJ = $(FKSS_HOME)/objects
DOC = $(FKSS_HOME)/doc

ifndef PROBLEM
   PROBLEM = tests/examples/hs14.f90
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

# Generate the main FKSS library
lib: base
	$(MAKE) -C $(SOL) install
	$(MAKE) -C $(SRC) all install

# User-defined executable
fkss: all
	$(FC) $(foreach i,$(SOLVERLIB) $(LIB),-L$(i) ) \
	$(FCC) $(PROBLEM) $(SOL)/*.o $(LOPTS) -o $(BIN)/$@

# Hock-Schittkowski test set executable
hstests: all
	$(MAKE) -C $(TES) hs

	$(FC) $(foreach i,$(SOLVERLIB) $(LIB),-L$(i) ) \
	$(FCC) tests/hs/hstests.f -lhs $(SOL)/*.o $(LOPTS) -o $(BIN)/$@

# Documentation
doc:
	$(MAKE) -C $(DOC) all

clean:
	rm -vfr *~ $(LIB)/* $(BIN)/* $(OBJ)/*
	$(foreach i,$(SRC) $(SOL) $(TES) $(DOC),$(MAKE) -C $(i) clean;)

.PHONY: lib all clean doc
