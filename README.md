FIRD - Filter Inexact Restoration Derivative-free algorithm
===========================================================

This algorithm is based on the working paper of [Ferreira, Karas,
Sachine,
Sobral](http://www.optimization-online.org/DB_FILE/2015/07/5024.pdf)
(*Submitted*).

`FIRD` is able to solve general derivative-free constrained
optimization problems in which the derivatives of the objective
function are not available.

How to install
--------------

In the provided version, `FIRD` uses a modified version of the
algorithm [`TRDF`](https://github.com/fsobral/trdf) for building and
updating the trust region models. To solve the nonlinear trust-region
problems, it uses [`ALGENCAN`](http://www.ime.usp.br/~egbirgin/tango).

### Installing `ALGENCAN`

To install and build `ALGENCAN` library

1. [Download `ALGENCAN`](http://www.ime.usp.br/~egbirgin/tango);

1. Unpack it and, in the main folder, e.g `$ALGSRC`, simply type
`make` in order to generate `$ALGSRC/lib/libalgencan.a` file;

### Installing and building `FIRD`

1. Download and unpack `FIRD`;

1. In the main subdirectory, edit file `Makefile` with the path to
`ALGENCAN` library

        SOLVERLIB = $ALGSRC/lib

    where `$ALGSRC` is the name of the full path where `ALGENCAN` has been
installed;

1. Type

        make fird

    to generate test problem HS14 from the Hock-Schittkowski
    collection. The executable file will be generated in `/bin`
    subdirectory. Or, type


        make fird PROBLEM=tests/examples/toyprob.f90

    for a small example.

1. Run it

        ./bin/fird


Advanced users
--------------

`FIRD` can be modified in order to use different solvers for

* Feasibility phase
* Construction and updating models (optimality phase)
* Solving nonlinear programming trust-region subproblems (optimality phase)

The manual for advanced information can be generated with the command

    make doc