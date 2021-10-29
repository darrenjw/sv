
# sv

## ANSI 'C' library for stochastic volatility models

**Please note that this is very old and unmaintained code. It is archived here mainly for historical interest.**

C Library for Bayesian (MCMC) analysis of Stochastic Volatility models.
C structs "isv" (for univariate SV) and "fsv" (for multivariate factor
SV) are declared in "isv.h" and "fsv.h" respectively. A (totally
inadequate) test can be found in "isv_test.c". Examples of use can be
found in "isvprog.c" and "fsvprog.c". The theory and notation are as
explained in:

* Pitt, M.K. & N. Shephard (1999), Time-varying covariances: a factor
stochastic volatility approach, Bayesian Statistics 6, OUP, 547-570.

Once you've read and understood this paper, you should be able to read
the source code for the "classes" ("isv.c" and "fsv.c"), together with
that of the examples, in order to figure out exactly what is going on.

There are however, one or two other things that need to be explained in
order to understand what is going on. Firstly, "make" is used to build
and run all of the programs --- you will need to study the Makefile
carefully in order to figure out what is going on, and may have to
"tweak" it slightly to sort it out for your system.

Type "make isvprog" to build the code, and "make isv.tab" to run it.
Type "make gendata_isv" to build the data simulation code and "make
isvdata.dat" to simulate some data. However, because make handles all
dependencies for you, just typing "make" will build and run everything
in the correct order. Typing "make test" will build and run a few tests.
Similarly, typing "make fsv.tab" will build and run a simple FSV
example.

The code relies heavily on the GNU Scientific Library (GSL) as well as
the standard C library. However, there will be easy-to-install binary
packages for the GSL for your system, so this shouldn't cause any
problems. The FSV part of the library makes quite a bit of use of BLAS.
By default, the Makefile uses the BLAS which ships with the GSL, so
performance could be improved slightly by instead using a BLAS optimised
for your system. That said, I don't think it will make much difference,
so I personally wouldn't bother...

The code files have the following purposes:

* isv.h Declaration of ISV type struct and function prototypes 
* isv.c Functions and procedures making up the ISV library 
* isvprog.c Example
main executable program (for the ISV stuff) 
* gendata_isv.c Program to
simulate some data (for the ISV stuff) 
* test_isv.c Some not very good
test code... ;-) 
* fsv.h Declaration of FSV type struct and function
prototypes
* fsv.c Functions and procedures making up the FSV library
* fsvprog.c Example main executable program (for the FSV stuff)
* gendata_fsv.c Program to simulate some data (for the FSV stuff)

All code is Copyright (C) 2004 Darren Wilkinson but is free in the sense
of the GNU General Public License.

Historic note
-------------

This library is a port of some classes I wrote for the Sather
programming language early in 2001. I've tried to mimic them in the
usual C way. I did the port of the ISV class at the end of 2002, as part
of the case study for my book chapter on "Parallel Bayesian
computation". The port of FSV was done mid 2004, when I decided that
Sather was dead and buried...

Last updated: 24/6/04

Darren Wilkinson

