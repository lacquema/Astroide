#! /bin/sh
#
# Makefile to build simulations
#

SHELL = /bin/sh
ARCH = $(shell uname)
DEBUGGING=NO
PARALLEL=NO
DEBUG = -g
FFLAGS = -O3
LDFLAGS = -O3

# dahu
ifeq ($(ARCH),Linux)
SWIFT_DIR = ./swift
SWIFT_DIR_LIB = ./lib
ifeq ($(DEBUGGING),YES)
FFLAGS += $(DEBUG)
LDFLAGS = $(DEBUG)
endif
ifeq ($(PARALLEL),YES)
FFLAGS += -fopenmp
LDFLAGS += -fopenmp
LOADLIBES = -L$(SWIFT_DIR_LIB) -lswiftp
else
LOADLIBES = -L$(SWIFT_DIR_LIB) -lswift
endif
FC = gfortran-mp-13
endif

%.o : %.f
	$(FC) -c $(FFLAGS) $< -o $@

% : %.o
ifeq ($(PARALLEL),YES)
	$(FC) $(LDFLAGS) -o $(SWIFT_DIR)/$(join $@,p) $? $(LOADLIBES)
else
	$(FC) $(LDFLAGS) -o $(SWIFT_DIR)/$@ $? $(LOADLIBES)
endif
	rm $?

clean:
	rm -f *.o *~ core

# Simulation

gen_tout_multi:gen_tout_multi.o
# gen-tout_multi_corr:gen_tout_multi_corr.o
# mbodies_multi: mbodies_multi.o
# followbodies: followbodies.o
# visee_multi: visee_multi.o





