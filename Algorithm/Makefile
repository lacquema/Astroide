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


BIN_DIR = ./bin

ifeq ($(DEBUGGING),YES)
FFLAGS += $(DEBUG)
LDFLAGS = $(DEBUG)
endif
ifeq ($(PARALLEL),YES)
FFLAGS += -fopenmp
LDFLAGS += -fopenmp
endif
FC = gfortran-mp-13

%.o : %.f
	$(FC) -c $(FFLAGS) $< -o $@

% : %.o
ifeq ($(PARALLEL),YES)
	$(FC) $(LDFLAGS) -o $(BIN_DIR)/$(join $@,p) $?
else
	$(FC) $(LDFLAGS) -o $(BIN_DIR)/$@ $?
endif
	rm $?

clean:
	rm -f *.o *~ core

# Simulation

COORD = swift/coord/
IO = swift/io
ORBEL = swift/orbel
UTIL = swift/util

gen_tout_multi: ${COORD}/coord_b2h.o ${IO}/io_dump_pl.o ${IO}/io_dump_tp.o ${IO}/io_init_param_hb.o ${IO}/io_open.o ${ORBEL}/orbel_el2xv.o ${UTIL}/util_hills.o gen_tout_multi.o
# gen-tout_multi_corr:gen_tout_multi_corr.o
# mbodies_multi: mbodies_multi.o
# followbodies: followbodies.o
# visee_multi: visee_multi.o





