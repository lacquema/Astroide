#! /bin/sh
#
# Makefile to build simulations
#

SHELL = /bin/sh
ARCH = $(shell uname)
DEBUGGING=NO
PARALLEL=YES
DEBUG = -g
FFLAGS = -O3
LDFLAGS = -O3

# dahu
ifeq ($(ARCH),Linux)
SWIFT_DIR = /home/lacquema/bin
SWIFT_DIR_LIB = /home/lacquema/lib
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
FC = gfortran
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

swift_hjs_tid_corr: swift_hjs_tid_corr.o
swift_hjs_corr: swift_hjs_corr.o
swift_mvs: swift_mvs.o
swift_mig: swift_mig.o
swift_mig_fix: swift_mig_fix.o
swift_mvs_par: swift_mvs_par.o
swift_mvs_mig: swift_mvs_mig.o
swift_mvs_mig_fd: swift_mvs_mig_fd.o
swift_mvs_mig_fd_2: swift_mvs_mig_fd_2.o
swift_mvs_fd: swift_mvs_fd.o
swift_mjs: swift_mjs.o
swift_mvs_coll: swift_mvs_coll.o
swift_rmvs3: swift_rmvs3.o
swift_rmvs3_corr: swift_rmvs3_corr.o
swift_rmvs3_par: swift_rmvs3_par.o
swift_rmvs4: swift_rmvs4.o
swift_hjs: swift_hjs.o
swift_hjs_hyp: swift_hjs_hyp.o
swift_hjs4: swift_hjs4.o
swift_hjs_tid: swift_hjs_tid.o
swift_av: swift_av.o
swift_symba5: swift_symba5.o
swift_symba5_dup: swift_symba5_dup.o
read: read.o




