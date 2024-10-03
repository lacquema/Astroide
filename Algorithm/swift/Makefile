#! /bin/sh
#
# Makefile to build simulations
#

SHELL = /bin/sh
ARCH = $(shell uname)
DEBUGGING=NO
PARALLEL=YES
DEBUG = -g
FFLAGS = -c 
OPTIMISATION = -O3

# dahu
ifeq ($(ARCH),Linux)
SWIFT_DIR = .
SWIFT_DIR_LIB = ../lib
ifeq ($(DEBUGGING),YES)
FFLAGS += $(DEBUG)
endif
ifeq ($(PARALLEL),YES)
FFLAGS += -fopenmp
endif
FFLAGS += $(OPTIMISATION)
FC = gfortran
endif


COLL = $(SWIFT_DIR)/coll
ANAL = $(SWIFT_DIR)/anal
BS = $(SWIFT_DIR)/bs
COORD = $(SWIFT_DIR)/coord
DISCARD = $(SWIFT_DIR)/discard
IO = $(SWIFT_DIR)/io
LYAP = $(SWIFT_DIR)/lyap
LYAP2 = $(SWIFT_DIR)/lyap2
MVITS = $(SWIFT_DIR)/mvits
MVS = $(SWIFT_DIR)/mvs
DRIFT = $(SWIFT_DIR)/mvs/drift
GETACCH = $(SWIFT_DIR)/mvs/getacch
KICKVH = $(SWIFT_DIR)/mvs/kickvh
STEP =  $(SWIFT_DIR)/mvs/step
ORBEL = $(SWIFT_DIR)/orbel
WIMPS = $(SWIFT_DIR)/wimps
RMVS = $(SWIFT_DIR)/rmvs
RMVS2 = $(SWIFT_DIR)/rmvs2
RMVS3 = $(SWIFT_DIR)/rmvs3
TU4 = $(SWIFT_DIR)/tu4
OBL = $(SWIFT_DIR)/obl
UTIL = $(SWIFT_DIR)/util
HJS = $(SWIFT_DIR)/hjs
HJS_TIDE = $(SWIFT_DIR)/hjs_tide
SYMBA5 = $(SWIFT_DIR)/symba5
SYMBA6 = $(SWIFT_DIR)/symba6
SYMBA5Q = $(SWIFT_DIR)/symba5q
SYMBATR = $(SWIFT_DIR)/symbatr
HELIO = $(SWIFT_DIR)/helio
HELIOQ = $(SWIFT_DIR)/helioq

clean:
	rm -f *.o *~ core

cleanall: clean
	rm -f *.a

%.a:
	-rm $(SWIFT_DIR_LIB)/$@
	$(FC) $(FFLAGS) $(ANAL)/*.f 
	$(FC) $(FFLAGS) $(BS)/*.f 
	$(FC) $(FFLAGS) $(COORD)/*.f   
	$(FC) $(FFLAGS) $(DISCARD)/*.f   
	$(FC) $(FFLAGS) $(IO)/*.f   
	$(FC) $(FFLAGS) $(LYAP)/*.f   
	$(FC) $(FFLAGS) $(LYAP2)/*.f   
	$(FC) $(FFLAGS) $(MVITS)/*.f   
	$(FC) $(FFLAGS) $(DRIFT)/*.f   
	$(FC) $(FFLAGS) $(GETACCH)/*.f   
	$(FC) $(FFLAGS) $(KICKVH)/*.f   
	$(FC) $(FFLAGS) $(STEP)/*.f   
	$(FC) $(FFLAGS) $(ORBEL)/*.f   
	$(FC) $(FFLAGS) $(HJS)/*.f
	$(FC) $(FFLAGS) $(HJS_TIDE)/*.f
#	$(FC) $(WIMPS)/*.f   
	$(FC) $(FFLAGS) $(RMVS)/*.f   
	$(FC) $(FFLAGS) $(RMVS2)/*.f   
	$(FC) $(FFLAGS) $(RMVS3)/*.f   
	$(FC) $(FFLAGS) $(TU4)/*.f   
	$(FC) $(FFLAGS) $(OBL)/*.f   
	$(FC) $(FFLAGS) $(UTIL)/*.f   
	$(FC) $(FFLAGS) $(SYMBA5)/*.f   
	$(FC) $(FFLAGS) $(SYMBA6)/*.f  
	$(FC) $(FFLAGS) $(SYMBA5Q)/*.f  
	$(FC) $(FFLAGS) $(SYMBATR)/*.f  
	$(FC) $(FFLAGS) $(HELIOQ)/*.f  
	$(FC) $(FFLAGS) $(HELIO)/*.f  

	ar -rv $(SWIFT_DIR_LIB)/$@  *.o
	rm *.o

