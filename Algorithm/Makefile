#! /bin/bash

# Makefile to build libraries and executables

# Parameters
F = gfortran
PARALLEL = NO
ADD_FLAGS =
LIB_FLAGS = -O3 -c
GEN_FLAGS = -O3
ALG_FLAGS = -O3
DIR = .

# Induced directories
SWIFT_DIR = $(DIR)/swift
LIB_DIR = $(DIR)/lib
BIN_DIR = $(DIR)/bin
GEN_DIR = $(SWIFT_DIR)/gen
ALG_DIR = $(SWIFT_DIR)/main

COLL_DIR = $(SWIFT_DIR)/coll
ANAL_DIR = $(SWIFT_DIR)/anal
BS_DIR = $(SWIFT_DIR)/bs
COORD_DIR = $(SWIFT_DIR)/coord
DISCARD_DIR = $(SWIFT_DIR)/discard
IO_DIR = $(SWIFT_DIR)/io
LYAP_DIR = $(SWIFT_DIR)/lyap
LYAP2_DIR = $(SWIFT_DIR)/lyap2
MVITS_DIR = $(SWIFT_DIR)/mvits
MVS_DIR = $(SWIFT_DIR)/mvs
DRIFT_DIR = $(SWIFT_DIR)/mvs/drift
GETACCH_DIR = $(SWIFT_DIR)/mvs/getacch
KICKVH_DIR = $(SWIFT_DIR)/mvs/kickvh
STEP_DIR =  $(SWIFT_DIR)/mvs/step
ORBEL_DIR = $(SWIFT_DIR)/orbel
WIMPS_DIR = $(SWIFT_DIR)/wimps
RMVS_DIR = $(SWIFT_DIR)/rmvs
RMVS2_DIR = $(SWIFT_DIR)/rmvs2
RMVS3_DIR = $(SWIFT_DIR)/rmvs3
TU4_DIR = $(SWIFT_DIR)/tu4
OBL_DIR = $(SWIFT_DIR)/obl
UTIL_DIR = $(SWIFT_DIR)/util
HJS_DIR = $(SWIFT_DIR)/hjs
HJS_TIDE_DIR = $(SWIFT_DIR)/hjs_tide
SYMBA5_DIR = $(SWIFT_DIR)/symba5
SYMBA6_DIR = $(SWIFT_DIR)/symba6
SYMBA5Q_DIR = $(SWIFT_DIR)/symba5q
SYMBATR_DIR = $(SWIFT_DIR)/symbatr
HELIO_DIR = $(SWIFT_DIR)/helio
HELIOQ_DIR = $(SWIFT_DIR)/helioq

# Parallelization option 
ifeq ($(PARALLEL),NO)
LIB = $(LIB_DIR)/libswift.a
GEN = $(BIN_DIR)/gen_tout_multi
else
GEN_FLAGS += -fopenmp
LIB_FLAGS += -fopenmp
ALG_FLAGS += -fopenmp
LIB = $(LIB_DIR)/libswift_par.a
GEN = $(BIN_DIR)/gen_tout_multi_par
endif

# Utilities
library:
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(ANAL_DIR)/*.f 
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(BS_DIR)/*.f 
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(COORD_DIR)/*.f   
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(DISCARD_DIR)/*.f   
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(IO_DIR)/*.f   
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(LYAP_DIR)/*.f   
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(LYAP2_DIR)/*.f   
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(MVITS_DIR)/*.f   
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(DRIFT_DIR)/*.f   
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(GETACCH_DIR)/*.f   
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(KICKVH_DIR)/*.f   
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(STEP_DIR)/*.f   
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(ORBEL_DIR)/*.f   
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(HJS_DIR)/*.f
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(HJS_TIDE_DIR)/*.f
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(RMVS_DIR)/*.f   
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(RMVS2_DIR)/*.f   
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(RMVS3_DIR)/*.f   
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(TU4_DIR)/*.f   
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(OBL_DIR)/*.f   
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(UTIL_DIR)/*.f   
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(SYMBA5_DIR)/*.f   
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(SYMBA6_DIR)/*.f  
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(SYMBA5Q_DIR)/*.f  
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(SYMBATR_DIR)/*.f  
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(HELIOQ_DIR)/*.f  
	$(F) $(LIB_FLAGS) $(ADD_FLAGS) $(HELIO_DIR)/*.f  
	ar -rv $(LIB) *.o
	rm *.o

gen_tout_multi:
	$(F) $(GEN_FLAGS) $(ADD_FLAGS) $(GEN_DIR)/$@.f $(LIB) -o $(GEN)

# Algorithms
swift_%:
ifeq ($(PARALLEL),NO)
	$(F) $(ALG_FLAGS) $(ADD_FLAGS) $(ALG_DIR)/$@.f $(LIB) -o $(BIN_DIR)/$@
else
	$(F) $(ALG_FLAGS) $(ADD_FLAGS) $(ALG_DIR)/$@.f $(LIB) -o $(BIN_DIR)/$@_par
endif

swift_whm:

swift_rmvs3:

swift_hjs:

all: 
	make library
	make gen_tout_multi
	make swift_whm
	make swift_rmvs3
	make swift_hjs

binaries: 
	make gen_tout_multi
	make swift_whm
	make swift_rmvs3
	make swift_hjs

clean:
	rm *.o

cleanall: 
	clean
	rm *.a


