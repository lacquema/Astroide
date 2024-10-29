#! /bin/bash

# Makefile to build libraries and executables

# Parameters to input
COMPILF = /opt/homebrew/bin/gfortran
PYTHON3 = /Users/lacquema/Astroide.env/bin/python3

# Others parameters
PARALLEL = NO
ADD_FLAGS =
LIB_FLAGS = -O3 -c
GEN_FLAGS = -O3
ALG_FLAGS = -O3
DIR = .

# Induced directories
CODE_DIR = $(DIR)/Code

LIB_DIR = $(CODE_DIR)/lib
BIN_DIR = $(CODE_DIR)/bin
SWIFT_DIR = $(CODE_DIR)/swift

SUB_DIR = $(SWIFT_DIR)/sub
ALG_DIR = $(SWIFT_DIR)/main

GEN_DIR = $(SUB_DIR)/gen
COLL_DIR = $(SUB_DIR)/coll
ANAL_DIR = $(SUB_DIR)/anal
BS_DIR = $(SUB_DIR)/bs
COORD_DIR = $(SUB_DIR)/coord
DISCARD_DIR = $(SUB_DIR)/discard
IO_DIR = $(SUB_DIR)/io
LYAP_DIR = $(SUB_DIR)/lyap
LYAP2_DIR = $(SUB_DIR)/lyap2
MVITS_DIR = $(SUB_DIR)/mvits
MVS_DIR = $(SUB_DIR)/mvs
DRIFT_DIR = $(SUB_DIR)/mvs/drift
GETACCH_DIR = $(SUB_DIR)/mvs/getacch
KICKVH_DIR = $(SUB_DIR)/mvs/kickvh
STEP_DIR =  $(SUB_DIR)/mvs/step
ORBEL_DIR = $(SUB_DIR)/orbel
WIMPS_DIR = $(SUB_DIR)/wimps
RMVS_DIR = $(SUB_DIR)/rmvs
RMVS2_DIR = $(SUB_DIR)/rmvs2
RMVS3_DIR = $(SUB_DIR)/rmvs3
TU4_DIR = $(SUB_DIR)/tu4
OBL_DIR = $(SUB_DIR)/obl
UTIL_DIR = $(SUB_DIR)/util
HJS_DIR = $(SUB_DIR)/hjs
HJS_TIDE_DIR = $(SUB_DIR)/hjs_tide
SYMBA5_DIR = $(SUB_DIR)/symba5
SYMBA6_DIR = $(SUB_DIR)/symba6
SYMBA5Q_DIR = $(SUB_DIR)/symba5q
SYMBATR_DIR = $(SUB_DIR)/symbatr
HELIO_DIR = $(SUB_DIR)/helio
HELIOQ_DIR = $(SUB_DIR)/helioq

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
	test ! -f $(LIB) || rm $(LIB)
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(ANAL_DIR)/*.f 
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(BS_DIR)/*.f 
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(COORD_DIR)/*.f   
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(DISCARD_DIR)/*.f   
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(IO_DIR)/*.f   
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(LYAP_DIR)/*.f   
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(LYAP2_DIR)/*.f   
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(MVITS_DIR)/*.f   
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(DRIFT_DIR)/*.f   
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(GETACCH_DIR)/*.f   
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(KICKVH_DIR)/*.f   
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(STEP_DIR)/*.f   
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(ORBEL_DIR)/*.f   
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(HJS_DIR)/*.f
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(HJS_TIDE_DIR)/*.f
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(RMVS_DIR)/*.f   
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(RMVS2_DIR)/*.f   
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(RMVS3_DIR)/*.f   
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(TU4_DIR)/*.f   
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(OBL_DIR)/*.f   
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(UTIL_DIR)/*.f   
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(SYMBA5_DIR)/*.f   
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(SYMBA6_DIR)/*.f  
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(SYMBA5Q_DIR)/*.f  
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(SYMBATR_DIR)/*.f  
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(HELIOQ_DIR)/*.f  
	$(COMPILF) $(LIB_FLAGS) $(ADD_FLAGS) $(HELIO_DIR)/*.f  
	ar -rv $(LIB) *.o
	rm *.o

packages:
	$(PYTHON3) -m pip install -r $(DIR)/requirements.txt

gen_tout_multi:
	$(COMPILF) $(GEN_FLAGS) $(ADD_FLAGS) $(GEN_DIR)/$@.f $(LIB) -o $(GEN)

# Algorithms
swift_%:
ifeq ($(PARALLEL),NO)
	test ! -f $(BIN_DIR)/$@ || rm $(BIN_DIR)/$@
	$(COMPILF) $(ALG_FLAGS) $(ADD_FLAGS) $(ALG_DIR)/$@.f $(LIB) -o $(BIN_DIR)/$@
else
	test ! -f $(BIN_DIR)/$@_par || rm $(BIN_DIR)/$@_par
	$(COMPILF) $(ALG_FLAGS) $(ADD_FLAGS) $(ALG_DIR)/$@.f $(LIB) -o $(BIN_DIR)/$@_par
endif

swift_whm:

swift_rmvs3:

swift_hjs:

all: 
	make packages
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


