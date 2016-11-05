# Do NOT edit NOR commit changes to this file for configuration purposes,
# create your own customized make.conf instead, see make.conf.sample.
include make.conf

KMER_LEN ?= 30

MIN_TC ?= 4
MAX_NC ?= 1

WMIN ?= 7
WLEN ?= 10

VERBOSE ?= 0
CC_0 = @echo "CC $@"; g++
CC_1 = g++
CC = $(CC_$(VERBOSE))

GSH_INC   ?= /usr/include/sparsehash
MCQ_INC   ?= /usr/include/concurrentqueue
FOLLY_INC ?= /usr/include/folly
BOOST_INC ?= /usr/include/boost
ROCKS_INC ?= /usr/include/rocksdb
BOOST_LIB ?= /usr/lib
ROCKS_LIB ?= /usr/lib

BIN = sm
SRC = $(wildcard src/*.cpp)
OBJ = $(SRC:.cpp=.o)
DEP = $(SRC:.cpp=.d)

DEF = -DKMER_LEN=$(KMER_LEN) \
      -DMIN_TC=$(MIN_TC) -DMAX_NC=$(MAX_NC) \
      -DWMIN=$(WMIN) -DWLEN=$(WLEN)
INC = -Isrc -I$(GSH_INC) -I$(MCQ_INC) -I$(FOLLY_INC) \
      -I$(BOOST_INC) -I$(ROCKS_INC)
LIB = -lz -lpthread -lrocksdb

CFLAGS += -std=c++11
LFLAGS += -L$(ROCKS_LIB)

all: $(BIN)

-include $(DEP)

.cpp.o:
	$(CC) $(CFLAGS) $(DEF) $(INC) -MMD -c -o $@ $<

$(BIN): $(OBJ)
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN) $(OBJ) $(LIB)

clean:
	rm -f $(BIN)

distclean: clean
	rm -f $(OBJ)
	rm -f $(DEP)
