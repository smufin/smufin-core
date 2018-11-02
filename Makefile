# Do NOT edit NOR commit changes to this file for configuration purposes,
# create your own customized make.conf instead, see make.conf.sample.
include make.conf

VERBOSE ?= 0
CC_0 = @echo "CC $@"; g++
CC_1 = g++
CC = $(CC_$(VERBOSE))

MAX_READ_LEN ?= 120

GSH_INC   ?= /usr/include/sparsehash
MCQ_INC   ?= /usr/include/concurrentqueue
RWQ_INC   ?= /usr/include/readerwriterqueue
BOOST_INC ?= /usr/include/boost
BF_INC    ?= /usr/include/libbf
BF_LIB    ?= /usr/lib
ROCKS_INC ?= /usr/include/rocksdb
ROCKS_LIB ?= /usr/lib
HTS_INC   ?= /usr/include/htslib
HTS_LIB   ?= /usr/lib

BIN = sm
SRC = $(wildcard src/*.cpp)
OBJ = $(SRC:.cpp=.o)
DEP = $(SRC:.cpp=.d)

INC = -Isrc -I$(GSH_INC) -I$(MCQ_INC) -I$(RWQ_INC) \
      -I$(BOOST_INC) -I$(BF_INC) -I$(ROCKS_INC) -I$(HTS_INC)
LIB = -lboost_iostreams -lz -lpthread -lbf -lrocksdb -lhts

CFLAGS += -std=c++11 -DMAX_READ_LEN=$(MAX_READ_LEN)
LFLAGS += -L$(BF_LIB) -L$(ROCKS_LIB) -L$(HTS_LIB)

all: $(BIN)

-include $(DEP)

.cpp.o:
	$(CC) $(CFLAGS) $(INC) -MMD -c -o $@ $<

$(BIN): $(OBJ)
	$(CC) $(CFLAGS) $(LFLAGS) -o $(BIN) $(OBJ) $(LIB)

clean:
	rm -f $(BIN)

distclean: clean
	rm -f $(OBJ)
	rm -f $(DEP)
