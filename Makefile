# Do NOT edit NOR commit changes to this file for configuration purposes,
# create your own customized make.conf instead, see make.conf.sample.
include make.conf

VERBOSE ?= 0
CC_0 = @echo "CC $@"; g++
CC_1 = g++
CC = $(CC_$(VERBOSE))

GSH_INC   ?= /usr/include/sparsehash
MCQ_INC   ?= /usr/include/concurrentqueue
FOLLY_INC ?= /usr/include/folly
BOOST_INC ?= /usr/include/boost
ROCKS_INC ?= /usr/include/rocksdb
ROCKS_LIB ?= /usr/lib

BIN = sm
SRC = $(wildcard src/*.cpp)
OBJ = $(SRC:.cpp=.o)
DEP = $(SRC:.cpp=.d)

INC = -Isrc -I$(GSH_INC) -I$(MCQ_INC) -I$(FOLLY_INC) \
      -I$(BOOST_INC) -I$(ROCKS_INC)
LIB = -lz -lpthread -lrocksdb

CFLAGS += -std=c++11
LFLAGS += -L$(ROCKS_LIB)

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
