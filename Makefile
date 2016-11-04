# Do NOT edit NOR commit changes to this file for configuration purposes,
# create your own customized make.conf instead, see make.conf.sample.
include make.conf

KMER_LEN ?= 30

MIN_TC ?= 4
MAX_NC ?= 1

WMIN ?= 7
WLEN ?= 10

GSH_INC   ?= /usr/include
BOOST_INC ?= /usr/include/boost
BOOST_LIB ?= /usr/lib
MCQ_INC   ?= $(HOME)/src/concurrentqueue
FOLLY_INC ?= $(HOME)/src/folly
RDB_INC   ?= /usr/include/rocksdb
RDB_LIB   ?= /usr/lib

MAIN_BIN = bin/sm

MAIN_SRC = src/db.cpp src/common.cpp src/count.cpp src/filter.cpp \
           src/merge.cpp src/group.cpp src/hash.cpp src/stage.cpp src/main.cpp

CFLAGS = $(FLAGS) -std=c++11 -DKMER_LEN=$(KMER_LEN) \
         -DMIN_TC=$(MIN_TC) -DMAX_NC=$(MAX_NC) \
         -DWMIN=$(WMIN) -DWLEN=$(WLEN)

all: $(MAIN_BIN)

$(MAIN_BIN): $(MAIN_SRC)
	g++ $(CFLAGS) -Isrc \
		-I$(GSH_INC) -I$(BOOST_INC) -L$(BOOST_LIB) -I$(MCQ_INC) \
		-I$(FOLLY_INC) -I$(RDB_INC) \
		-o $(MAIN_BIN) $(MAIN_SRC) \
		$(RDB_LIB) -lz -lpthread

clean:
	rm -f $(MAIN_BIN)
