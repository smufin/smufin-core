# Do NOT edit NOR commit changes to this file for configuration purposes,
# create your own customized make.conf instead, see make.conf.sample.
include make.conf

KMER_LEN  ?= 30

MIN_TC ?= 4
MAX_NC ?= 1

WMIN ?= 7
WLEN ?= 10

GSH_INC   ?= /usr/include
GPT_INC   ?= /usr/include
GPT_LIB   ?= /usr/lib
BOOST_INC ?= /usr/include/boost
BOOST_LIB ?= /usr/lib
MCQ_INC   ?= $(HOME)/src/concurrentqueue
FOLLY_INC ?= $(HOME)/src/folly
RDB_INC   ?= /usr/include/rocksdb
RDB_LIB   ?= /usr/lib

PROCESS_BIN = bin/sm-process

FILTER_BIN = bin/sm-filter
FILTER_SRC = src/common.cpp src/process.cpp src/filter.cpp \
             src/hash.cpp src/main_filter.cpp

GROUP_BIN = bin/sm-group
GROUP_SRC = src/common.cpp src/main_group.cpp

MERGE_BIN = bin/sm-merge
MERGE_SRC = src/db.cpp src/merge.cpp src/main_merge.cpp

CFLAGS = $(FLAGS) -std=c++11 -DKMER_LEN=$(KMER_LEN) \
         -DMIN_TC=$(MIN_TC) -DMAX_NC=$(MAX_NC) \
         -DWMIN=$(WMIN) -DWLEN=$(WLEN)

all: $(PROCESS_BIN) $(FILTER_BIN) $(GROUP_BIN) $(MERGE_BIN)

$(PROCESS_BIN): $(FILTER_SRC)
	g++ $(CFLAGS) -DENABLE_PROCESS -Isrc -I$(GSH_INC) -I$(GPT_INC) \
		-L$(GPT_LIB) -I$(BOOST_INC) -L$(BOOST_LIB) -I$(MCQ_INC) \
		-I$(FOLLY_INC) -o $(PROCESS_BIN) $(FILTER_SRC) -lz -lpthread

$(FILTER_BIN): $(FILTER_SRC)
	g++ $(CFLAGS) -DENABLE_FILTER -Isrc -I$(GSH_INC) -I$(GPT_INC) \
		-L$(GPT_LIB) -I$(BOOST_INC) -L$(BOOST_LIB) -I$(MCQ_INC) \
		-I$(FOLLY_INC) -o $(FILTER_BIN) $(FILTER_SRC) -lz -lpthread

$(GROUP_BIN): $(GROUP_SRC)
	g++ $(CFLAGS) -Isrc -I$(GSH_INC) -I$(BOOST_INC) -L$(BOOST_LIB) \
		-I$(BOOST_INC) -L$(BOOST_LIB) -I$(MCQ_INC) -I$(FOLLY_INC) \
		-o $(GROUP_BIN) $(GROUP_SRC) \
		-lz -lboost_iostreams

$(MERGE_BIN): $(MERGE_SRC)
	g++ $(CFLAGS) -Isrc -I$(MCQ_INC) -I$(FOLLY_INC) -I$(RDB_INC) \
		-o $(MERGE_BIN) $(MERGE_SRC) $(RDB_LIB)

clean:
	rm -f $(PROCESS_BIN) $(FILTER_BIN) $(GROUP_BIN) $(MERGE_BIN)
