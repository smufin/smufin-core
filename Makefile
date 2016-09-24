# Do NOT edit NOR commit changes to this file for configuration purposes,
# create your own customized make.conf instead, see make.conf.sample.
include make.conf

KMER_LEN  ?= 30

MIN_TC ?= 4
MAX_NC ?= 1

WMIN ?= 7
WLEN ?= 10

GSH_INC   ?= /usr/include
BOOST_INC ?= /usr/include/boost
BOOST_LIB ?= /usr/lib
MCQ_INC   ?= $(HOME)/src/concurrentqueue
FOLLY_INC ?= $(HOME)/src/folly

MAIN_BIN = bin/sm
MAIN_SRC = src/common.cpp src/count.cpp src/filter.cpp src/hash.cpp \
           src/stage.cpp src/main.cpp

GROUP_BIN = bin/sm-group
GROUP_SRC = src/common.cpp src/main_group.cpp

JOINF_BIN = bin/sm-join-fq
JOINF_SRC = src/main_join_fq.cpp

CFLAGS = $(FLAGS) -std=c++11 -DKMER_LEN=$(KMER_LEN) \
         -DMIN_TC=$(MIN_TC) -DMAX_NC=$(MAX_NC) \
         -DWMIN=$(WMIN) -DWLEN=$(WLEN)

all: $(MAIN_BIN)

$(MAIN_BIN): $(MAIN_SRC)
	g++ $(CFLAGS) -Isrc \
		-I$(GSH_INC) -I$(BOOST_INC) -L$(BOOST_LIB) -I$(MCQ_INC) -I$(FOLLY_INC) \
		-o $(MAIN_BIN) $(MAIN_SRC) \
		-lz -lpthread

$(GROUP_BIN): $(GROUP_SRC)
	g++ $(CFLAGS) -Isrc -I$(GSH_INC) -I$(BOOST_INC) -L$(BOOST_LIB) \
		-I$(BOOST_INC) -L$(BOOST_LIB) -I$(MCQ_INC) -I$(FOLLY_INC) \
		-o $(GROUP_BIN) $(GROUP_SRC) \
		-lz -lboost_iostreams

$(JOINF_BIN): $(JOINF_SRC)
	g++ $(CFLAGS) -Isrc -I$(GSH_INC) -I$(BOOST_INC) -L$(BOOST_LIB) \
		-o $(JOINF_BIN) $(JOINF_SRC)

clean:
	rm -f $(MAIN_BIN) $(GROUP_BIN) $(JOINF_BIN)
