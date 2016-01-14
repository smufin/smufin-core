# Do NOT edit NOR commit changes to this file for configuration purposes,
# create your own customized make.conf instead, see make.conf.sample.
include make.conf

READ_LEN  ?= 100
KMER_LEN  ?= 30

GSH_INC   ?= /usr/include
GPT_INC   ?= /usr/include
GPT_LIB   ?= /usr/lib
BOOST_INC ?= /usr/include/boost
BOOST_LIB ?= /usr/lib
MCQ_INC   ?= $(HOME)/src/concurrentqueue
FOLLY_INC ?= $(HOME)/src/folly

FILTER_BIN = bin/sm-filter
FILTER_SRC = src/common.cpp src/process.cpp src/filter.cpp \
             src/hash.cpp src/main_filter.cpp

GROUP_BIN = bin/sm-group
GROUP_SRC = src/main_group.cpp

CFLAGS = $(FLAGS) -std=c++11 -DREAD_LEN=$(READ_LEN) -DKMER_LEN=$(KMER_LEN)

all: $(FILTER_BIN) $(GROUP_BIN)

$(FILTER_BIN): $(FILTER_SRC)
	g++ $(CFLAGS) -Isrc -I$(GSH_INC) -I$(GPT_INC) -L$(GPT_LIB) \
		-I$(BOOST_INC) -L$(BOOST_LIB) -I$(MCQ_INC) -I$(FOLLY_INC) \
		-o $(FILTER_BIN) $(FILTER_SRC) \
		-lprofiler -lz -lpthread

$(GROUP_BIN): $(GROUP_SRC)
	g++ $(CFLAGS) -Isrc -I$(GSH_INC) -I$(BOOST_INC) -L$(BOOST_LIB) \
		-o $(GROUP_BIN) $(GROUP_SRC) \
		-lz -lboost_iostreams

clean:
	rm -f $(FILTER_BIN) $(GROUP_BIN)
