# Do NOT edit NOR commit changes to this file for configuration purposes,
# create your own customized make.conf instead, see make.conf.sample.
include make.conf

KMER_LEN  ?= 30

GSH_INC   ?= /usr/include
GPT_INC   ?= /usr/include
GPT_LIB   ?= /usr/lib
BOOST_INC ?= /usr/include/boost
BOOST_LIB ?= /usr/lib
MCQ_INC   ?= $(HOME)/src/concurrentqueue
FOLLY_INC ?= $(HOME)/src/folly
SIMD_INC  ?= $(HOME)/src/simdcomp/include
SIMD_LIB  ?= $(HOME)/src/simdcomp

BINARY = bin/sm-standalone
SOURCES = src/sm_*.cpp

CFLAGS = $(FLAGS) -std=c++11 -DKMER_LEN=$(KMER_LEN)

all: $(BINARY)

$(BINARY): $(SOURCES)
	g++ $(CFLAGS) -Isrc -I$(GSH_INC) -I$(GPT_INC) -L$(GPT_LIB) \
		-I$(BOOST_INC) -L$(BOOST_LIB) -I$(MCQ_INC) -I$(FOLLY_INC) \
		-I$(SIMD_INC) -L$(SIMD_LIB) \
		-o $(BINARY) $(SOURCES) \
		-lprofiler -lz -lpthread -lsimdcomp

clean:
	rm -f $(BINARY)
