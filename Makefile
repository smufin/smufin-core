BINARIES = bin/sm-standalone

SOURCES = src/sm_*.cpp

GSH_INCLUDE = $(HOME)/opt/sparsehash/include
GPT_INCLUDE = $(HOME)/opt/gperftools/include
GPT_LIB     = $(HOME)/opt/gperftools/lib

BOOST_INCLUDE = $(HOME)/src/boost/boost_1_57_0
BOOST_LIB     = $(HOME)/src/boost/boost_1_57_0/stage/lib

CQ_INCLUDE  = $(HOME)/src/concurrentqueue
RWQ_INCLUDE = $(HOME)/src/readerwriterqueue
FF_INCLUDE  = $(HOME)/src/folly

# CFLAGS = -O2 -std=c++11 -mcpu=power8 -DNDEBUG -DPROFILE
CFLAGS = -O3 -std=c++11 -mcpu=power8 -DNDEBUG

all: $(BINARIES)

bin/sm-standalone: $(SOURCES)
	g++ $(CFLAGS) -Isrc -I$(GSH_INCLUDE) -I$(GPT_INCLUDE) -L$(GPT_LIB) \
		-I$(BOOST_INCLUDE) -L$(BOOST_LIB) -I$(CQ_INCLUDE) -I$(FF_INCLUDE)\
		-o bin/sm-standalone $(SOURCES) \
		-lprofiler -lz -lpthread

INPUT = var/input-chr22
MAP = var/map-4char-2proc-16srvs

launch: bin/sm-standalone
	@seq 0 1 | xargs -P 2 -n 1 -I ID bash -c \
		'/usr/bin/time -o ID.time ./bin/sm-standalone -i $(INPUT) -m $(MAP) -p ID > ID.out'

clean:
	rm -f $(BINARIES)
