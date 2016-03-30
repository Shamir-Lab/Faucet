CFLAGS+= --std=c++0x -g -O4 -D_FILE_OFFSET_BITS=64# needed to handle files > 2 GB on 32 bits systems
SRC_UTILS =../utils/Cap.cpp ../utils/DoubleKmer.cpp ../utils/ReadKmer.cpp ../utils/Bloom.cpp ../utils/Kmer.cpp \
../utils/JChecker.cpp ../utils/JunctionMap.cpp ../utils/Junction.cpp ../utils/JuncPairs.cpp 
SRC_READSCAN= ReadScanner.cpp Contig.cpp ContigNode.cpp ContigGraph.cpp ContigIterator.cpp
SRC_MINK= Mink.cpp
SRC_TEST= PairFinderTest.cpp
TEST_MAIN=PairFinderTest.cpp

EXEC=mink
OBJ_BOTH= $(SRC_UTILS:.cpp=.o) $(SRC_READSCAN:.cpp=.o) 
OBJ_MINK= $(SRC_MINK:.cpp=.o) 
OBJ_TEST= $(SRC_TEST:.cpp=.o)

all:
	$(MAKE) $(EXEC)

mink: $(OBJ_BOTH) $(OBJ_MINK) Mink.cpp 
	g++ --std=c++0x $(SRC_READSCAN) $(SRC_UTILS) $(SRC_MINK) -o mink $(CFLAGS)

minkTest: $(OBJ_BOTH) $(OBJ_TEST) $(TEST_MAIN)
	g++ --std=c++0x $(SRC_READSCAN) $(SRC_UTILS) $(SRC_TEST) -o minkTest $(CFLAGS)

%.o: %.cpp %.h
	g++ --std=c++0x -o $@ -c $< $(CFLAGS)

install:
	cp mink /usr/local/bin
	cp minia /usr/local/bin
