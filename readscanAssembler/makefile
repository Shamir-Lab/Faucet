CFLAGS+= -O4 -D_FILE_OFFSET_BITS=64 # needed to handle files > 2 GB on 32 bits systems
SRC_UTILS =../utils/Cap.cpp ../utils/DoubleKmer.cpp ../utils/ReadKmer.cpp ../utils/Bloom.cpp ../utils/Kmer.cpp ../utils/JChecker.cpp ../utils/JunctionMap.cpp ../utils/Junction.cpp
SRC_MINK =Graph.cpp Node.cpp ReadScanner.cpp Mink.cpp

EXEC=mink
OBJ= $(SRC_UTILS:.cpp=.o) $(SRC_MINK:.cpp=.o)

all:
	$(MAKE) $(EXEC)

mink: $(OBJ) Mink.cpp 
	g++ --std=c++0x $(SRC_MINK) $(SRC_UTILS) -o mink $(CFLAGS)

%.o: %.cpp %.h
	g++ --std=c++0x -o $@ -c $< $(CFLAGS)

install:
	cp mink /usr/local/bin
	cp minia /usr/local/bin
