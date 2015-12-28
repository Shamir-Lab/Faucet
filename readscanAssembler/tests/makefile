CFLAGS+= -O4 -D_FILE_OFFSET_BITS=64 # needed to handle files > 2 GB on 32 bits systems
SRC= PairFinderTest.cpp ../Contig.cpp

EXEC=minkTests
OBJ= $(SRC:.cpp=.o)

all:
	$(MAKE) $(EXEC)

minkTests: $(OBJ) PairFinderTest.cpp
	cd .. && $(MAKE)
	g++ --std=c++0x $(SRC) -o minkTests $(CFLAGS) 

%.o: %.cpp %.h
	g++ --std=c++0x -o $@ -c $< $(CFLAGS)

install:
	cp minkTests /usr/local/bin
