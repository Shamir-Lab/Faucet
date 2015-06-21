CFLAGS+= -O4 -D_FILE_OFFSET_BITS=64 # needed to handle files > 2 GB on 32 bits systems
SRC_UTILS      =../utils/Bloom.cpp ../utils/Kmer.cpp ../utils/JChecker.cpp ../utils/JunctionMap.cpp ../utils/Junction.cpp
SRC_MINIA=KpomerScanner.cpp Minia.cpp

EXEC=minia
OBJ= $(SRC_UTILS:.cpp=.o) $(SRC_MINIA:.cpp=.o) 

all:
	$(MAKE) $(EXEC)

minia: $(OBJ) Minia.cpp 
	g++ --std=c++0x $(SRC_MINIA) $(SRC_UTILS) -o minia $(CFLAGS)

%.o: %.cpp %.h
	g++ --std=c++0x -o $@ -c $< $(CFLAGS)

install:
	cp mink /usr/local/bin
	cp minia /usr/local/bin
