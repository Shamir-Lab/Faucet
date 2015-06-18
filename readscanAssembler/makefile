CFLAGS+= -O4 -D_FILE_OFFSET_BITS=64 # needed to handle files > 2 GB on 32 bits systems
SRC      =Bloom.cpp Kmer.cpp JChecker.cpp JunctionMap.cpp Junction.cpp
SRC_MINIA=KpomerScanner.cpp Minia.cpp
SRC_MINK =ReadScanner.cpp Mink.cpp

EXEC=mink minia
OBJ= $(SRC:.cpp=.o) $(SRC_MINIA:.cpp=.o) $(SRC_MINK:.cpp=.o)

all:
	$(MAKE) $(EXEC)

mink: $(OBJ) Mink.cpp 
	g++ --std=c++0x $(SRC_MINK) $(SRC) -o mink $(CFLAGS)

minia: $(OBJ) Minia.cpp 
	g++ --std=c++0x $(SRC_MINIA) $(SRC) -o minia $(CFLAGS)

%.o: %.cpp %.h
	g++ --std=c++0x -o $@ -c $< $(CFLAGS)

install:
	cp mink /usr/local/bin
	cp minia /usr/local/bin
