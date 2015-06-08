CFLAGS+= -O4 -D_FILE_OFFSET_BITS=64 # needed to handle files > 2 GB on 32 bits systems
SRC=Bloom.cpp Kmer.cpp Debloom.cpp Minia.cpp
EXEC=minia
OBJ= $(SRC:.cpp=.o)

all:
	$(MAKE) $(EXEC)

minia: $(OBJ) Minia.cpp
	g++ --std=c++0x $(SRC) -o minia 

%.o: %.cpp %.h
	g++ --std=c++0x -o $@ -c $< $(CFLAGS)

install:
	cp minia /usr/local/bin
