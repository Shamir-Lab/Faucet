CFLAGS+= --std=c++0x -g -O4 -D_FILE_OFFSET_BITS=64# needed to handle files > 2 GB on 32 bits systems

# Prefix for readscan files and util files
READSCAN_PREFIX =./
UTIL_PREFIX =../utils/

# List of just filenames for utils and readscanAssembler
UTIL_FILES =Cap.cpp DoubleKmer.cpp ReadKmer.cpp Bloom.cpp Kmer.cpp JChecker.cpp JunctionMap.cpp Junction.cpp JuncPairs.cpp
READSCAN_FILES= ReadScanner.cpp Contig.cpp ContigNode.cpp ContigGraph.cpp ContigIterator.cpp

# Full path to files
SRC_READSCAN= $(addprefix $(READSCAN_PREFIX),$(READSCAN_FILES))
SRC_UTILS= $(addprefix $(UTIL_PREFIX),$(UTIL_FILES))

SRC_MINK= Mink.cpp

EXEC=mink

#Object files for readscan utils and mink
OBJ_BOTH= $(SRC_UTILS:.cpp=.o) $(SRC_READSCAN:.cpp=.o) 
OBJ_MINK= $(SRC_MINK:.cpp=.o) 
$(info OBJ_BOTH is $(OBJ_BOTH))

mink: $(OBJ_BOTH) $(OBJ_MINK)
	g++ --std=c++0x $(SRC_READSCAN) $(SRC_UTILS) $(SRC_MINK) -o mink $(CFLAGS)

%.o: %.cpp %.h
	g++ -o $@ -c $< $(CFLAGS)

install:
	cp mink /usr/local/bin
	cp minia /usr/local/bin

clean:
	\rm *.o ../utils/*.o mink

# MAKE DEPEND OUTPUT
# To regenerate this when dependencies change, run "bash make_depend.sh" from this directory.
# DO NOT DELETE

../utils/Cap.o: ../utils/Cap.h ../utils/Kmer.h /usr/include/stdint.h
../utils/Cap.o: /usr/include/features.h /usr/include/stdc-predef.h
../utils/DoubleKmer.o: ../utils/DoubleKmer.h ../utils/Kmer.h
../utils/DoubleKmer.o: /usr/include/stdint.h /usr/include/features.h
../utils/DoubleKmer.o: /usr/include/stdc-predef.h
../utils/ReadKmer.o: ../utils/ReadKmer.h ../utils/Kmer.h
../utils/ReadKmer.o: /usr/include/stdint.h /usr/include/features.h
../utils/ReadKmer.o: /usr/include/stdc-predef.h ../utils/DoubleKmer.h
../utils/Bloom.o: /usr/include/stdio.h /usr/include/features.h
../utils/Bloom.o: /usr/include/stdc-predef.h /usr/include/libio.h
../utils/Bloom.o: /usr/include/_G_config.h /usr/include/wchar.h
../utils/Bloom.o: /usr/include/stdlib.h /usr/include/alloca.h
../utils/Bloom.o: /usr/include/string.h /usr/include/xlocale.h
../utils/Bloom.o: /usr/include/inttypes.h /usr/include/stdint.h
../utils/Bloom.o: /usr/include/math.h ../utils/Bloom.h ../utils/Kmer.h
../utils/Bloom.o: ../utils/ReadKmer.h ../utils/DoubleKmer.h
../utils/Bloom.o: ../utils/JuncPairs.h
../utils/Kmer.o: /usr/include/stdio.h /usr/include/features.h
../utils/Kmer.o: /usr/include/stdc-predef.h /usr/include/libio.h
../utils/Kmer.o: /usr/include/_G_config.h /usr/include/wchar.h
../utils/Kmer.o: /usr/include/assert.h ../utils/Kmer.h /usr/include/stdint.h
../utils/Kmer.o: ../utils/lut.h
../utils/JChecker.o: ../utils/JChecker.h ../utils/Bloom.h
../utils/JChecker.o: /usr/include/stdlib.h /usr/include/features.h
../utils/JChecker.o: /usr/include/stdc-predef.h /usr/include/alloca.h
../utils/JChecker.o: /usr/include/inttypes.h /usr/include/stdint.h
../utils/JChecker.o: ../utils/Kmer.h ../utils/ReadKmer.h
../utils/JChecker.o: ../utils/DoubleKmer.h ../utils/JuncPairs.h
../utils/JChecker.o: /usr/include/stdio.h /usr/include/libio.h
../utils/JChecker.o: /usr/include/_G_config.h /usr/include/wchar.h
../utils/JunctionMap.o: ../utils/JunctionMap.h ../utils/Kmer.h
../utils/JunctionMap.o: /usr/include/stdint.h /usr/include/features.h
../utils/JunctionMap.o: /usr/include/stdc-predef.h ../utils/Junction.h
../utils/JunctionMap.o: ../utils/Cap.h ../utils/ReadKmer.h
../utils/JunctionMap.o: ../utils/DoubleKmer.h ../utils/Bloom.h
../utils/JunctionMap.o: /usr/include/stdlib.h /usr/include/alloca.h
../utils/JunctionMap.o: /usr/include/inttypes.h ../utils/JuncPairs.h
../utils/JunctionMap.o: ../utils/JChecker.h ../readscanAssembler/Contig.h
../utils/JunctionMap.o: ../utils/Kmer.h ../utils/JuncPairs.h ../utils/Bloom.h
../utils/JunctionMap.o: ContigNode.h Contig.h ../utils/Junction.h
../utils/JunctionMap.o: /usr/include/string.h /usr/include/xlocale.h
../utils/JunctionMap.o: ../readscanAssembler/ContigNode.h
../utils/JunctionMap.o: ../readscanAssembler/ContigGraph.h
../utils/JunctionMap.o: ../utils/JunctionMap.h ../utils/JChecker.h
../utils/JunctionMap.o: BfSearchResult.h ContigIterator.h ContigGraph.h
../utils/JunctionMap.o: ../readscanAssembler/BfSearchResult.h
../utils/Junction.o: ../utils/Junction.h ../utils/Kmer.h
../utils/Junction.o: /usr/include/stdint.h /usr/include/features.h
../utils/Junction.o: /usr/include/stdc-predef.h /usr/include/limits.h
../utils/JuncPairs.o: ../utils/JuncPairs.h ../utils/Kmer.h
../utils/JuncPairs.o: /usr/include/stdint.h /usr/include/features.h
../utils/JuncPairs.o: /usr/include/stdc-predef.h
ReadScanner.o: ReadScanner.h /usr/include/stdio.h /usr/include/features.h
ReadScanner.o: /usr/include/stdc-predef.h /usr/include/libio.h
ReadScanner.o: /usr/include/_G_config.h /usr/include/wchar.h
ReadScanner.o: /usr/include/stdlib.h /usr/include/alloca.h
ReadScanner.o: /usr/include/time.h /usr/include/xlocale.h
ReadScanner.o: /usr/include/string.h /usr/include/inttypes.h
ReadScanner.o: /usr/include/stdint.h /usr/include/unistd.h
ReadScanner.o: /usr/include/getopt.h ../utils/Bloom.h ../utils/Kmer.h
ReadScanner.o: ../utils/JChecker.h ../utils/JunctionMap.h ../utils/Junction.h
ReadScanner.o: ../utils/ReadKmer.h ../utils/Kmer.h ../utils/DoubleKmer.h
ReadScanner.o: ../utils/Cap.h ../utils/JuncPairs.h
Contig.o: Contig.h
ContigNode.o: ContigNode.h Contig.h ../utils/Kmer.h ../utils/JuncPairs.h
ContigNode.o: ../utils/Junction.h /usr/include/time.h /usr/include/features.h
ContigNode.o: /usr/include/stdc-predef.h /usr/include/xlocale.h
ContigGraph.o: ContigGraph.h /usr/include/math.h /usr/include/features.h
ContigGraph.o: /usr/include/stdc-predef.h /usr/include/time.h
ContigGraph.o: /usr/include/xlocale.h /usr/include/stdlib.h
ContigGraph.o: /usr/include/alloca.h
ContigIterator.o: ContigIterator.h ContigGraph.h Contig.h ContigNode.h
ContigIterator.o: ../utils/Kmer.h ../utils/JuncPairs.h ../utils/Junction.h
Mink.o: /usr/include/stdio.h /usr/include/features.h
Mink.o: /usr/include/stdc-predef.h /usr/include/libio.h
Mink.o: /usr/include/_G_config.h /usr/include/wchar.h /usr/include/stdlib.h
Mink.o: /usr/include/alloca.h /usr/include/time.h /usr/include/xlocale.h
Mink.o: /usr/include/string.h /usr/include/unistd.h /usr/include/getopt.h
Mink.o: /usr/include/inttypes.h /usr/include/stdint.h ../utils/Bloom.h
Mink.o: ../utils/Kmer.h ReadScanner.h ../utils/JChecker.h
Mink.o: ../utils/JunctionMap.h ../utils/Junction.h ../utils/ReadKmer.h
Mink.o: ../utils/Kmer.h ../utils/DoubleKmer.h ../utils/Cap.h
Mink.o: ../utils/JuncPairs.h
