#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <inttypes.h>
#include <cmath> // for log2f
#include <algorithm> // for max
#include <unistd.h> // for truncate
#include <string>
#include <unordered_map>
#include <set>
using std::string;
using std::unordered_map;
using std::set;
using std::ofstream;

#ifndef READSCAN_H
#define READSCAN_H

#include "../utils/Bloom.h"
#include "../utils/Kmer.h"
#include "../utils/JChecker.h"
#include "../utils/JunctionMap.h"
#include "../utils/Junction.h"

#define DEBUGE(a)  //printf a

class ReadScanner{

private:
    void allocateJunctionMap(uint64_t size);
    int j;
    Bloom* bloom;
    int spacerDist;
    set<kmer_type> jcheckedSet;
    set<kmer_type> nextRealSet;
    Junction * juncInfo;
    uint64_t hash0, hash1, nextHash0, nextHash1;

    int readLength;
    string reads_file;

    uint64_t NbCandKmer, NbRawCandKmer, NbJCheckKmer, NbNoJuncs, 
        NbSkipped, NbProcessed, readsProcessed, NbSolidKmer, NbSpacers;

    JChecker* jchecker;
    JunctionMap* junctionMap;
public:
    JunctionMap* getJunctionMap();
    void smart_traverse_read(string read);
    bool find_next_junction(int* pos, kmer_type * kmer, string read);
    void setJ(int j);
    void scanReads();
    void printScanSummary();
    ReadScanner(string readFile, Bloom* bloom, JChecker* jchecker);

};
#endif