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
#include "../utils/DoubleKmer.h"
#include "../utils/Cap.h"
#include "../utils/DoubleKmer.h"

#define DEBUGE(a)  //printf a

class ReadScanner{

private:
    int j;
    Bloom* bloom;
    kmer_type firstReadJunc, lastReadJunc; //stores kmers of last and first junction in a read from the read scan.  For use by spacers.
    //int spacerDist;
    set<kmer_type> jcheckedSet;
    set<kmer_type> nextRealSet;
    uint64_t hash0, hash1,
    nextHash0, nextHash1;
    int readLength;
    string reads_file;

    uint64_t NbCandKmer, NbRawCandKmer, NbJCheckKmer, NbNoJuncs, 
        NbSkipped, NbProcessed, readsProcessed, NbSolidKmer, NbSpacers;

    JChecker* jchecker;
    JunctionMap* junctionMap;
public:
    void resetHashes(kmer_type kmer);//for testing
    JunctionMap* getJunctionMap();
    void setJ(int j);
    void scanReads();
    void printScanSummary();
    
    bool testForJunction(DoubleKmer kmer);
    void scan_forward(string read); 
    bool find_next_junction(DoubleKmer * kmer);//adjusts position and kmer and returns the junction

    ReadScanner(string readFile, Bloom* bloom, JChecker* jchecker);
};
#endif