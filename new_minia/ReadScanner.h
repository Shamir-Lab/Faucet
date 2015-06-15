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

#include "Bloom.h"
#include "Kmer.h"

#define DEBUGE(a)  //printf a

struct junction{
    unsigned char ext[4];
    unsigned char cov[4];
    unsigned char back;
};

class ReadScanner{

private:
    int numPathsOut(junction j);
    int getNumComplexJunctions();
    void allocateJunctionMap(uint64_t size);
    void writeJunction(ofstream* jFile, junction toPrint);
    void updateJunction(junction* junc, int nucExt, int lengthFor, int lengthBack);
    void createJunction(kmer_type kmer, int nucExt);
    int j;
    Bloom* bloom;
    int spacerDist;
    set<kmer_type> jcheckedSet;
    unordered_map<kmer_type,junction> junctionMap;
    set<kmer_type> nextRealSet;
    unordered_map<kmer_type, junction>::iterator juncIt;
    junction * juncInfo;
    uint64_t hash0, hash1, nextHash0, nextHash1;

    int readLength;
    string reads_file;

    uint64_t ** lastHashes;
    uint64_t ** nextHashes;
    uint64_t ** tempor;


    uint64_t NbCandKmer, NbRawCandKmer, NbJCheckKmer, NbNoJuncs, 
        NbSkipped, NbProcessed, readsProcessed, NbSolidKmer, NbSpacers;

public:
    unordered_map<kmer_type,junction>  getJunctionMap();
    void smart_traverse_read(string read);
    bool find_next_junction(int* pos, kmer_type * kmer, string read);
    bool jcheck(char* kmerSeq, uint64_t nextH0, uint64_t nextH1);
    void setJ(int j);
    void junctionMapToFile(string filename);
    void scanReads( int genome_size);
    void printScanSummary();
    ReadScanner(string readFile, Bloom* bloom);

};
#endif