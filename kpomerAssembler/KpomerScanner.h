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
#include "../utils/Junction.h"
#include "../utils/JunctionMap.h"

#define DEBUGE(a)  //printf a

struct junction{
    unsigned char ext[4];
    unsigned char cov[4];
    unsigned char back;
};

class KpomerScanner{

private:
    Bloom* bloom;
    string kmers_file;
    JunctionMap* junctionMap;

    uint64_t NbCandKmer, NbRawCandKmer, NbJCheckKmer, NbNoJuncs, 
        NbSkipped, NbProcessed, readsProcessed, NbSolidKmer, NbSpacers;

        JChecker* jchecker;

public:
    JunctionMap*  getJunctionMap();
    void scanKpomer();
    void printScanSummary();
    void scan_kpomer(string kpomer);
    int scan_all_kpomers(string kmers_file);
    KpomerScanner(string readFile, Bloom* bloom, JChecker * jchecker);

};
#endif