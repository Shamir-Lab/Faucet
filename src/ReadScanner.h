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
// #include <unordered_map>
#include <set>
using std::string;
// using std::unordered_map;
using std::set;
using std::ofstream;

#ifndef READSCAN_H
#define READSCAN_H

#include "../utils/Bloom.h"
#include "../utils/Kmer.h"
#include "../utils/JChecker.h"
#include "../utils/JunctionMap.h"
#include "../utils/Junction.h"
#include "../utils/ReadKmer.h"
#include "../utils/Cap.h"
#include "../utils/JuncPairs.h"

#define DEBUGE(a)  //printf a

class ReadScanner{

private:
    int maxSpacerDist; // maximum distance in bases between two junctions (spacers added to bridge gaps)
    Bloom* bloom;
    Bloom* short_pair_filter;
    Bloom* long_pair_filter;
    set<kmer_type> jcheckedSet;
    set<kmer_type> nextRealSet;
    set<std::pair<kmer_type, kmer_type> > juncPairSet;
    set<kmer_type> backwardSet;
    uint64_t hash0, hash1,
    nextHash0, nextHash1;
    string reads_file;

    uint64_t NbCandKmer, NbRawCandKmer, NbJCheckKmer, NbNoJuncs, 
        NbSkipped, NbProcessed, readsProcessed, NbSolidKmer,readsNoErrors,
         NbJuncPairs, unambiguousReads;

    JChecker* jchecker;
    JunctionMap* junctionMap;

    //Should only be called on a read with no real junctions
    //Adds a fake junction in the middle and points it to the two ends.  This ensures we have coverage of long linear regions, and that we capture
    //sinks at the end of such regions.
    kmer_type add_fake_junction(string read);
    
public:
    JunctionMap* getJunctionMap();

    //Scans one input read; breaks into small segments and calls scan_forward
    //Returns back junctions along read from beginning to end
    std::list<kmer_type> scanInputRead(string read, bool no_cleaning);

    void scanReads(bool fastq, bool paired_ends, bool no_cleaning); //scans all the reads.  Fastq if fastq, otherwise fasta
    void printScanSummary(); //prints statistics from the readscan
    
    //Determines if the given ReadKmer is a junction.
    //If it's on the middle of the read, just verifies alternate extensions.
    //Special logic is needed to handle kmers that are near the ends, if j is not 0, to ensure that 
    //the real extension seen on the read represents a valid, jcheckable option, and not a tip shorter than j.
    bool testForJunction(ReadKmer kmer);

    //Starting from the given kmer, scans forward until a junction is found or the end of the read is hit.
    //Returns true if a junction was found.  The supplied ReadKmer is also adjusted to the position of the new junction,\
    //or to the end of the read.
    bool find_next_junction(ReadKmer * kmer, int lastJuncPos);

    //Returns substrings of the read that are valid with BF and longer than sizeKmer
    std::list<string> getValidReads(string read);

    //Scans a read. 
    //Identifies all junctions on the read, and links adjacent junctions to each other.
    //Also updates the relevant distance field on the first junction to point to the start of the read, and on the last
    //Junction to point to the end of the read.
    //If there are no junctions, add_fake_junction is called
    //Returns back junctions along the read from beginning to end
    std::list<kmer_type> scan_forward(string read, bool no_cleaning); 

    ReadScanner(JunctionMap* juncMap, string readFile, Bloom* bloom, Bloom* short_pair_filter, Bloom* long_pair_filter, JChecker* jchecker, int maxSpacerDist);
};
#endif