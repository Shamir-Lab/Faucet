#ifndef JUNCTION_MAP
#define JUNCTION_MAP

#include <unordered_map>
#include <set>
#include <string>
#include "Kmer.h"
#include "Junction.h"
#include "Cap.h"
#include "ReadKmer.h"
#include "Bloom.h"
#include <fstream>
using std::ofstream;
using std::unordered_map;
using std::string;
using std::set;

class JunctionMap{

private: 
    Bloom* bloom;

public:
    unordered_map<kmer_type,Junction> junctionMap;  
    set<kmer_type> getCFPs();// returns the set of all cFPs. Destroys junctions that are no longer needed along the way.
    set<kmer_type> getSinks(); //returns a set of all sinks as determined by scanning the bloom filter
    int getSkipDist(ReadKmer* readKmer, bool direction);
    void directLinkJunctions(ReadKmer* kmer1, ReadKmer* kmer2);//directly links two junctions on the same read
    int getNumComplexJunctions();
    int getNumCaps();
    void createJunction(kmer_type kmer);
    void createJunction(ReadKmer* readKmer);
    void writeToFile(string filename);
    int getNumJunctions();
    bool contains(kmer_type kmer);
    bool contains(ReadKmer* readKmer);
    Junction* getJunction(ReadKmer kmer);
    Junction* getJunction(kmer_type kmer);
    kmer_type * findSink(Junction junc, kmer_type kmer ,int i);
    void killJunction(kmer_type kmer);
    void finishGraphWriteBasicContigs();
    JunctionMap(Bloom* bloo);
};
#endif