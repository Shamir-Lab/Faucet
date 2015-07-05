#ifndef JUNCTION_MAP
#define JUNCTION_MAP

#include <unordered_map>
#include <string>
#include "Kmer.h"
#include "Junction.h"
#include "Cap.h"
#include "DoubleKmer.h"
#include "Bloom.h"
#include <fstream>
using std::ofstream;
using std::unordered_map;
using std::string;

class JunctionMap{

private: 
    Bloom* bloom;

public:
    unordered_map<kmer_type,Junction> junctionMap;
    int getSkipDist(DoubleKmer* doubleKmer, bool direction);
    void directLinkJunctions(DoubleKmer* kmer1, DoubleKmer* kmer2);//directly links two junctions on the same read
    int getNumComplexJunctions();
    int getNumCaps();
    void createJunction(kmer_type kmer);
    void createJunction(DoubleKmer* doubleKmer);
    void writeToFile(string filename);
    int getNumJunctions();
    bool contains(kmer_type kmer);
    bool contains(DoubleKmer* doubleKmer);
    Junction* getJunction(DoubleKmer kmer);
    Junction* getJunction(kmer_type kmer);
    void killJunction(kmer_type kmer);
    void finishGraphWriteBasicContigs();
    JunctionMap(Bloom* bloo);
};
#endif