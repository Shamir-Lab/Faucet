#ifndef JUNCTION_MAP
#define JUNCTION_MAP

#include <unordered_map>
#include <string>
#include "Kmer.h"
#include "Junction.h"
#include "Cap.h"

using std::unordered_map;
using std::string;

class JunctionMap{

private: 
    unordered_map<kmer_type,Junction> junctionMap;
    unordered_map<kmer_type,Cap> capMap;

public:

    void linkJunctions(kmer_type kmer1, int ext1, kmer_type kmer2, int ext2, int dist);
    int getNumComplexJunctions();
    int getNumCaps();
    void createJunction(kmer_type kmer);
    void addCap(kmer_type kmer, Cap cap);
    void writeToFile(string filename);
    int getNumJunctions();
    bool contains(kmer_type kmer);
    Junction* getJunction(kmer_type kmer);
    Cap* getCap(kmer_type kmer);
    JunctionMap();

};
#endif