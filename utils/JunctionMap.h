#ifndef JUNCTION_MAP
#define JUNCTION_MAP

#include <unordered_map>
#include <string>
#include "Kmer.h"
#include "Junction.h"

using std::unordered_map;
using std::string;

class JunctionMap{

private: 
    unordered_map<kmer_type,Junction> junctionMap;

public:

    int getNumComplexJunctions();
    void createJunction(kmer_type kmer, int nucExt);
    void writeToFile(string filename);
    int getNumJunctions();
    bool contains(kmer_type kmer);
    Junction* getJunction(kmer_type kmer);
    JunctionMap();

};
#endif