#ifndef JUNCTION_MAP
#define JUNCTION_MAP

#include <unordered_map>
#include <string>
#include "Kmer.h"
#include "Junction.h"
#include "Cap.h"
#include "DoubleKmer.h"

using std::unordered_map;
using std::string;

class JunctionMap{

private: 
    unordered_map<kmer_type,Junction> junctionMap;
    unordered_map<kmer_type,Cap> capMap;

public:
    int getSkipDist(DoubleKmer* doubleKmer, bool direction);
    void directLinkJunctions(DoubleKmer* kmer1, DoubleKmer* kmer2);//directly links two junctions on the same read
    void extendJunctionCap(DoubleKmer* kmer, bool dir);//extends the junction at the kmer to cap at the end of the read.  Deletes old cap.
    void superLinkJunctionToCap(DoubleKmer kmer, Cap cap);//For when a cap points the other way and the junction needs to be pointed to the other end
    void superLinkCapToCap(Cap cap1, Cap cap2);//when two caps point away from each other and need to be merged
    int getNumComplexJunctions();
    int getNumCaps();
    void createJunction(kmer_type kmer);
    void writeToFile(string filename);
    int getNumJunctions();
    bool contains(kmer_type kmer);
    Junction* getJunction(kmer_type kmer);
    JunctionMap();
};
#endif