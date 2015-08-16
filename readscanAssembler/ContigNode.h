#ifndef CONTIGNODE
#define CONTIGNODE

 class Contig; // forward declare to avoid circ. dependency
// following http://www.cplusplus.com/forum/articles/10627/#msg49679

#include <iostream>
#include <vector>
#include "Node.h"
#include "Contig.h"
#include "../utils/Kmer.h"
using std::ofstream;


class ContigNode{
    
    unsigned char cov[4];
public:
    Contig * contigs[5];

    ContigNode(Junction junction);
    ContigNode(Node node);
    ContigNode();
    int numPathsOut();
    int getCoverage(int nucExt);
    void setCoverage(Junction junc);
    void update(int nucExt, Contig * contig);

    //removes the given path out of this node.
    //Removes contig pointer, set coverage to 0
    void breakPath(int nucExt);

    //for traversal
    bool hasNeighbor(int index);
    ContigNode* getNeighbor(int index);
    std::string getString();
    kmer_type getKmer();
};

#endif