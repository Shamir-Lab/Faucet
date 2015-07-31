#ifndef CONTIGNODE
#define CONTIGNODE

// class Contig; // forward declare to avoid circ. dependency
// following http://www.cplusplus.com/forum/articles/10627/#msg49679

#include <iostream>
#include <vector>
#include "Node.h"
#include "Contig.h"
using std::ofstream;


class ContigNode{
public:
    unsigned char cov[5];
    // std::vector<Contig*> contigs;
    Contig * contigs[5];

    ContigNode(Node node);
    ContigNode();
    void update(int nucExt, Contig * contig);

};

#endif