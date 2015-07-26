#ifndef CONTIGNODE
#define CONTIGNODE

#include <iostream>
#include "Node.h"
using std::ofstream;

class ContigNode{
public:
    unsigned char cov[5];
    Contig * contigs[5];

    ContigNode(Node node);
    ContigNode();

};

#endif