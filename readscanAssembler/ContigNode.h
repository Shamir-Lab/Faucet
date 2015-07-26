#ifndef CONTIGNODE
#define CONTIGNODE

#include <iostream>
#include "Node.h"
//#include "Contig.h"
using std::ofstream;

class Contig; // forward declare to avoid circ. dependency

class ContigNode{
public:
    unsigned char cov[5];
    Contig * contigs[5];

    ContigNode(Node node);
    ContigNode();

};

#endif