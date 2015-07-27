#ifndef CONTIGNODE
#define CONTIGNODE

class Contig; // forward declare to avoid circ. dependency
// following http://www.cplusplus.com/forum/articles/10627/#msg49679

#include <iostream>
#include "Node.h"
//#include "Contig.h"
using std::ofstream;


class ContigNode{
public:
    unsigned char cov[5];
    Contig * contigs[5];

    ContigNode(string cstr);
    ContigNode();

};

#endif