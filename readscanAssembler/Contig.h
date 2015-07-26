#ifndef CONTIG
#define CONTIG

#include <iostream>
#include "ContigNode.h"
using std::ofstream;

class Contig{ // modeled after node implementation
public:
    // coverage can be obtained by end nodes
    // length from sequence
    ContigNode * node1; // lex. lesser end k-mer
    ContigNode * node2;
    char ind1;
    char ind2;
    string seq;
    Contig(ContigNode* n1, char ind1, ContigNode* n2, char ind2);
    Contig();

};

#endif