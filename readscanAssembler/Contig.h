#ifndef CONTIG
#define CONTIG

#include <iostream>
#include "../utils/Node.h"
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
    Contig(Node* n1, char ind1, Node* n2, char ind2);
    Contig();

};

#endif