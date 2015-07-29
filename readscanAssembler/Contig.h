#ifndef CONTIG
#define CONTIG

class ContigNode; // forward declaration

#include <iostream>
#include <string.h>

using std::ofstream;

class Contig{ // modeled after node implementation
public:
    // coverage can be obtained by end nodes
    // length from sequence
    ContigNode * node1; // lex. lesser end k-mer
    ContigNode * node2;
    int ind1;
    int ind2;
    std::string seq;
    Contig(ContigNode* n1, int ind1, ContigNode* n2, int ind2, std::string seq);
    Contig();

};

#endif