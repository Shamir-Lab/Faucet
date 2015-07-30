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
    const ContigNode * node1_p; // lex. lesser end k-mer
    const ContigNode * node2_p;

    int ind1;
    int ind2;
    std::string * seq_p;
    // Contig(const ContigNode* n1, int ind1, const ContigNode* n2, int ind2, const std::string& seq);
    // Contig(const Contig& contig);
    void setEnds(const ContigNode* n1, int i1, const ContigNode* n2, int i2);
    void setSeq(const std::string& cont);
    Contig();

};

#endif