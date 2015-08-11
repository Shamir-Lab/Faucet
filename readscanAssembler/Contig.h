#ifndef CONTIG
#define CONTIG

class ContigNode; // forward declaration

#include <iostream>
#include <string.h>

using std::ofstream;

class Contig{
public:
    // coverage can be obtained by end nodes
    // length from sequence
    ContigNode * node1_p; 
    ContigNode * node2_p;

    ContigNode* otherEndNode(ContigNode * oneEnd);//returns a pointer to the node at the other end

    int ind1;
    int ind2;
    std::string * seq_p;
    void setEnds(ContigNode* n1, int i1, ContigNode* n2, int i2);
    void setSeq(const std::string& cont);
    int getMinIndex();
    Contig();

};

#endif