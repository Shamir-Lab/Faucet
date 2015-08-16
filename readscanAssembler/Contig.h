#ifndef CONTIG
#define CONTIG

class ContigNode; // forward declaration

#include "../utils/Kmer.h"
#include "ContigNode.h"
#include <iostream>
#include <string.h>

using std::ofstream;

class Contig{
private:
    //utility for linking if they're both facing forward
    //Glues end 2 of this contig to end 1 of the other
    //Doesn't change the value of either contig- just returns the concatenation.
    Contig* concatenate(Contig* otherContig);

public:
    // coverage can be obtained by end nodes
    // length from sequence
    ContigNode * node1_p; 
    ContigNode * node2_p;
    int ind1;
    int ind2;
    std::string seq;
    std::list<unsigned char> juncDistances;
    int coverageSum;

    //Concatenates the two contigs, gluing together the specified sides
    Contig* concatenate(Contig* otherContig, int thisSide, int otherSide);

    void reverse(); //reverses the contig orientation
    ContigNode* otherEndNode(ContigNode * oneEnd);//returns a pointer to the node at the other end
    void setEnds(ContigNode* n1, int i1, ContigNode* n2, int i2);
    void setIndices(int i1, int i2);
    void setSeq(std::string cont);
    void setJuncDistances(std::list<unsigned char> juncDists);
    void setCoverage(int cov);
    void addCoverage(int cov);
    float getAvgCoverage();
    int getMinIndex();
    kmer_type getNodeKmer(ContigNode * contigNode);    //Assumes the given contig node points to one end of this contig
    kmer_type getSideKmer(int side);    //either 1 or 2
    int getSide(ContigNode* node);
    std::string getStringRep();
    Contig();

};

#endif