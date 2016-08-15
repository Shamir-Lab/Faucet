#ifndef CONTIG
#define CONTIG

#include "../utils/Kmer.h"
#include "../utils/JuncPairs.h"
#include "../utils/Bloom.h"
#include "../utils/JuncPairs.h"
#include "../utils/ContigJuncList.h"
#include "ContigNode.h"
#include <iostream>
#include <string.h>

using std::ofstream;

class ContigNode; // forward declaration

class Contig{
private:
    //utility for linking if they're both facing forward
    //Glues end 2 of this contig to end 1 of the other
    //Doesn't change the value of either contig- just returns the concatenation.
    Contig* concatenate(Contig* otherContig);

    std::vector<std::pair<Contig*, bool> > getNeighbors(bool forward);

public:
    ~Contig();
    std::string getFastGName(bool RC);
    std::string getFastGHeader(bool RC);

    // length can be obtained from sequence
    ContigNode * node1_p;  //adjacent node on side 1
    ContigNode * node2_p; //adjacent node on side 2
    int ind1; //index on which it connects to the node, on side 1
    int ind2; //index on which it connects to the node, on side 2
    
    //list of coverage and distance for interior junctions along this contig- since we can use for pair BF and coverage info
    ContigJuncList contigJuncs;

    bool checkValidity();
    bool isDegenerateLoop();//returns true if both sides have same node and same index

    //Concatenates the two contigs, gluing together the specified sides
    Contig* concatenate(Contig* otherContig, int thisSide, int otherSide);

    void printPairStatistics(Bloom* pair_filter);
    int length(); //returns length of sequence
    void reverse(); //reverses the contig orientation
    ContigNode* otherEndNode(ContigNode * oneEnd);//returns a pointer to the node at the other end
    void setEnds(ContigNode* n1, int i1, ContigNode* n2, int i2);
    void setIndices(int i1, int i2);
    void setSeq(std::string cont){contigJuncs.setSeq(cont);}
    std::string getSeq(){return contigJuncs.getSeq();}
    double getAvgCoverage();
    double getAvgCoverage(std::list<JuncResult> results);

    double getCoverageSampleVariance();
    double getCoverageSampleVariance(std::list<JuncResult> results);
    void setContigJuncs(ContigJuncList juncList){ contigJuncs = juncList;}
    std::list<JuncResult> getJuncResults(int side, int startDist, int maxDist);

    //gets the node coverages on each end, returns minimum as base line for how much this should be covered.
    //If the contig is isolated, returns 0
    int getTotalDistance(){ return contigJuncs.getTotalDistance(); }

    float getMass();
    int getMinIndex();
    kmer_type getNodeKmer(ContigNode * contigNode);    //Assumes the given contig node points to one end of this contig
    kmer_type getSideKmer(int side);    //either 1 or 2
    int getSide(ContigNode* node);
    int getSide(ContigNode* node, int index);
    ContigNode* getNode(int side);
    int getIndex(int side);
    bool isIsolated();//return true if both sides point to null
    void setSide(int side, ContigNode* node);
    std::string getStringRep();
    Contig();

};

#endif