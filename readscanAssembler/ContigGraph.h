#ifndef CONTIG_GRAPH
#define CONTIG_GRAPH

class ContigGraph; //forward declare

#include <fstream>
#include <iostream> 
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <vector>
#include "../utils/Kmer.h"
#include "../utils/Junction.h"
#include "../utils/Bloom.h"
#include "../utils/JunctionMap.h"
#include "../utils/JChecker.h"
#include "Contig.h"
#include "ContigNode.h"
#include "BfSearchResult.h"
#include "ContigIterator.h"
using std::unordered_set;
using std::unordered_map;
using std::string;

class ContigGraph 
{
std::vector<ContigNode> nodeVector;
std::vector<Contig> isolated_contigs;
unordered_map<kmer_type,ContigNode> nodeMap;

public: 
    unordered_map<kmer_type, ContigNode> * getNodeMap();
    void switchToNodeVector();

    //Returns a list of kmers that could be in a junction pair that would help to
    //disentangle a contig.  
    //Specifically all recorded branch junctions on the contig on index index of node node
    std::list<kmer_type> getPairCandidates(ContigNode* node, int index);

    //Gets tail bound for binomia ldistribution with number of trials, probability,
    //and for result specified
    double getTailBound(int numTrials, double p, int result);

    //Gets number of supporting pairs given candidate list
    //TODO: normalize by expected FP rate of filter
    double getScore(std::list<kmer_type> leftCand, std::list<kmer_type> rightCand, Bloom* pair_filter, double fpRate);
   

    //a,b are on backNode, c,d are on forwardNode
    //a pairs with c, b pairs with d
    //Does not go ahead with the operation if degeneracies are detected
    //Returns true if it goes ahead with disentanglement
    bool disentanglePair(Contig* contig, ContigNode* backNode, ContigNode* forwardNode, int a, int b, int c, int d);
    void addIsolatedContig(Contig contig); 
    bool isErrorContig(Contig* contig);
    void deleteContig(Contig* contig);
    bool cleanGraph(); //Cleans graph and returns true if any changes were made

    void checkGraph();
    void printContigFastG(ofstream* fastgFile, Contig * contig);
    void printContigs(string filename); //prints the contigs raw
     void printGraph(string fileName); //prints graph : TBD print format- fastg?
    ContigGraph();

    //Creates a contig node if it doesn't already exist
    //If it exists, does nothing and returns the existing one.
    //Otherwise, returns the new one
    ContigNode * createContigNode(kmer_type kmer, Junction junction);    
    int disentangle(Bloom* pair_filter);
private:
    int deleteErrorContigs();   //remove tips, chimeras, and bubbles. Return number of deleted contigs.
    int collapseDummyNodes(); //removes nodes with only one real extension, merges forward and back contigs
    int destroyDegenerateNodes();// Removes nodes with no back contig or no forward contigs

    unordered_map<kmer_type, ContigNode> contigNodeMap; // maps kmers to ContigNodes after contigs constructed
    void collapseNode(ContigNode * node);
    void cutPath(ContigNode* node, int index); //used on nodes with no backward contig
};

#endif
