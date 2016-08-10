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
#include "../utils/JuncPairs.h"
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
unordered_map<kmer_type,ContigNode>::iterator it;
// unordered_set<kmer_type> expiredKmers;

int read_length;

public: 
    unordered_map<kmer_type, ContigNode> * getNodeMap();
    void switchToNodeVector();
    void setReadLength(int length);

    //Gets tail bound for binomia ldistribution with number of trials, probability,
    //and for result specified
    double getTailBound(int numTrials, double p, int result);

    //Gets number of supporting pairs given candidate list
    //TODO: normalize by expected FP rate of filter
    double getScore(std::list<JuncResult> leftCand, std::list<JuncResult> rightCand, Bloom* pair_filter, double fpRate, int insertSize);
   

    //a,b are on backNode, c,d are on forwardNode
    //a pairs with c, b pairs with d
    //Does not go ahead with the operation if degeneracies are detected
    //Returns true if it goes ahead with disentanglement
    void disentanglePair(Contig* contig, ContigNode* backNode, ContigNode* forwardNode, int a, int b, int c, int d);
    void addIsolatedContig(Contig contig);
    std::vector<int> getUnsupportedExtensions(ContigNode* node, Bloom* pair_filter, int insertSize);
    bool isLowCovContig(Contig* contig);
    bool isLowMassContig(Contig* contig);
    bool isTip(ContigNode* node, int i);
    bool isBubble(ContigNode* node);
    std::list<Contig*> getPathIfSimpleBulge(ContigNode* node, int max_dist);

    void deleteContig(Contig* contig);
    bool cleanGraph(Bloom* short_pair_filter, Bloom* long_pair_filter, int insertSize); //Cleans graph and returns true if any changes were made

    bool checkGraph();
    void printContigFastG(ofstream* fastgFile, Contig * contig);
    void printContigs(string filename); //prints the contigs raw
    void printGraph(string fileName); //prints graph : TBD print format- fastg?
    ContigGraph();

    Contig* getLongestContig();

    //Creates a contig node if it doesn't already exist
    //If it exists, does nothing and returns the existing one.
    //Otherwise, returns the new one
    ContigNode * createContigNode(kmer_type kmer, Junction junction);    
    int disentangle(Bloom* pair_filter, int insertSize);
    int removeChimericExtensions(int insertSize);
    int collapseBulges(int max_dist);
    bool deleteTipsAndClean();
    bool breakPathsAndClean(Bloom* pair_filter, int insertSize);
    bool disentangleAndClean(Bloom* pair_filter, int insertSize);
    bool areEquivalentContigCoverages(Contig* contig_a, Contig* contig_b, 
        ContigNode * node_a, ContigNode * node_b, double frac, int insertSize);

private:
    int popBubblesByCoverageRatio();
    int deleteTipsAndLowCoverageContigs();   //remove tips, chimeras, and bubbles. Return number of deleted contigs.
    int deleteTips();
    int deleteIsolatedContigs();
    bool testAndCutIfDegenerate(ContigNode* node);
    int breakUnsupportedPaths(Bloom* pair_filter, int insertSize); //removes extensions of junctions not supported by paired ends
    int collapseDummyNodes(); //removes nodes with only one real extension, merges forward and back contigs
    int destroyDegenerateNodes();// Removes nodes with no back contig or no forward contigs
    // int cutIfDegenerate(ContigNode* node, kmer_type kmer, auto it);

    unordered_map<kmer_type, ContigNode> contigNodeMap; // maps kmers to ContigNodes after contigs constructed
    void collapseNode(ContigNode * node, kmer_type kmer);
    void cutPath(ContigNode* node, int index); //used on nodes with no backward contig
};

#endif
