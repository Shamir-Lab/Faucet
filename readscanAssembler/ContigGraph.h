#ifndef CONTIG_GRAPH
#define CONTIG_GRAPH

class ContigGraph; //forward declare
#include <ostream>
#include <fstream>
#include <iostream> 
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <vector>
#include "../utils/Kmer.h"
#include "../utils/Junction.h"
#include "../utils/BloomFilter.hpp"
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
#include "../utils/sparsepp.h"
using std::string;
using spp::sparse_hash_map;


class ContigGraph 
{
std::vector<ContigNode> nodeVector;
std::vector<Contig> isolated_contigs;
unordered_map<kmer_type,ContigNode> nodeMap;
// sparse_hash_map<kmer_type,ContigNode> nodeMap;
unordered_map<kmer_type,ContigNode>::iterator it;
// sparse_hash_map<kmer_type,ContigNode>::iterator it;


int read_length;

public: 
    std::vector<Contig> * getIsolatedContigs();
    unordered_map<kmer_type, ContigNode> * getNodeMap();
    // sparse_hash_map<kmer_type, ContigNode> * getNodeMap();

    void switchToNodeVector();
    void setReadLength(int length);

    //Gets tail bound for binomia ldistribution with number of trials, probability,
    //and for result specified
    double getTailBound(int numTrials, double p, int result);
    bool isCollapsible(ContigNode * node);

    //Gets number of supporting pairs given candidate list
    //TODO: normalize by expected FP rate of filter
    double getScore(std::list<JuncResult> leftCand, std::list<JuncResult> rightCand, Bloom* pair_filter, double fpRate, int insertSize);
    std::pair <Contig*,Contig*> getMinMaxForwardExtensions(ContigNode * node, std::string trait);


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
    bool cleanGraph(Bloom* short_pair_filter, Bloom* long_pair_filter); //Cleans graph and returns true if any changes were made

    bool checkGraph();
    void printContigFastG(std::ostream* fastgFile, Contig * contig);
    void printContigs(string filename); //prints the contigs raw
    void printGraph(string fileName); //prints graph : TBD print format- fastg?
    ContigGraph();

    Contig* getLongestContig();

    //Creates a contig node if it doesn't already exist
    //If it exists, does nothing and returns the existing one.
    //Otherwise, returns the new one
    ContigNode * createContigNode(kmer_type kmer, Junction junction);    
    int disentangleParallelPaths(Bloom* pair_filter, double insertSize, double std);
    int disentangleLoopPaths(Bloom* pair_filter, double insertSize, double std);
    int removeChimericExtensions(int insertSize);
    int collapseBulges(int max_dist);
    bool deleteTipsAndClean();
    bool breakPathsAndClean();
    bool disentangleAndClean(Bloom* pair_filter, double insertSize, double std);
    bool areEquivalentContigCoverages(ContigJuncList A, ContigJuncList B, double frac);
    bool areDifferentialContigCoverages(ContigJuncList A, ContigJuncList B);

private:
    int deleteTips();
    int deleteIsolatedContigs();
    bool testAndCutIfDegenerate(ContigNode* node);
    int collapseDummyNodes(); //removes nodes with only one real extension, merges forward and back contigs
    int destroyDegenerateNodes();// Removes nodes with no back contig or no forward contigs
    // int cutIfDegenerate(ContigNode* node, kmer_type kmer, auto it);

    unordered_map<kmer_type, ContigNode> contigNodeMap; // maps kmers to ContigNodes after contigs constructed
    // sparse_hash_map<kmer_type, ContigNode> contigNodeMap; // maps kmers to ContigNodes after contigs constructed

    void collapseNode(ContigNode * node, kmer_type kmer);
    void cutPath(ContigNode* node, int index); //used on nodes with no backward contig
};

// class NodeContigData{ //contains all info for an entry in the queue for a node BFS
// public:
//     Contig * contig;
//     ContigNode * far_node;
//     int len;
//     double cov;
//     // std::list<JuncResult> juncs;

//     NodeContigData(ContigNode * node, int index, int insertSize);
//     NodeContigData();

// };

#endif
