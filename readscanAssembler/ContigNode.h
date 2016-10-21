#ifndef CONTIGNODE
#define CONTIGNODE

class Contig; // forward declare to avoid circ. dependency
// following http://www.cplusplus.com/forum/articles/10627/#msg49679

#include <iostream>
#include <vector>
#include <set>
#include <queue>
#include <unordered_map>
#include <list>
#include "Contig.h"
#include "../utils/Kmer.h"
#include "../utils/JuncPairs.h"
#include "../utils/Junction.h" 
using std::ofstream;


class ContigNode{
    
public:
    unsigned char cov[4];
    Contig * contigs[5];

    ContigNode(Junction junction);
    ContigNode();
    bool isInvertedRepeatNode();


    bool checkValidity();
    // returns 0 if target node not reached at up to max_dist, otherwise returns distance on branching path
    std::list<Contig*> doPathsConvergeNearby(int max_ind, int min_ind, int max_dist);


    //gets the neighbors of the specified contig- if contigIndex is 4, returns all forward neighbors
    //If contig index isn't 4, only returns back contig as a neighbor
    std::vector<std::pair<Contig*, bool> > getFastGNeighbors(int contigIndex);



    std::list<JuncResult> getPairCandidates(int index, int maxDist);
      
    void replaceContig(Contig* oldContig, Contig* newContig);
    int numPathsOut();
    int indexOf(Contig* contig);
    kmer_type getForwardExtension(int index);
    std::vector<int> getIndicesOut();
    int getCoverage(int nucExt);
    int getTotalCoverage();//returns getCoverage(4)
    void setCoverage(Junction junc);
    void setCoverage(int nucExt, int coverage);
    void update(int nucExt, Contig * contig);
    kmer_type getUniqueKmer(int index);//returns base kmer for backward, or extension for forward index

    //removes the given path out of this node.
    //Removes contig pointer, set coverage to 0
    void breakPath(int nucExt);
    void clearNode();

    //for traversal
    bool hasNeighbor(int index);
    ContigNode* getNeighbor(int index);
    std::string getString();
    kmer_type getKmer();
};

class NodeQueueEntry{ //contains all info for an entry in the queue for a node BFS
public:
    ContigNode* node;
    int index;
    int startDist;

    NodeQueueEntry(ContigNode* n, int i, int s);
    NodeQueueEntry();

    std::list<JuncResult> getJuncResults(int m); //returns immediate junc results from contig along this index

    void addNeighbors(std::deque<NodeQueueEntry> & queue); // , bool to_back); //searches forward one step, adds relevant nodes to the queue
    void recordParents(std::unordered_map<NodeQueueEntry, NodeQueueEntry>& parents);
    std::list<Contig*> reconstructPathFromParents(std::unordered_map<NodeQueueEntry, NodeQueueEntry>& parents);
    friend bool operator==(NodeQueueEntry a, NodeQueueEntry b) { 
        return a.node == b.node && a.index == b.index && a.startDist == b.startDist; 
    };

};


// struct MyHash {
//   size_t operator()(const NodeQueueEntry& x) const { return std::hash<uint64_t>()(x.node->getUniqueKmer(x.index)); }
// };
namespace std {
  template <> struct hash<NodeQueueEntry>
  {
    size_t operator()(const NodeQueueEntry & x) const
    {
        return std::hash<uint64_t>()(x.node->getUniqueKmer(x.index)); }
    };
}
#endif