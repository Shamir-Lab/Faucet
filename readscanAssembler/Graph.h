#ifndef GRAPH
#define GRAPH

#include <unordered_map>
#include <unordered_set>
#include <string>
#include "../utils/Kmer.h"
#include "../utils/Junction.h"
#include "../utils/Bloom.h"
#include "../utils/JunctionMap.h"
#include "../utils/JChecker.h"
#include "Node.h"
using std::unordered_map;
using std::unordered_set;
using std::string;

//Stores all the info involved in moving from one node to another in the graph 
struct GraphSearchResult{
    kmer_type kmer;
    bool isNode; //either node or sink
    int index; //index from it that points back to the start point
    int distance; //how far away it was
};

//Stores all the info involved in moving from one node to another in the graph by doing a bloom scan
struct BfSearchResult{
    kmer_type kmer;
    bool isNode; //either node or sink
    int index; //index from it that points back to the start point
    int distance; //how far away it was
    string contig; //the contig!
};

class Graph
{
private:
    string contigFile; //just the file name
    unordered_map<kmer_type, Node> nodeMap; //maps kmers to nodes
    unordered_map<kmer_type, int>* realExtensions; //acts like a FP set- maps kmers at FP junctions to their one and only real extension
    unordered_set<kmer_type>* sinks; //set of sink kmers
    Bloom* bloom;
    JChecker* jchecker;
    Node * getNode(kmer_type kmer);

    //Assumes cFPs and sinks were already handled, only complex junctions should remain.  
    //Gets junctions from juncMap and builds equivalent Nodes for graph. Kills junctions as it goes.
    void getNodesFromJunctions(JunctionMap* juncMap); 

    //Uses the bloom filter to move from the given node to the next node, along the given index.
    //Useful for linking nodes initially and generating contigs
    BfSearchResult findNeighborBf(Node node, kmer_type startKmer , int index);

    //Uses only the graph to move from the given node to the next node, along the given index
    //Useful for graph simplification and assembly traverrsal
    GraphSearchResult findNeighborGraph(Node node, kmer_type startKmer , int index);

    bool isSink(kmer_type kmer);
    bool isNode(kmer_type kmer);
    bool isRealExtension(kmer_type kmer, int ext);

    //Gets the valid extension of the given kmer based on the bloom filter and cFPs.  JChecks! so this cuts off tips
    //Assume its not a junction
    //Returns -1 if there is no valid extension
    //Assumes full cFP set and uses it as reference
    int getValidJExtension(DoubleKmer kmer, int dist, int max);

    //puts in distances and junction IDs to link the two nodes
    void directLinkNodes(kmer_type kmer1, int index1, kmer_type kmer2, int index2, int distance);
    void linkNodeToSink(kmer_type nodeKmer, int index, kmer_type sinkKmer, int distance);

    //for linking a node to a neighbor found in a bloom scan
    void linkNeighbor(kmer_type startKmer, int index, BfSearchResult result);

    //Traverses every contig using the bloom filter, in order to either link the nodes to each other or print the contigs
    void traverseContigs(bool linkNodes, bool printContigs);
    
    //Returns true if the node at that index leads to a sink, and the path to the sink is shorter than maxTipLength
    bool isTip(Node node, int index, int maxTipLength);
    int getNumTips(Node node, int maxTipLength);

    //deletes the node associated with the given kmer, and replaces it with a realExtension entry
    //Prints an error if the node has multiple valid paths or if the kmer doesn't correspond to a node
    void deleteNode(kmer_type kmer);

    //returns the index from startNode that leads to the destination kmer
    int getPathIndex(Node startNode, kmer_type destinationKmer);

public: 
    int cutTips(int maxTipLength);   //remove all short tips, return number cut
    void linkNodes(); //use BF to link nodes
    void printContigs(string filename); 
    void buildGraph(JunctionMap* juncMap); //makes nodes out of complex junctions, replaces others with real extensions
    void printGraph(string fileName);
    Graph(Bloom* bloom, JChecker* jcheck);
};

#endif
