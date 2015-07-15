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

class Graph
{
private:

    unordered_map<kmer_type, Node> nodeMap;
    unordered_map<kmer_type, int>* realExtensions;
    unordered_set<kmer_type>* sinks;
    Bloom* bloom;
    JChecker* jchecker;
    Node * getNode(kmer_type kmer);
    void getNodesFromJunctions(JunctionMap* juncMap);
    void findAndLinkNeighbor(Node node, kmer_type startKmer , int index, ofstream* cFile);
    bool isSink(kmer_type kmer);
    bool isNode(kmer_type kmer);
    bool isRealExtension(kmer_type kmer, int ext);
    int getValidJExtension(DoubleKmer kmer, int dist, int max);
    void directLinkNodes(kmer_type kmer1, int index1, kmer_type kmer2, int index2, int distance);
    void linkNodeToSink(kmer_type nodeKmer, int index, kmer_type sinkKmer, int distance);

public: 
    void buildGraph(JunctionMap* juncMap);
    void linkNodesPrintContigs(string fileName);
    void printGraph(string fileName);
    Graph(Bloom* bloom, JChecker* jcheck);
};

#endif
