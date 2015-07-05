#ifndef GRAPH
#define GRAPH

#include <unordered_map>
#include <set>
#include <string>
#include "../utils/Kmer.h"
#include "../utils/Junction.h"
#include "../utils/Bloom.h"
#include "../utils/JunctionMap.h"
#include "Node.h"
using std::unordered_map;
using std::set;
using std::string;

class Graph
{
private:

    unordered_map<kmer_type, Node> nodeMap;
    set<kmer_type> cFPs;
    set<kmer_type> sinks;
    Bloom* bloom;
public: 
    void buildGraph(JunctionMap* juncMap);
    void linkAllNodes();
    void printGraph(string fileName);
    Graph(Bloom* bloom);
};

#endif
