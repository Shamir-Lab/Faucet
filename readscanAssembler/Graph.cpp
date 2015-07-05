#include "Graph.h"

using std::string;

Graph::Graph(Bloom* bloo1){
    bloom = bloo1;
    nodeMap = {};
    cFPs = {};
}

void Graph::buildGraph(JunctionMap* juncMap){
    kmer_type kmer;
    Junction junction;
    Node node;
    for(auto it = juncMap->junctionMap.begin(); it != juncMap->junctionMap.end(); it++){
        kmer = it->first;
        junction = it->second;
        if(junction.numPathsOut()>1){
            node = Node(junction);
            nodeMap.insert(std::pair<kmer_type, Node>(kmer, node));
        }
        else if (junction.numPathsOut() == 1){
            for(int i = 0; i < 4; i++){
                if(junction.cov[i] == 0){
                    kmer_type nextKmer = next_kmer(kmer, i, FORWARD);
                    if(bloom->oldContains(nextKmer)){
                        cFPs.insert(nextKmer);
                    }
                }
            }
        }
    }
}

void Graph::linkAllNodes(){

}

void Graph::printGraph(string fileName){

}