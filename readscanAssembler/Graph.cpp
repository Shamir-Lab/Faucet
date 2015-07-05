#include "Graph.h"

using std::string;

Graph::Graph(Bloom* bloo1){
    bloom = bloo1;
    nodeMap = {};
    cFPs = {};
    sinks = {};
}


//IMPLEMENTED
//NOT TESTED
//Simply replaces all of the junctions of juncMap with nodes in the graph.  After this the juncMap
//Should be empty and cFPs, sinks, and nodes should exist everywhere they need to.
void Graph::buildGraph(JunctionMap* juncMap){
    //get sinks and cFPs
    sinks = juncMap->getSinks();
    cFPs = juncMap->getCFPs();

    //After getCFPs juncMap should only contains complex junctions
    kmer_type kmer;
    Junction junction;
    Node node;
    for(auto it = juncMap->junctionMap.begin(); it != juncMap->junctionMap.end(); it++){
        kmer = it->first;
        junction = it->second;
        if(junction.numPathsOut()>1){
            node = Node(junction);
            nodeMap.insert(std::pair<kmer_type, Node>(kmer, node));
            juncMap->killJunction(kmer);
        }
        else{
            fprintf(stderr, "ERROR: found junction with one or no paths out.\n");
        }
    }
}


//NOT IMPLEMENTED
//NOT TESTED
//Assumes the graph has all nodes, sinks, and cFPs properly initialized. 
//Iterates through the nodes and links all adjacent nodes and sinks.
void Graph::linkAllNodes(){
    kmer_type kmer;
    Node node;
    for(auto it = nodeMap.begin(); it != nodeMap.end(); it++){
        kmer = it->first;
        node = it->second;
        for(int i = 0; i < 5; i++){
            if(node.cov[i]  > 0 && node.nextJunc[i] == -1){ //if there is coverage but not yet a link
                //search along the path i out of this node until a junction or sink is hit.
                //link the node to the node or sink at the end of the path
            }
        }
    }
}

void Graph::printGraph(string fileName){

}