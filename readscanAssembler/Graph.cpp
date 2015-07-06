#include "Graph.h"
#include <time.h>

using std::string;

Graph::Graph(Bloom* bloo1){
    bloom = bloo1;
    nodeMap = {};
    cFPs = {};
    sinks = {};
}

//IMPLEMENTED
//NOT TESTED
//Assumes cFPs and sinks were already handled, only complex junctions should remain
//Gets junctions from juncMap and builds Nodes for graph. Kills junctions as it goes.
void Graph::getNodesFromJunctions(JunctionMap* juncMap){
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

//IMPLEMENTED
//NOT TESTED
//Simply replaces all of the junctions of juncMap with nodes in the graph.  After this the juncMap
//Should be empty and cFPs, sinks, and nodes should exist everywhere they need to.
void Graph::buildGraph(JunctionMap* juncMap){
    //get sinks and cFPs
    time_t start, stop;

    time(&start);
    sinks = juncMap->getSinks();
    time(&stop);

    printf("Sinks found: %d\n", sinks->size());
    printf("Time to find sinks: %f\n", difftime(stop,start));

    time(&start);
    cFPs = juncMap->getCFPs();
    time(&stop);

    printf("cFPs found: %d\n", cFPs->size());
    printf("Time to find cFPs: %f\n", difftime(stop,start));

    time(&start);
    getNodesFromJunctions(juncMap);
    time(&stop);
    printf("Time to get nodes from junctions: %f\n", difftime(stop,start));
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
    kmer_type kmer;
    Node node;
    ofstream jFile;
    jFile.open(fileName);
    printf("Writing graph to file.\n");
    printf("Writing graph nodes to file.\n");
    jFile << "Graph nodes: \n";
    for(auto it = nodeMap.begin(); it != nodeMap.end(); it++){
        kmer = it->first;
        node = it->second;
        jFile << print_kmer(kmer) << " ";
        node.writeToFile(&jFile);
        jFile << "\n";
    }
    printf("Done writing graph nodes to file \n");

    printf("Writing sinks to file.\n");
    jFile <<"Sinks: \n";
    for(auto it = sinks->begin(); it != sinks->end(); it++){
        kmer = *it;
        jFile << print_kmer(kmer) << "\n";
    }
    printf("Done writing sinks to file.\n");

    printf("Writing cFPs to file.\n");
    jFile <<"cFPs: \n";
    for(auto it = cFPs->begin(); it != cFPs->end(); it++){
        kmer = *it;
        jFile << print_kmer(kmer) << "\n";
    }
    printf("Done writing cFPs to file.\n");
    printf("Done writing graph to file.\n");
    jFile.close();
}