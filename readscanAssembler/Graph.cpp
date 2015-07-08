#include "Graph.h"
#include <time.h>
#include <fstream>

using std::ofstream;
using std::string;

Graph::Graph(Bloom* bloo1, JChecker* jcheck){
    bloom = bloo1;
    jchecker = jcheck;
    nodeMap = {};
    realExtensions = {};
    sinks = {};
}

Node * Graph::getNode(kmer_type kmer){
    return &(nodeMap.find(kmer)->second);
}

bool Graph::isSink(kmer_type kmer){
    return sinks->find(kmer) != sinks->end();
}

bool Graph::isNode(kmer_type kmer){
    return nodeMap.find(kmer) != nodeMap.end();
}

bool Graph::isRealExtension(kmer_type kmer, int ext){
    auto it = realExtensions->find(kmer);
    if(it == realExtensions->end()){
        return false;
    }
    return it->second == ext;
}

//IMPLEMENTED
//NOT TESTED
//Note: as tested, there is never a case with ambiguous options.  So that can't be cause of bugs.
//Gets the valid extension of the given kmer based on the bloom filter and cFPs.  Must JCheck! so this cuts off tips
//Assume its not a junction
//Returns -1 if there is no valid extension
//Assumes full cFP set and uses it as reference
int Graph::getValidJExtension(DoubleKmer kmer, int dist, int max){
    //printf("Getting valid extension of %s\n", print_kmer(kmer.kmer));
    kmer_type nextKmer;
    int answer = -1;
    for(int i = 0; i < 4; i++){
        nextKmer = kmer.getExtension(i, FORWARD);
        //printf("Testing extension %s\n", print_kmer(nextKmer));
        if(bloom->oldContains(get_canon(nextKmer))){
            if(jchecker->jcheck(nextKmer)){
                if(answer != -1){

                            //REMOVE AFTER TESTING!
                            auto it = realExtensions->find(kmer.kmer);
                            if(it == realExtensions->end()){
                                printf("Base kmer %s at %d/%d is not in the real extension map.\n", print_kmer(kmer.kmer), dist, max);
                                printf("Real extensions: %d, %d\n", answer, i);
                            }

                    //second correct answer! check real extensions
                    if(isRealExtension(kmer.kmer, i)){
                        return i;
                    }
                            
                }
                else{
                    answer = i; //set answer to first option found
                }
            }
        }
    }
    return answer;
}

//IMPLEMENTED
//TESTED
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
//TESTED
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
    realExtensions = juncMap->getRealExtensions();
    time(&stop);

    printf("Real extensions found: %d\n", realExtensions->size());
    printf("Time to find real extensions: %f\n", difftime(stop,start));

    time(&start);
    getNodesFromJunctions(juncMap);
    time(&stop);
    printf("Time to get nodes from junctions: %f\n", difftime(stop,start));
}

//IMPLEMENTED
//NOT TESTED
//puts in distances and junction IDs to link the two nodes
void Graph::directLinkNodes(kmer_type kmer1, int index1, kmer_type kmer2, int index2, int distance){
    Node* node1 = getNode(kmer1);
    Node* node2 = getNode(kmer2); 
    node1->update(index1, distance, kmer2);
    node2->update(index2, distance, kmer1);
}

//IMPLEMENTED
//NOT TESTED
//puts in distances and junction IDs to link the node to the sink
void Graph::linkNodeToSink(kmer_type nodeKmer, int index, kmer_type sinkKmer, int distance){
    Node* node = getNode(nodeKmer);
    node->update(index, distance, sinkKmer);
}

//IMPLEMENTED
//NOT TESTED
//Note: tested that no backwards kmer is ever a sink.  So that's good.
//finds the neighbor of the given node on the given extension, if it exists.
//If it does exist, the two are linked with references to each other and the appropriate distance.
//Also, the resulting contig is printed to the cFile, with NO PRINT LINE
void Graph::findAndLinkNeighbor(Node node, kmer_type startKmer , int index, ofstream* cFile){
    DoubleKmer doubleKmer(startKmer);
    int dist = 0;
    int lastNuc;//stores the last nuc so we know which extension we came from

    //get to the first kmer from which we need to bloom scan.  This is different for forwards and backward extensions
    if(index == 4){
        //in this case that's the reverse kmer 
        doubleKmer.reverse(); 
        for(int i = 0; i < sizeKmer; i++){
            *cFile << getNucChar(code2nucleotide(doubleKmer.kmer, i));
        }
        if(isNode(doubleKmer.kmer)){
            directLinkNodes(startKmer, index, doubleKmer.kmer, 4, 1);
            return;
        }
        if(isSink(doubleKmer.kmer)){
            linkNodeToSink(startKmer, index, doubleKmer.kmer, 1);
            return;
        }
        dist = 2;//set it up for second position scan below
    }
    else{
        //in this case thats the forward version of the kmer extension specified by the index
        lastNuc =first_nucleotide(doubleKmer.revcompKmer); 
        doubleKmer.forward(index);
        for(int i = 0; i < sizeKmer; i++){
            *cFile << getNucChar(code2nucleotide(doubleKmer.kmer, i));
        }
        if(isNode(doubleKmer.revcompKmer)){
            directLinkNodes(startKmer, index, doubleKmer.revcompKmer, lastNuc, 1);
            return;
        }
        if(isNode(doubleKmer.kmer)){
            directLinkNodes(startKmer, index, doubleKmer.kmer, 4, 2);
            return;
        }
        if(isSink(doubleKmer.kmer)){
            linkNodeToSink(startKmer, index, doubleKmer.kmer, 2);
            return;
        }
        dist = 3;//set it up for third position scan below
    }

    //Scan forward until there's no chance of finding a junction that indicates an overlapping kmer 
    while(true){ 
        //move forward if possible
        int validExtension = getValidJExtension(doubleKmer, dist,getNode(startKmer)->dist[index] );
        if(validExtension == -1){
            //Note: tested, never at a sink
            printf("ERROR: No valid extension of %s while searching from ",  print_kmer(doubleKmer.kmer));
            printf("%s on index %d at dist %d/%d\n",  print_kmer(startKmer), index, dist, getNode(startKmer)->dist[index]);
            *cFile << " SINKERROR";
            return;
            //should not happen since we checked for sinks already
        }
        lastNuc = first_nucleotide(doubleKmer.revcompKmer);
        doubleKmer.forward(validExtension); 

        //handle backward junction case
        if(isNode(doubleKmer.revcompKmer)){
            directLinkNodes(startKmer, index, doubleKmer.revcompKmer, lastNuc, dist);
            return;
        }
        dist++;

        //this happens in the middle of rev and forward directions since we don't want to include 
        //this nuc in the contig if the other junction is facing backward
        *cFile << getNucChar(validExtension);

        //handle forward junction case
        if(isNode(doubleKmer.kmer)){
            directLinkNodes(startKmer, index, doubleKmer.kmer, 4, dist);
            return; 
        }
        if(isSink(doubleKmer.kmer)){
            linkNodeToSink(startKmer, index, doubleKmer.kmer, dist);
            return;
        }

        dist++;
    }
}

//IMPLEMENTED
//NOT TESTED
//Assumes the graph has all nodes, sinks, and cFPs properly initialized. 
//Iterates through the nodes and links all adjacent nodes and sinks.
void Graph::linkNodesPrintContigs(string fileName){
    time_t start,stop;
    time(&start);
    printf("Linking nodes and printing contigs.\n");
    kmer_type kmer;
    Node node;
    ofstream cFile(fileName);
    for(auto it = nodeMap.begin(); it != nodeMap.end(); it++){
        kmer = it->first;
        node = it->second;
        for(int i = 0; i < 5; i++){
            if((node.cov[i]  > 0 || i == 4) && node.nextJunc[i] == -1){ //if there is coverage but not yet a link
                findAndLinkNeighbor(node, kmer, i, &cFile);
                cFile << "\n";
            }
        }
    }
    cFile.close();
    printf("Done linking nodes and printing contigs.\n");
    time(&stop);
    printf("Time: %f\n", difftime(stop, start));
}

void Graph::printGraph(string fileName){
    kmer_type kmer;
    Node node;
    int ext;
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

    printf("Writing real extensions to file.\n");
    jFile <<"Real extensions: \n";
    for(auto it = realExtensions->begin(); it != realExtensions->end(); it++){
        kmer = it->first;
        ext = it->second;
        jFile << print_kmer(kmer) << " " << getNucChar(ext) << "\n";
    }
    printf("Done writing real extensions to file.\n");
    printf("Done writing graph to file.\n");
    jFile.close();
}