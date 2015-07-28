#include "Graph.h"
#include <time.h>
#include <fstream>
#include <iostream> 
#include <sstream>
#include <algorithm>

using std::ofstream;
using std::string;
using std::stringbuf;
using std::min;

Graph::Graph(Bloom* bloo1, JChecker* jcheck){
    bloom = bloo1;
    jchecker = jcheck;
    nodeMap = {};
    contigNodeMap = {};
    realExtensions = {};
    sinks = {};
}

//returns the index from startNode that leads to the destination kmer
int Graph::getPathIndex(Node startNode, kmer_type destinationKmer){
    int index = -1;
    for(int i = 0; i < 5; i++){
        if(startNode.nextJunc[i] == destinationKmer){
            return i;
        }
    }
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

//Gets the valid extension of the given kmer based on the bloom filter and cFPs.  JChecks! so this cuts off tips
//Assume its not a junction
//Returns -1 if there is no valid extension
//Assumes full cFP set and uses it as reference
int Graph::getValidJExtension(DoubleKmer kmer, int dist, int max){
    //printf("Getting valid extension of %s\n", print_kmer(kmer.kmer));
    kmer_type nextKmer;
    int answer = -1;
    for(int i = 0; i < 4; i++){
        nextKmer = kmer.getExtension(i, FORWARD);
        if(bloom->oldContains(get_canon(nextKmer))){
            if(jchecker->jcheck(nextKmer)){
               
                if (answer == -1 ){ //if this is the first correct option found
                    answer = i; 
                } 
                else { //second correct answer! check real extensions
                    if(isRealExtension(kmer.kmer, i)){
                        return i;
                    }
                            //REMOVE AFTER TESTING!
                            auto it = realExtensions->find(kmer.kmer);
                            if(it == realExtensions->end()){
                                printf("Base kmer %s at %d/%d is not in the real extension map.\n", print_kmer(kmer.kmer), dist, max);
                                printf("Real extensions: %d, %d\n", answer, i);
                                return -1;
                            }   
                }
            }
        }
    }
    return answer;
}

//Assumes cFPs and sinks were already handled, only complex junctions should remain.  
//Gets junctions from juncMap and builds equivalent Nodes for graph. Kills junctions as it goes.
void Graph::getNodesFromJunctions(JunctionMap* juncMap){
    kmer_type kmer;
    Junction junction;
    Node node;
    for(auto it = juncMap->junctionMap.begin(); it != juncMap->junctionMap.end(); ){
        kmer = it->first;
        junction = it->second;
        it++; //must update iterator before killing the kmer
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

//Simply replaces all of the junctions of juncMap with nodes in the graph.  After this the juncMap
//Should be empty and cFPs, sinks, and nodes should exist everywhere they need to.
void Graph::buildNodeGraph(JunctionMap* juncMap){
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

//puts in distances and junction IDs to link the two nodes
void Graph::directLinkNodes(kmer_type kmer1, int index1, kmer_type kmer2, int index2, int distance){
    Node* node1 = getNode(kmer1);
    Node* node2 = getNode(kmer2); 
    node1->update(index1, distance, kmer2, index2);
    node2->update(index2, distance, kmer1, index1);
}

//puts in distances and junction IDs to link the node to the sink
void Graph::linkNodeToSink(kmer_type nodeKmer, int index, kmer_type sinkKmer, int distance){
    Node* node = getNode(nodeKmer);
    node->update(index, distance, sinkKmer,-1);
}

void Graph::linkNeighbor(kmer_type startKmer, int index, BfSearchResult result){
    if(result.isNode){
        directLinkNodes(startKmer, index, result.kmer, result.index, result.distance);
    }
    else{
        linkNodeToSink(startKmer, index, result.kmer, result.distance);
    }
}

//Finds the neighbor of the given node on the given extension using the graph.
GraphSearchResult Graph::findNeighborGraph(Node node, kmer_type startKmer, int index){
    kmer_type kmer = node.nextJunc[index];
    bool isNode = nodeMap.find(kmer) != nodeMap.end();
    int otherIndex = -1;
    if(isNode){
        otherIndex = node.backIndex[index];
    }
    int distance = node.dist[index];
    return GraphSearchResult { kmer, isNode, otherIndex, distance };
}

//Finds the neighbor of the given node on the given extension.  
//Returns a result with kmer -1 if there is an error
BfSearchResult Graph::findNeighborBf(Node node, kmer_type startKmer, int index){
    DoubleKmer doubleKmer(startKmer);
    int dist = 1;
    int lastNuc;//stores the last nuc so we know which extension we came from
    string contig = "";

    //First, get to the first kmer from which we can do a proper bloom scan.  This is different for forwards and backward extensions
    if(index == 4){
        //in this case that's the reverse kmer 
        doubleKmer.reverse(); 
        for(int i = 0; i < sizeKmer; i++){
            contig += getNucChar(code2nucleotide(doubleKmer.kmer, i));
        }
        if(isNode(doubleKmer.kmer)){
            return BfSearchResult { doubleKmer.kmer, true, 4, 1, contig };
        }
        if(isSink(doubleKmer.kmer)){
            return BfSearchResult { doubleKmer.kmer, false, -1, 1, contig } ;
        }
        dist = 2;//set it up for second position scan below
    }
    else{
        //in this case thats the next forward kmer- but since we're at a junction we must get there manually using the given index, no bloom scan possible 
        lastNuc =first_nucleotide(doubleKmer.revcompKmer); 
        doubleKmer.forward(index);
        for(int i = 0; i < sizeKmer; i++){
            contig += getNucChar(code2nucleotide(doubleKmer.kmer, i));
        }
        if(isNode(doubleKmer.revcompKmer)){
            return BfSearchResult { doubleKmer.revcompKmer, true, lastNuc, 1, contig };
        }
        dist = 2;
        if(isNode(doubleKmer.kmer)){
            return BfSearchResult { doubleKmer.kmer, true, 4, 2, contig };
        }
        if(isSink(doubleKmer.kmer)){
            return BfSearchResult { doubleKmer.kmer, false, -1, 2, contig };
        }
        dist = 3;//set it up for third position scan below
    }

    //Now, bloom scan forward until there's no chance of finding a junction that indicates an overlapping kmer 
    while(true){ 
        //move forward if possible
        int validExtension = getValidJExtension(doubleKmer, dist,getNode(startKmer)->dist[index] );
        if(validExtension == -1){
            //Note: tested, never at a sink
            printf("ERROR: No valid extension of %s while searching from ",  print_kmer(doubleKmer.kmer));
            printf("%s on index %d at dist %d/%d\n",  print_kmer(startKmer), index, dist, getNode(startKmer)->dist[index]);
            contig += "SINKERROR";
            return BfSearchResult { -1, false, -1, -1, contig };
            //should not happen since we checked for sinks already
        }
        lastNuc = first_nucleotide(doubleKmer.revcompKmer);
        doubleKmer.forward(validExtension); 

        //handle backward junction case
        if(isNode(doubleKmer.revcompKmer)){
            return BfSearchResult { doubleKmer.revcompKmer, true, lastNuc, dist, contig };
        }
        dist++;

        //this happens in the middle of rev and forward directions since we don't want to include 
        //this nuc in the contig if the other junction is facing backward
        contig += getNucChar(validExtension);

        //handle forward junction case
        if(isNode(doubleKmer.kmer)){
            return BfSearchResult { doubleKmer.kmer, true, 4, dist, contig };
        }
        if(isSink(doubleKmer.kmer)){
            return BfSearchResult { doubleKmer.kmer, false, -1, dist, contig };
        }
        dist++;
    }
}

//private utility for doing a BF traversal and either printing contigs or linking nodes or both
void Graph::traverseContigs(bool linkNodes, bool printContigs){
    kmer_type kmer;
    Node node;
    BfSearchResult result;
    ofstream cFile(contigFile);
    for(auto it = nodeMap.begin(); it != nodeMap.end(); it++){
        kmer = it->first;
        node = it->second;
        for(int i = 0; i < 5; i++){
            if(node.cov[i]  > 0 || i == 4){ //if there is coverage or its the backwards direction
                result = findNeighborBf(node, kmer, i);
                if(result.kmer == -1){ //if the search function returned an error, print the error
                    std::cout << "Error occured while searching from " << print_kmer(kmer) << "\n";
                    std::cout << "Search was on index " << i << "\n";
                }
                if(linkNodes && node.nextJunc[i] == -1){ //if we're supposed to link nodes and its not already linked
                    linkNeighbor(kmer, i, result);
                }
                if(printContigs){ //if we're supposed to print contigs, print them
                    if(result.contig <= revcomp_string(result.contig)){
                        cFile << result.contig << "\n";
                    }
                }
            }
        }
    }
    cFile.close();
}

ContigNode * getContigOppositeEnd(ContigNode cnode1, Contig contig){

}


// given Node graph (having all cFPs, sinks, k-mer extensions out of nodes), changes
// to nodes as ends of explicit contigs representation 
void Graph::buildContigGraph(){
    kmer_type kmer;
    Node node;
    Node * far_node;
    BfSearchResult result;
    ContigNode near_end;
    ContigNode far_end;
    Contig * contig;
    string cstr;

    // iterate through original node map
    for(auto it = nodeMap.begin(); it != nodeMap.end(); it++){
        kmer = it->first;
        node = it->second;

        // std::cout << "k-mer: " << print_kmer(kmer) << "\n";

        // if kmer not in contigNodeMap, add cnode:
        if(contigNodeMap.find(kmer) == contigNodeMap.end()){
            near_end = ContigNode(node);
            contigNodeMap.insert(std::pair<kmer_type, ContigNode>(kmer, near_end));
            // std::cout << "inserted k-mer: " << print_kmer(kmer) << "\n";

        }
        for(int i = 0; i < 5; i++){
            //if there is coverage or its the backwards direction, and the contig hasn't been captured yet
            if((node.cov[i]  > 0 || i == 4) && near_end.contigs[i] == NULL){
                // std::cout << "entered if block"<< "\n";

                result = findNeighborBf(node, kmer, i);
                if(result.kmer == -1){ //if the search function returned an error, print the error
                    std::cout << "Error occured while searching from " << print_kmer(kmer) << "\n";
                    std::cout << "Search was on index " << i << "\n";
                }
                cstr = min(result.contig, revcomp_string(result.contig));
                // std::cout << "cstr: "<< cstr << "\n";

                if(contigNodeMap.find(result.kmer) == contigNodeMap.end() && result.isNode){
                    far_node = getNode(result.kmer);
                    // std::cout << "called getNode, i is "<< i << " kmer is "<< print_kmer(result.kmer) << "\n";
                    if (result.kmer==-1){
                        std::cout << "no result.kmer"<< "\n";                    
                    }
                    if (far_node!=NULL){
                        // std::cout << "far_node not null"<< "\n";
                        for(int j = 0; j < 5; j++){
                            // std::cout << far_node->dist[j] << " ";
                        }
                        // std::cout << "\n";
                    }
                    far_end = ContigNode(*far_node);
                    // std::cout << "called ContigNode"<< "\n";                    
                    contigNodeMap.insert(std::pair<kmer_type, ContigNode>(kmer, far_end));
                    // std::cout << "called insert"<< "\n";
                    contig = new Contig(&near_end, i, &far_end, result.index, cstr);
                    // std::cout << "contig seq: "<< contig.seq << "\n";

                    if(far_end.contigs[result.index]==NULL){
                        // far_end->contigs[result.index] = &contig;
                        far_end.update(result.index, contig);
                    }
                    // near_end->contigs[i] = &contig;
                    near_end.update(i, contig);

                }
                
                // std::cout << "contigs[i]->seq: "<< near_end.contigs[i]->seq << "\n";

            }
        }
    }
    // iterate through contig node map to verify it has been loaded
    for(auto it = contigNodeMap.begin(); it != contigNodeMap.end(); it++){
        kmer = it->first;
        near_end = it->second;
        // print kmer
        std::cout << "contigNode k-mer: " << print_kmer(kmer) << "\n";

        for(int i = 0; i < 5; i++){
            // print coverage, contig
            // std::cout << "coverage: " << near_end->cov[i] << "\n";
            // if(near_end->cov[i]  > 0 || i == 4) {
            if( (int) near_end.cov[i] > 0){
                std::cout << (int) near_end.cov[i] << "\n";
            }
            

        }
    }


}


//Assumes the graph has all nodes, sinks, and cFPs properly initialized. 
//Iterates through the nodes and prints every contig path
void Graph::printContigs(string filename){
    contigFile = filename;
    time_t start,stop;
    time(&start);
    printf("Printing contigs.\n");

    traverseContigs(false, true);
    
    printf("Done printing contigs.\n");
    time(&stop);
    printf("Time: %f\n", difftime(stop, start));
}

//Assumes the graph has all nodes, sinks, and cFPs properly initialized. 
//Iterates through the nodes and links all adjacent nodes and sinks.
void Graph::linkNodes(){
    time_t start,stop;
    time(&start);
    printf("Linking nodes.\n");

    traverseContigs(true, false);

    printf("Done linking nodes.\n");
    time(&stop);
    printf("Time: %f\n", difftime(stop, start));
}

bool Graph::isTip(Node node, int index, int maxTipLength){
    return node.cov[index] > 0 && sinks->find(node.nextJunc[index]) != sinks->end() && node.dist[index] <= maxTipLength;
}

int Graph::getNumTips(Node node, int maxTipLength){
    int numTips = 0;
    for(int i = 0; i < 4; i++){
        if(isTip(node, i, maxTipLength)){
            numTips++;
        }
    }
    return numTips;
}

//Deletes a node that only has one real extension, replacing it with an entry in the realExtension map
//Does not link bridging nodes, that must be done beforehand.
void Graph::deleteNode(kmer_type kmer){ 
    if(nodeMap.find(kmer) == nodeMap.end()){
        printf("ERROR: tried to delete nonexistant node.\n");
    }
    Node node = nodeMap.find(kmer)->second;
    int numPaths = 0;
    int realPath = -1;
    for(int i = 0; i < 4; i++){
        if(node.cov[i] > 0){
            numPaths++;
            realPath = i;
        }
    }
    if(numPaths != 1){
        printf("ERROR: Tried to delete a node that did not have exactly one real path.\n");
    }
    (*realExtensions)[kmer] = realPath;
    nodeMap.erase(kmer);
}

int Graph::cutTips(int maxTipLength){
    time_t start,stop;
    time(&start);
    printf("Cutting tips.\n");

    int numTipsCut = 0, numNodesRemoved = 0;
    kmer_type kmer;
    Node* node;
    for(auto it = nodeMap.begin(); it != nodeMap.end(); ){
        kmer = it->first;
        node = &it->second;
        int numTips = getNumTips(*node, maxTipLength);
        if(numTips > 0 && numTips < node->numPathsOut()){
            for(int i = 0; i < 4; i++){
                if(isTip(*node, i, maxTipLength)){
                    if(sinks->find(node->nextJunc[i]) == sinks->end()){
                        printf("ERROR: trimmed a tip but it didn't lead to a sink. \n");
                        printf("Supposed sink kmer %s\n", print_kmer(node->nextJunc[i]));
                    }
                    sinks->erase(node->nextJunc[i]); 
                    node->deletePath(i);
                    Node newNode = nodeMap.find(kmer)->second;
                    if(newNode.cov[i] != 0 || newNode.dist[i] != 0 || newNode.nextJunc[i] != -1){
                        printf("ERROR: delete path doesn't work.\n");
                    }
                    numTipsCut++;
                }
            }

            if(node->numPathsOut() == 1){
                GraphSearchResult last, next;
                for(int i = 0; i < 4; i++){
                    if(node->cov[i] > 0){
                        next = findNeighborGraph(*node, kmer, i);
                    }
                }
                last = findNeighborGraph(*node, kmer, 4);

                if(last.isNode && next.isNode){ //back kmer is a node!
                    directLinkNodes(last.kmer, last.index, next.kmer, next.index, last.distance + next.distance);
                }
                else if (!last.isNode && next.isNode){
                    linkNodeToSink(next.kmer, next.index, last.kmer, last.distance + next.distance);
                }   
                else if (last.isNode && !next.isNode){
                    linkNodeToSink(last.kmer, last.index, next.kmer, last.distance + next.distance);
                }          

                it++;
                deleteNode(kmer);
                numNodesRemoved++;
                continue;
            }
        }
        it++;
    }
    printf("Done cutting %d tips and removing %d nodes.\n", numTipsCut, numNodesRemoved);
    time(&stop);
    printf("Time: %f\n", difftime(stop, start));
    return numTipsCut;
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
