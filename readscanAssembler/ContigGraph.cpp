#include "ContigGraph.h"

ContigNode * ContigGraph::getContigNode(kmer_type kmer){

}
    
bool ContigGraph::isContigNode(kmer_type kmer){

}

bool ContigGraph::isErrorContig(Contig* contig){
   // printf("Testing for error contig\n");
    if(!contig->juncDistances){
        //printf("No junc distances\n");
    }
    //printf("Contig length: %d\n", contig->seq.length());
    //printf("Coverage sum: %d\n", contig->coverageSum);
    //printf("Distances size: %d\n", contig->juncDistances->size());
    //std::cout << "Average coverage: " << contig->getAvgCoverage() << "\n";
    bool answer = contig->getAvgCoverage() < 3;
   // printf("Got result of test. Returning.\n");
    return answer;
}

void ContigGraph::deleteContig(Contig* contig){
    if(contig->node1_p){
        cutPath(contig->node1_p, contig->ind1);
    }
    if(contig->node2_p){
        cutPath(contig->node2_p, contig->ind2);
    }
    delete(contig);
}

void ContigGraph::deleteErrorContigs(){
    printf("Deleting error contigs.\n");
    int numDeleted = 0;

    printf("Deleting node mapped contigs.\n");
    //looks through all contigs adjacent to nodes
    for(auto it = nodeMap.begin(); it != nodeMap.end(); it++){
        ContigNode* node = &it->second;
        //printf("Got node.\n");
        for(int i = 0; i < 5; i++){
            if(node->contigs[i]){
                 //printf("Checking contig %d.\n", i);
                Contig* contig = node->contigs[i];
                //printf("Got contig. \n");
                if(isErrorContig(contig)){
                   // printf("Is error contig.\n");
                    numDeleted++;
                    deleteContig(contig);
                   // printf("Deleted contig.\n");
                }
                else{
                   // printf("No error contig\n");
                }
            }
        }
    }

    printf("Deleting isolated contigs.\n");
    //prints isolated contigs
    for(auto it = isolated_contigs.begin(); it != isolated_contigs.end();){
        Contig* contig = &*it;
        if(isErrorContig(contig)){
            numDeleted++;
            it++;
            isolated_contigs.erase(it);
        }
        else{
            it++;
        }
    }

    printf("Done deleting %d error contigs.\n", numDeleted);
}   

int ContigGraph::collapseDummyNodes(){
   printf("Collapsing dummy nodes.\n");
    int numCollapsed = 0;

    //looks through all contigs adjacent to nodes
    for(auto it = nodeMap.begin(); it != nodeMap.end(); ){
        ContigNode* node = &it->second;
        if(node->numPathsOut() == 1){
            numCollapsed++, it++;
            collapseNode(node);
        }
        else it++;
    }

    printf("Done collapsing %d nodes.\n", numCollapsed);
    return numCollapsed;
}

void ContigGraph::cutPath(ContigNode* node, int index){
    if(index == 4){ //if a node has no back path, it turns into sinks instead
        for(int i = 0; i < 4; i++){ //so find all contigs that it points to in the forward direciton
            if(node->contigs[i]){
                Contig* contig = node->contigs[i];
                contig->setSide(contig->getSide(node, i), nullptr); //set them to point to null instead of the node
                if(contig->isIsolated()){ //and add them to isolated_contigs if they're now isolated
                    isolated_contigs.push_back(*contig);
                }
            }
        }
        nodeMap.erase(node->getKmer());// finally erase the node- since it's backpath wasn't touched yet getKmer shouldn't fail
    }
    else{ // if its not the backpath, simply break the path and move on
        node->breakPath(index);
    }
}

void ContigGraph::collapseNode(ContigNode * node){
    Contig* backContig = node->contigs[4];
    Contig* frontContig;
    for(int i = 0; i < 4; i++){
        if(node->contigs[i]) {
            frontContig = node->contigs[i];
        }
    }
    if(frontContig == backContig){ //circular sequence with a redundant node
        addIsolatedContig(*frontContig);
    }
    else{ //normal case of collapsing a node between two other nodes
        ContigNode* backNode = backContig->otherEndNode(node);
        ContigNode* frontNode = frontContig->otherEndNode(node);


        int backSide = backContig->getSide(node);
        int frontSide = frontContig->getSide(node);

        Contig* newContig = backContig->concatenate(frontContig, backSide, frontSide);
        if(backNode){
               backNode->contigs[newContig->ind1] = newContig;
        }
        if(frontNode){
            frontNode->contigs[newContig->ind2] = newContig;
        }
        if(!backNode && !frontNode){
            addIsolatedContig(*newContig);
            delete newContig; 
        }
    }
    kmer_type toErase = node->getKmer();
    if(!nodeMap.erase(toErase)) printf("ERROR: tried to erase node %s but there was no node.\n", print_kmer(toErase));
    delete(backContig);
    delete(frontContig);
}

void ContigGraph::printGraph(string fileName){
    printf("Printing graph from contig graph.\n");
    ofstream jFile;
    jFile.open(fileName);
    int lineNum = 1;
    //printf("Printing contigs from contig graph of %d nodes.\n", nodeMap.size());

    //prints contigs that are adjacent to nodes
    for(auto it = nodeMap.begin(); it != nodeMap.end(); it++){
        ContigNode* node = &it->second;
        for(int i = 0; i < 5; i++){
            if(node->contigs[i]){
                Contig* contig = node->contigs[i];
                if(!contig->node1_p || !contig->node2_p){ //if one side is a sink, always print
                    jFile << contig->getStringRep();
                }
                else if(contig->node1_p == contig->node2_p){ //if the contig attaches to the same node twice, print when you see lower index
                    if(i == contig->getMinIndex()){
                        jFile << contig->getStringRep();;
                    }
              }
                else if(contig->getSide(node) == 1){ //If it attaches to two distinct nodes, print when you're on side 1
                        jFile << contig->getStringRep();
                }
            }
        }
    }

    //prints isolated contigs
    for(auto it = isolated_contigs.begin(); it != isolated_contigs.end(); it++){
        Contig contig = *it;
        //printf("Printing isolated contig.\n");
         jFile << contig.getStringRep();
    }

    //printf("Done printing contigs from contig graph.\n");
    jFile.close();
    printf("Done printing graph from contig graph.\n");
}

void ContigGraph::addIsolatedContig(Contig contig){
    isolated_contigs.push_back(contig);
}

void ContigGraph::printContigs(string fileName){
    printf("Printing contigs from contig graph.\n");
    ofstream jFile;
    jFile.open(fileName);
    int lineNum = 1;
    //printf("Printing contigs from contig graph of %d nodes.\n", nodeMap.size());

    //prints contigs that are adjacent to nodes
    for(auto it = nodeMap.begin(); it != nodeMap.end(); it++){
        ContigNode* node = &it->second;
        for(int i = 0; i < 5; i++){
            if(node->contigs[i]){
                Contig* contig = node->contigs[i];
                if(contig->getSide(node) == 1){
                    //printf("Printing from node at index %d\n", i);
                    // std::cout << "Contig " << lineNum << ":\n";
                    // std::cout << contig->seq << "\n";
                    jFile << canon_contig(contig->seq) << "\n";
                }
            }
        }
    }

    //prints isolated contigs
    for(auto it = isolated_contigs.begin(); it != isolated_contigs.end(); it++){
        Contig contig = *it;
        //printf("Printing isolated contig.\n");
        jFile << canon_contig(contig.seq) << "\n";
    }

    //printf("Done printing contigs from contig graph.\n");
    jFile.close();
    printf("Done printing contigs from contig graph.\n");
}

ContigGraph::ContigGraph(){
    nodeMap = {};
    isolated_contigs = {};
}

//Creates a contig node and returns it or returns the already extant one if it's already extant
//Uses coverage info from junction to create the node
ContigNode * ContigGraph::createContigNode(kmer_type kmer, Junction junc){
    return &(nodeMap.insert(std::pair<kmer_type, ContigNode>(kmer, ContigNode(junc))).first->second);
}