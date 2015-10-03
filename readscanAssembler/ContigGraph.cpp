#include "ContigGraph.h"

//returns true if there are no two identical non-null nodes in the list
bool allDistinct(std::vector<ContigNode*> nodes){
    for(int i = 0; i < nodes.size(); i++){
        for(int j = i+1; j < nodes.size(); j++){
            if(nodes[i] == nodes[j] && nodes[i]){
                return false;
            }
        }
    }
    return true;
}

void ContigGraph::switchToNodeVector(){
    for(auto it = nodeMap.begin(); it != nodeMap.end(); ){
        kmer_type kmer = it->first;
        ContigNode* originalNode = &it->second;
        nodeVector.push_back(*originalNode);
        ContigNode* newNode = &nodeVector.back();
        for(int i = 0; i < 5; i++){
            if(newNode->contigs[i]){
                if(newNode->contigs[i]->node1_p == originalNode){
                    newNode->contigs[i]->node1_p = newNode;
                } 
                if(newNode->contigs[i]->node2_p == originalNode){
                    newNode->contigs[i]->node2_p = newNode;
                }
            }
        }
        it++;
        nodeMap.erase(kmer);
    }
}

bool ContigGraph::cleanGraph(){
    bool result = false;
    if(deleteErrorContigs() > 0){
        result = true;
    }
    if(destroyDegenerateNodes() > 0){
        result = true;
    }
    if(collapseDummyNodes() > 0){
        result = true;
    }
    return result;
}

void ContigGraph::checkGraph(){
    printf("Checking node mapped contigs.\n");
    for(auto it = nodeMap.begin(); it != nodeMap.end(); it++){
        kmer_type kmer = it->first;
        ContigNode* node = &it->second;
        if(kmer != node->getKmer()){
            printf("GRAPHERROR: node kmer didn't match nodemap kmer.\n");
        }
        if(!node->contigs[4]){
            printf("GRAPHERROR: node lacks backcontig.\n");
        }
        for(int i = 0; i < 5; i++){
            if(node->contigs[i]){
                Contig* contig = node->contigs[i];
                int side = contig->getSide(node, i);
                if(side == 1){
                    if(contig->ind1 != i){
                        printf("GRAPHERROR: contig has wrong index.\n");
                    }
                    if(contig->node1_p != node){
                        printf("GRAPHERROR: contig points to wrong node.\n");
                    }
                }
                if(side == 2){
                    if(contig->ind2 != i){
                        printf("GRAPHERROR: contig has wrong index.\n");
                    }
                    if(contig->node2_p != node){
                        printf("GRAPHERROR: contig points to wrong node.\n");
                    }
                }
            }
        }
    }

    printf("Checking isolated contigs.\n");
    //prints isolated contigs
    for(auto it = isolated_contigs.begin(); it != isolated_contigs.end();){
        Contig* contig = &*it;
        if(contig->node1_p || contig->node2_p){
            printf("GRAPHERROR: isolated contig has pointer to at least one node.\n");
        }  
    }

}

bool ContigGraph::isErrorContig(Contig* contig){
    if(contig->getAvgCoverage()<3){
        return true;
    }
    if (20*contig->getAvgCoverage() < contig->getMinAdjacentCoverage()) { //more deletion for chimeras
        return true;
    }   
    return false;
}

void ContigGraph::deleteContig(Contig* contig){
    if(contig->node1_p){
        cutPath(contig->node1_p, contig->ind1);
    }
    if(contig->node2_p){
        cutPath(contig->node2_p, contig->ind2);
    }
    if(contig){ //shouldn't be necessary but we'll see
        delete contig;
    }
}

int ContigGraph::deleteErrorContigs(){
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
            }
        }
    }

    printf("Done deleting node mapped contigs.\n");
    printf("Deleting isolated contigs.\n");
    //prints isolated contigs
    for(auto it = isolated_contigs.begin(); it != isolated_contigs.end();){
        Contig* contig = &*it;
        if(isErrorContig(contig)){
            numDeleted++, it++;
            isolated_contigs.erase(it);
        }
        else it++;
    }

    printf("Done deleting %d error contigs.\n", numDeleted);
    return numDeleted;
}   

int ContigGraph::destroyDegenerateNodes(){
    printf("Destroying degenerate nodes.\n");
    int numDegen = 0;

    //looks through all contigs adjacent to nodes
    for(auto it = nodeMap.begin(); it != nodeMap.end(); ){
        ContigNode* node = &it->second;
        kmer_type kmer = it->first;
        it++;
        if(node->numPathsOut() == 0){
            if(node->contigs[4]){
                cutPath(node, 4);
            }
            numDegen++;
            nodeMap.erase(kmer);
        }
        else if(!node->contigs[4]){
            for(int i = 0; i < 4; i++){
                if(node->contigs[i]){
                    cutPath(node,i);
                }
            }
            numDegen++;
            nodeMap.erase(kmer);
        }
    }

    printf("Done destroying %d nodes.\n", numDegen);
    return numDegen;
}

int ContigGraph::disentangle(Bloom* pair_filter){
    int disentangled = 0;

    //looks through all contigs adjacent to nodes
    for(auto it = nodeMap.begin(); it != nodeMap.end(); ){
        //printf("-1\n");
        ContigNode* node = &it->second;
        //printf("Testing %s\n", print_kmer(node->getKmer()));
        kmer_type kmer = it->first;
       // printf("0\n");
        if(node->numPathsOut() == 2){
            //printf("1");
            Contig* contig = node->contigs[4];
            ContigNode* backNode = contig->otherEndNode(node);
          // printf(".2");
            if(backNode && node != backNode && backNode->numPathsOut() == 2 && backNode->indexOf(contig) == 4){
                // printf("Found an outward facing pair: %s,", print_kmer(node->getKmer()));
                // printf("%s\n", print_kmer(backNode->getKmer()));
                int a,b,c,d;
                a = backNode->getIndicesOut()[0], b = backNode->getIndicesOut()[1];
                c = node->getIndicesOut()[0], d = node->getIndicesOut()[1];
                kmer_type A,B,C,D;
                A = backNode->getForwardExtension(a), B = backNode->getForwardExtension(b);
                C = node->getForwardExtension(c), D = node->getForwardExtension(d);

                // printf("A: %s\n", print_kmer(A));
                // printf("B: %s\n", print_kmer(B));
                // printf("C: %s\n", print_kmer(C));
                // printf("D: %s\n", print_kmer(D));

                // printf("Pair filter says:\n");
                // printf("A/C: %d\n", pair_filter->containsPair(A,C));
                // printf("B/D: %d\n", pair_filter->containsPair(B,D));
                // printf("A/D: %d\n", pair_filter->containsPair(A,D));
                // printf("B/C: %d\n", pair_filter->containsPair(B,C));
                if(pair_filter->containsPair(A,C) && pair_filter->containsPair(B,D) && 
                    ! pair_filter->containsPair(A,D) && !pair_filter->containsPair(B,C)){
                    //printf("Found pair in filter in first orientation.\n");
                    
                    if(disentanglePair(contig, backNode, node, a, b, c, d)){
                        //printf("Disentangled pair.\n");
                        it = nodeMap.erase(it);
                        //printf("1\n");
                        if(it != nodeMap.end()){
                            //printf("2\n");
                            if(backNode->getKmer() == it->first){
                                //printf("3\n");
                                it++;
                            }
                        }
                        //printf("4\n");
                        nodeMap.erase(backNode->getKmer());
                        //printf("5\n");
                        disentangled++;
                        //printf("6\n");
                        continue;
                    }
                }
                if(pair_filter->containsPair(A,D) && pair_filter->containsPair(B,C) && 
                    ! pair_filter->containsPair(A,C) && !pair_filter->containsPair(B,D)){
                    //printf("Found pair in filer in second orientation\n");
                    if(disentanglePair(contig, backNode, node, a,b,d,c)){
                        //printf("Disentangled pair.\n");
                        it = nodeMap.erase(it);
                        if(it != nodeMap.end()){
                            if(backNode->getKmer() == it->first){
                                it++;
                            }
                        }
                        nodeMap.erase(backNode->getKmer());
                        disentangled++;
                        continue;
                    } 
                }
            }
       }
       it++;
    }

    printf("Done disentangling %d pairs of nodes.\n", disentangled);
    return disentangled;
}

//a,b are on backNode, c,d are on forwardNode
//a pairs with c, b pairs with d
/*
       a\                       /c              a--------c
         --------contig---------      ---->
       b/                       \d              b--------d
*/
bool ContigGraph::disentanglePair(Contig* contig, ContigNode* backNode, ContigNode* forwardNode, 
    int a, int b, int c, int d){
    Contig* contigA = backNode->contigs[a];
    Contig* contigB = backNode->contigs[b];
    Contig* contigC = forwardNode->contigs[c];
    Contig* contigD = forwardNode->contigs[d];

    ContigNode* nodeA = contigA->otherEndNode(backNode);
    ContigNode* nodeB = contigB->otherEndNode(backNode);
    ContigNode* nodeC = contigC->otherEndNode(forwardNode);
    ContigNode* nodeD = contigD->otherEndNode(forwardNode);

    if(!allDistinct({backNode, forwardNode, nodeA, nodeB, nodeC, nodeD})){
        return false;
    }
    Contig* contigAC = contigA->concatenate(contig, contigA->getSide(backNode), contig->getSide(backNode));
    contigAC = contigAC->concatenate(contigC, contigAC->getSide(forwardNode), contigC->getSide(forwardNode));
    if(nodeA){
        nodeA->replaceContig(contigA, contigAC);
    }
    if(nodeC){
        nodeC->replaceContig(contigC, contigAC);
    }
    if(!nodeA && !nodeC){
        isolated_contigs.push_back(*contigAC);
    }

    Contig* contigBD = contigB->concatenate(contig, contigB->getSide(backNode), contig->getSide(backNode));
    contigBD = contigBD->concatenate(contigD, contigBD->getSide(forwardNode), contigD->getSide(forwardNode));
    if(nodeB){
        nodeB->replaceContig(contigB, contigBD);
    }
    if(nodeD){
        nodeD->replaceContig(contigD, contigBD);
    }
    if(!nodeB && !nodeD){
        isolated_contigs.push_back(*contigBD);
    }
    return true;
}
    
int ContigGraph::collapseDummyNodes(){
   printf("Collapsing dummy nodes.\n");
    int numCollapsed = 0;

    //looks through all contigs adjacent to nodes
    for(auto it = nodeMap.begin(); it != nodeMap.end(); ){
        ContigNode* node = &it->second;
        kmer_type kmer = it->first;
        it++;
        if(node->numPathsOut() == 1){
            numCollapsed++;
            collapseNode(node);
        }
    }

    printf("Done collapsing %d nodes.\n", numCollapsed);
    return numCollapsed;
}

void ContigGraph::cutPath(ContigNode* node, int index){
    if(!node->contigs[index]){
        printf("ERROR: tried to cut nonexistant path.");
    }
    Contig* contig = node->contigs[index];
    int side = contig->getSide(node, index);
    int otherSide = 3 - side;
    if(contig->node1_p == contig->node2_p && contig->ind1 == contig->ind2){ //to handle hairpins
        contig->setSide(side, nullptr); //set to point to null instead of the node
        contig->setSide(otherSide, nullptr);
    }
    else{
        contig->setSide(side, nullptr);
    }
    node->breakPath(index);
}

void ContigGraph::collapseNode(ContigNode * node){
    Contig* backContig = node->contigs[4];
    Contig* frontContig;
    int fronti = 0;
    for(int i = 0; i < 4; i++){
        if(node->contigs[i]) {
            frontContig = node->contigs[i];
            fronti = i;
        }
    }
    if(frontContig == backContig){ //circular sequence with a redundant node
        addIsolatedContig(*frontContig);
    }
    else{ //normal case of collapsing a node between two other nodes
        ContigNode* backNode = backContig->otherEndNode(node);
        ContigNode* frontNode = frontContig->otherEndNode(node);


        int backSide = backContig->getSide(node, 4);
        int frontSide = frontContig->getSide(node, fronti);

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
                else if(contig->getSide(node,i) == 1){ //If it attaches to two distinct nodes, print when you're on side 1
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
                if(contig->getSide(node,i) == 1){
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
    nodeVector = {};
}

//Creates a contig node and returns it or returns the already extant one if it's already extant
//Uses coverage info from junction to create the node
ContigNode * ContigGraph::createContigNode(kmer_type kmer, Junction junc){
    return &(nodeMap.insert(std::pair<kmer_type, ContigNode>(kmer, ContigNode(junc))).first->second);
}