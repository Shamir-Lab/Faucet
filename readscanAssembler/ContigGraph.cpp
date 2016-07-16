#include "ContigGraph.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>

unordered_map<kmer_type, ContigNode> *  ContigGraph::getNodeMap(){
    return &nodeMap;
}

//returns true if there are no two identical non-null nodes in the list
template<class T>
bool allDistinct(const std::vector<T> & nodes)
{ 
    for(int i = 0; i < nodes.size(); i++){
        for(int j = i+1; j < nodes.size(); j++){
            if(nodes[i] == nodes[j] && nodes[i]){
                return false;
            }
        }
    }
    return true; 
} 


//returns true if all nodes in the list are non-null and identical
bool allSame(std::vector<ContigNode*> nodes){
    for(int i = 0; i < nodes.size(); i++){
        for(int j = i+1; j < nodes.size(); j++){
            if(nodes[i] != nodes[j] || !nodes[i]){
                return false;
            }
        }
    }
    return true;
}

void ContigGraph::setReadLength(int length){
    read_length = length;
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

bool ContigGraph::deleteTipsAndClean(){
    bool result = false;
    if(deleteTipsAndLowCoverageContigs() > 0){
        result = true;
    } 
    // if(popBubblesByCoverageRatio() > 0){
    //     result = true;
    // }  
    // if(destroyDegenerateNodes() > 0){
    //     result = true;
    // }
    // if(collapseDummyNodes() > 0){
    //     result = true;
    // }
    return result;
}

bool ContigGraph::breakPathsAndClean(Bloom* pair_filter, int insertSize){
    bool result = false;

    // if(breakUnsupportedPaths(pair_filter, insertSize) > 0){
    //     result = true;
    // }
    if(collapseBulges(150) > 0){
        result = true;
    }
    // if(destroyDegenerateNodes() > 0){
    //     result = true;
    // }
    // if(collapseDummyNodes() > 0){
    //     result = true;
    // }

    return result;
}

bool ContigGraph::disentangleAndClean(Bloom* pair_filter, int insertSize){
    bool result = false;
    if(disentangle(pair_filter, insertSize, true) > 0){
        result = true;
    }
    // if(disentangle(pair_filter, insertSize, false) > 0){
    //     result = true;
    // }
    // if(collapseBulges(150) > 0){
    //     result = true;
    // }
    if(destroyDegenerateNodes() > 0){
        result = true;
    }
    if(collapseDummyNodes() > 0){
        result = true;
    }

    return result;
}

bool ContigGraph::cleanGraph(Bloom* short_pair_filter, Bloom* long_pair_filter, int insertSize){

    bool result = false;
    deleteTipsAndClean();
    if(breakPathsAndClean(short_pair_filter, read_length)){
        result = true;
    }
    if(disentangleAndClean(short_pair_filter, read_length)){
        result = true;
    }
    if(disentangleAndClean(long_pair_filter, 500)){
        result = true;
    }
    
   

    return result;
}

bool ContigGraph::checkGraph(){
    printf("Checking node mapped contigs.\n");
    for(auto it = nodeMap.begin(); it != nodeMap.end(); it++){
        kmer_type kmer = it->first;
        ContigNode* node = &it->second;

        //since it uses the back contig to generate the kmer, can't run this check without one
        if(node->contigs[4]){
            if(kmer != node->getKmer()){
                printf("GRAPHERROR: node kmer didn't match nodemap kmer.\n");
                return false;
            }
        }

        //Want to run this check during cleaning, so can't look for this -it'll be violated till we destroy degenerate nodes
        // if(!node->contigs[4]){
        //     printf("GRAPHERROR: node lacks backcontig.\n");
        // }
        
        if(!node->checkValidity()){
            return false;
        }
        for(int i = 0; i < 5; i++){
            if(node->contigs[i]){
                if(!node->contigs[i]->checkValidity()){
                    return false;
                }
            }
        }                
    }

    printf("Checking %d isolated contigs.\n", isolated_contigs.size());
    //prints isolated contigs
    for(auto it = isolated_contigs.begin(); it != isolated_contigs.end(); it++){
        Contig* contig = &*it;
        if(contig->node1_p || contig->node2_p){
            printf("GRAPHERROR: isolated contig has pointer to at least one node.\n");
            return false;
        }  
    }

    printf("Done checking graph\n");
    return true;

}

//Returns a list of extensions which are there but have no support
std::vector<int> ContigGraph::getUnsupportedExtensions(ContigNode* node, Bloom* pair_filter, int insertSize){
    if(!node->contigs[4]){
        printf("GRAPHERROR: tried to test support at node with no backcontig.\n");
        return {};
    }
    // int insertSize = 500;
    std::list<JuncResult> backResults = node->getPairCandidates(4, insertSize);
    std::list<JuncResult> forResults;
    std::vector<int> results = {};
    for (int i = 0; i < 4; i++){
        if(node->contigs[i]){
            Contig* contig = node->contigs[i];
            forResults = node->getPairCandidates(i, insertSize);
            if(getScore(backResults, forResults, pair_filter, .01, insertSize)==0){
                results.push_back(i);
            }
        }   
    }
    return results;
}

bool ContigGraph::isTip(ContigNode* node, int index){
    Contig* contig = node->contigs[index];
    if(contig->getSeq().length() < read_length && contig->otherEndNode(node) == nullptr){
        return true;
    }
    return false;
}

bool ContigGraph::isLowCovContig(Contig* contig){
    if(contig->getAvgCoverage() < 3){ //* contig->getSeq().length() < 250){
        return true;
    }
    // if (20*contig->getAvgCoverage() < contig->getMinAdjacentCoverage()) { //more deletion for chimeras
    //     return true;
    // }   
    return false;
}

bool ContigGraph::isLowMassContig(Contig* contig){
    if(contig->getAvgCoverage() * contig->getSeq().length() < 500){
        return true;
    }
    return false;
}

void ContigGraph::deleteContig(Contig* contig){
    //std::cout << contig << "\n";
    ContigNode * node1 = contig->node1_p;
    ContigNode * node2 = contig->node2_p;


    if(node1){
        node1->breakPath(contig->ind1); // remove ptr, cov on node
    }
    if(node2){
        node2->breakPath(contig->ind2);        
    }

    delete contig;
    contig = nullptr;
}


void ContigGraph::cutPath(ContigNode* node, int index){
    //std::cout << "Node: " << node << "\n";
    //std::cout << "Cov: \n" << node->getCoverage(0) << "\n";
    if(!node->contigs[index]){
        printf("ERROR: tried to cut nonexistant path.");
    }
    //printf("Did error check \n");
    Contig* contig = node->contigs[index];
    //printf("Getting side\n");
    int side = contig->getSide(node, index);
    int otherSide = 3 - side;
    //printf("A\n");
    if(contig->node1_p == contig->node2_p && contig->ind1 == contig->ind2){ //to handle hairpins
       // printf("A1\n");
        contig->setSide(side, nullptr); //set to point to null instead of the node
        contig->setSide(otherSide, nullptr);
    }
    else{
        //printf("A2\n");
        contig->setSide(side, nullptr);
    }
    //printf("To break path.\n");
    node->breakPath(index);
}



int ContigGraph::breakUnsupportedPaths(Bloom* pair_filter, int insertSize){
    printf("Breaking unsupported paths contigs. Starting with %d nodes. \n", nodeMap.size());
    int numBroken = 0;
    clock_t t = clock();
    //looks through all contigs adjacent to nodes
    int j = 1;
    for(auto it = nodeMap.begin(); it != nodeMap.end(); it++){
        if(j % 1000 == 0){
            std::cout << j << " nodes processed. Time: "<< clock()-t << "\n";
            t = clock();
        }
        j++;
        ContigNode* node = &it->second;
        //printf("Got node.\n");
        std::vector<int> unsupported = getUnsupportedExtensions(node, pair_filter,insertSize);
        for(auto it = unsupported.begin(); it != unsupported.end(); it++){
                    numBroken++;
                    // printf("Deleting contig %d \n" , numDeleted);
                    // std::cout << "Sequence: " << contig->seq << "\n";
                    // std::cout << "Header: " << contig->getFastGHeader(true) << "\n";
                    cutPath(node,*it);

                    // if(!checkGraph()){
                    //     std::cout << " Graph Check failed. Exiting.\n";
                    //     exit(1);
                    // }
                    // printGraph("lastValidGraph.fastg");
                    // printf("Deleted contig.\n");
                
            
        }
    }

    printf("Done breaking %d unsupported paths.\n", numBroken);
    return numBroken;
}   

int ContigGraph::deleteTipsAndLowCoverageContigs(){
    printf("Deleting tips and low coverage contigs. Starting with %d nodes. \n", nodeMap.size());
    int numDeleted = 0;
    std::set<kmer_type> seenKmers = {};

    //looks through all contigs adjacent to nodes
    int k = 1;
    for(auto it = nodeMap.begin(); it != nodeMap.end(); it++){
        if(k % 10000 == 0){
            std::cout << k << " nodes processed.\n";
        }
        k++;
        ContigNode* node = &it->second;
        kmer_type kmer = it->first;
        if (seenKmers.find(kmer) == seenKmers.end()){
            for(int i = 0; i < 5; i++){
                Contig* contig = node->contigs[i];
                if(contig){
                    if(isTip(node, i)){
                        numDeleted++;
                        deleteContig(contig);

                    }
                    
                    else if(isLowCovContig(contig)){
                        numDeleted++;
                        ContigNode* far_node = contig->otherEndNode(node);
                        deleteContig(contig);
                        kmer_type far_kmer;
                        if (far_node){
                            far_kmer = far_node->getKmer();
                        }

                        if (far_node){
                            if (testAndCutIfDegenerate(far_node)) seenKmers.insert(far_kmer);

                            if (far_node == node){ 
                                // hairpin - break to avoid revisiting collapsed node
                                // treating node (below) will make sure it is collapsed 
                                break;
                            }
                            else if(far_node->numPathsOut() == 1){
                                collapseNode(far_node, far_kmer);         
                                seenKmers.insert(far_kmer);
                       
                            }
                        }
                    }
                } 
            }     
            
            if (testAndCutIfDegenerate(node)) seenKmers.insert(kmer);

            if(node->numPathsOut() == 1 && seenKmers.find(kmer) == seenKmers.end()){
                collapseNode(node, kmer);         
                seenKmers.insert(kmer);
       
            }
            
        }

    }
    for (auto k_it = seenKmers.begin(); k_it != seenKmers.end(); k_it++){
        nodeMap.erase(*k_it);
    }
    printf("Deleted %d tips and low coverage contigs. %d nodes remain\n", numDeleted, nodeMap.size());
    numDeleted = 0;
    printf("Deleting isolated low mass contigs. Starting with %d.\n", isolated_contigs.size());
    //prints isolated contigs
    for(auto it = isolated_contigs.begin(); it != isolated_contigs.end();){
        Contig* contig = &*it;
        if(isLowMassContig(contig)){
            numDeleted++;
            it = isolated_contigs.erase(it);
        }
        else it++;
    }

    printf("Deleted %d isolated low mass contigs. %d nodes remain\n", numDeleted, isolated_contigs.size());
    return numDeleted;
}   

bool ContigGraph::testAndCutIfDegenerate(ContigNode* node){
    if(node->numPathsOut() == 0){
        if(node->contigs[4]){
            cutPath(node, 4);
        }
        return true;
    }
    else if(!node->contigs[4]){
        for(int j = 0; j < 4; j++){
            if(node->contigs[j]){
                cutPath(node,j);
            }
        }
        return true;
    }
    return false;
}


int ContigGraph::destroyDegenerateNodes(){
    printf("Destroying degenerate nodes. Starting with %d nodes.\n", nodeMap.size());
    int numDegen = 0;

    //looks through all contigs adjacent to nodes
    // for(auto it = nodeMap.begin(); it != nodeMap.end(); ){
    it = nodeMap.begin();
    while(it != nodeMap.end()){
        ContigNode* node = &it->second;
        if (testAndCutIfDegenerate(node)) {
            numDegen++;
            it = nodeMap.erase(it);
        }
        else{
            it++;
        }
    }

    printf("Done destroying %d nodes.\n", numDegen);
    return numDegen;
}


double ContigGraph::getTailBound(int numTrials, double p, int result){
    double mean = 1.0*numTrials*p;
    double delta = (1.0*result / mean) - 1;
    double chernoffBound;
    if( delta < 0){
        chernoffBound = 1.0;
    }
    else if(delta < 1){
        chernoffBound = exp(-delta*delta*mean/3);
    }
    else if (delta > 1){
        chernoffBound = exp(-delta*mean/3);
    }
    //std::cout <<" Chernoff: " << chernoffBound << " \n";
    double markovBound = 1.0;
    if (result >  mean){
        markovBound = mean/result;
    }
    //std::cout <<" Markov: " << markovBound << " \n";
    //std::cout << "Min: " << std::min(chernoffBound, markovBound) << "\n";
    return std::min(chernoffBound, markovBound);
}

//returns a score based on how many pairs of kmers from the first and second lists are in the filter,
//relative to the FP rate of the filter
double ContigGraph::getScore(std::list<JuncResult> leftCand, std::list<JuncResult> rightCand, Bloom* pair_filter, double fpRate, int insertSize){
    double score = 0;
    // int insertSize = 500;
    // int readLength = read_length; //100;
    int counter = 0;
    std::unordered_set<JuncPair> seenPairs = {};
    
    //Looks for junction pairs from a single read, within readlength of the junction either way
    for(auto itL = leftCand.begin(); itL != leftCand.end() && itL->distance < insertSize; itL++){
        for(auto itR = rightCand.begin(); itR != rightCand.end()  && itR->distance < insertSize; itR++){
            counter++;
            JuncPair pair = JuncPair(itL->kmer, itR->kmer);
            JuncPair rev_pair = JuncPair(itR->kmer, itL->kmer);
            if(pair_filter->containsPair(pair) && seenPairs.find(pair)==seenPairs.end() 
                && seenPairs.find(rev_pair)==seenPairs.end()){
                //std::cout << "Distance: " << itL->distance + itR->distance << "\n";
                score += 1;
                seenPairs.insert(pair);
            } 
        }
    }

    //Looks for junction pairs one full insert size apart
    //Only searches for pairs at distances [IS-readLength, IS+readLength]
    //TODO: come up with way of gauging variance in IS and redefine this based on that
    // std::reverse(leftCand.begin(), leftCand.end());
    // auto rightStart = rightCand.begin();
    // for(auto itL = leftCand.begin(); itL != leftCand.end(); itL++){
    //     while(rightStart != rightCand.end() && rightStart->distance + itL->distance < insertSize - readLength){
    //         rightStart++;
    //     }
    //     if(rightStart == rightCand.end()) {break;}
    //     for(auto itR = rightStart; itR != rightCand.end() 
    //         && itR->distance + itL->distance < insertSize + readLength; itR++){
    //         counter++;
    //         if(pair_filter->containsPair(JuncPair(itL->kmer, itR->kmer))){
    //             //std::cout << "Distance: " << itL->distance + itR->distance << "\n";
    //             score += 1;
    //         } 
    //     }
    // }

    //std::cout << "Total pairs: " << rightCand.size()*leftCand.size() << ", counter: " << counter << "\n";
    //std::cout << "Junctions: " << leftCand.size() << "," << rightCand.size() << ". Score: " << score << "\n";
    return score; //getTailBound(leftCand.size()*rightCand.size(),fpRate, score);
}

bool ContigGraph::isBubble(ContigNode* node){
    if (node->numPathsOut() == 2){ // TODO: generalize to more than 2
        std::vector<int> inds = node->getIndicesOut();
        std::vector<ContigNode*> far_nodes;
        for(int i = 0; i != inds.size(); i++) {
            Contig * tig = node->contigs[inds[i]];
            far_nodes.push_back(tig->otherEndNode(node));
        }
        if (allSame(far_nodes)){
            return true;
        }
    }
    return false;
}

std::list<Contig*> ContigGraph::getPathIfSimpleBulge(ContigNode* node, int max_dist){
    std::list<Contig*>  path = {};
    std::list<Contig*>  cand_path = {};
    if (node->numPathsOut() == 2){

        std::vector<int> inds = node->getIndicesOut();
        // target is far end of longer out contig
        ContigNode* target = node->contigs[inds[1]]->otherEndNode(node);
        int dist = node->contigs[inds[1]]->getSeq().length();
        int max_ind = 1;
        int min_ind = 0;
        if (node->contigs[inds[0]]->getSeq().length() > node->contigs[inds[1]]->getSeq().length()){
            target = node->contigs[inds[0]]->otherEndNode(node);
            dist = node->contigs[inds[0]]->getSeq().length();
            max_ind = 0;
            min_ind = 1;
        }
        if (dist > max_dist ){ return path; }

        // try to reach target starting from other contig
        Contig * start = node->contigs[inds[min_ind]];
        if (start->otherEndNode(node) == target){ 
            path.push_back(start); //NodeQueueEntry(node, min_ind, 0));
            return path;
        }
        // BFS from start up to d or max_dist
        else{
            cand_path = node->doPathsConvergeNearby(inds[max_ind], inds[min_ind], max_dist);
            if (cand_path.size()==0) return path;
            else{
                // NodeQueueEntry entry = *cand_path.end(); 
                int target_dist = 0;
                for (auto it = cand_path.begin(); it!= cand_path.end(); ++it){
                    Contig * tig = (*it);
                    target_dist += tig->getSeq().length();
                }
                // int target_dist = entry.startDist + entry.node->contigs[entry.index].getSeq().length();
                int diff = std::abs(target_dist - dist);
                if ( diff < 4 || diff <= 0.05 * dist) return cand_path;    
            }
            
        }
      
    }
    return path;

}


int ContigGraph::collapseBulges(int max_dist){
    printf("Collapsing simple bulges. Starting with %d nodes. \n", nodeMap.size());
    int numDeleted = 0;
    std::set<kmer_type> seenKmers = {};

    // for(auto it = nodeMap.begin(); it != nodeMap.end(); ){
    it = nodeMap.begin();
    while(it!=nodeMap.end()){
        ContigNode* node = &it->second;
        ContigNode* far_node;
        kmer_type kmer = it->first;
        std::list<Contig*> path = getPathIfSimpleBulge(node, max_dist);
        if (path.size() > 0 && seenKmers.find(kmer) == seenKmers.end()){
            // P to be merged into Q
            // if lengths differ but have same next node, P is lower coverage extension
            // if lengths differ but have different next node, P is longer extension
            // P is lower coverage extension if lengths equal
            // ignore if both coverage and length identical
            std::vector<int> inds = node->getIndicesOut();
            std::vector<double> covs;
            std::vector<int> lengths;
            for(int i = 0; i != inds.size(); i++) {
                Contig * tig = node->contigs[inds[i]];
                covs.push_back(tig->getAvgCoverage());
                lengths.push_back(tig->getSeq().length());
            }
            int P_index, Q_index;
            
            if (lengths[0]==lengths[1]){
                if (covs[0]==covs[1]){
                    it++;
                    continue; 
                }
                auto result = std::minmax_element(covs.begin(), covs.end());
                P_index = inds[(result.first - covs.begin())];
                Q_index = inds[(result.second - covs.begin())];
                
            }else{
                if (node->contigs[inds[0]]->otherEndNode(node) == node->contigs[inds[1]]->otherEndNode(node)){
                    auto result = std::minmax_element(covs.begin(), covs.end());
                    P_index = inds[(result.first - covs.begin())];
                    Q_index = inds[(result.second - covs.begin())];
                }
                else{
                    auto result = std::minmax_element(lengths.begin(), lengths.end());   
                    Q_index = inds[(result.first - lengths.begin())];
                    P_index = inds[(result.second - lengths.begin())];  
                }
            }

            // From here on we break stuff...
            // INCOMPLETE - TODO - reassign coverage to max_contig (Q as a path)
            Contig* P = node->contigs[P_index];
            Contig* Q = node->contigs[Q_index];
            // printf("P cov %f, length %d : Q cov %f, length %d\n", P->getAvgCoverage(), P->getSeq().length(), Q->getAvgCoverage(), Q->getSeq().length());            
            far_node = P->otherEndNode(node);
            kmer_type far_kmer;
            if (far_node){
                far_kmer = far_node->getKmer();
            }
            int P_cov = std::round(P->getAvgCoverage()); 
            // for (std::list<Contig*>::iterator *it = path.begin(); it != path.end(); ++it){
            for(auto it = path.begin(); it != path.end(); it++){
                Contig* contig = *it;
                // if(contig->node1_p || contig->node2_p){
                if(contig->node1_p){
                    contig->node1_p->cov[contig->ind1] += P_cov;
                }
                if(contig->node2_p){
                    contig->node2_p->cov[contig->ind2] += P_cov;
                }
            }
            deleteContig(P);

            if (testAndCutIfDegenerate(node)) seenKmers.insert(kmer);
            if(node->numPathsOut() == 1){
                collapseNode(node, kmer);         
                seenKmers.insert(kmer);     
            }
            
            if (far_node){
                if (testAndCutIfDegenerate(far_node)) seenKmers.insert(far_kmer);
                if(far_node->numPathsOut() == 1){
                    collapseNode(far_node, far_kmer);         
                    seenKmers.insert(far_kmer);
                }             
            }
            numDeleted++;
            it++;
        }
        else{
            it++;
        }
    }

    for (auto k_it = seenKmers.begin(); k_it != seenKmers.end(); k_it++){
        nodeMap.erase(*k_it);
    }

    printf("Collapsed %d bulge contigs.\n", numDeleted);
    printf("%d nodes left in node map \n", nodeMap.size());

    return numDeleted; 
}

int ContigGraph::popBubblesByCoverageRatio(){
    printf("Popping bubbles. Starting with %d nodes. \n", nodeMap.size());

    int numDeleted = 0;
    for(auto it = nodeMap.begin(); it != nodeMap.end(); it++){
        
        ContigNode* node = &it->second;
        if (isBubble(node)){
            std::vector<int> inds = node->getIndicesOut();
            std::vector<double> covs;
            for(int i = 0; i != inds.size(); i++) {
                Contig * tig = node->contigs[inds[i]];
                covs.push_back(tig->getAvgCoverage());
            }
            auto result = std::minmax_element(covs.begin(), covs.end());
            int min_contig_index = inds[(result.first - covs.begin())];
            int max_contig_index = inds[(result.second - covs.begin())];
            double min_val = covs[(result.first - covs.begin())];
            double max_val = covs[(result.second - covs.begin())];
            // std::cout << "min element at: " << (result.first - covs.begin()) << ", value is "<< covs[(result.first - covs.begin())] <<'\n';
            // std::cout << "max element at: " << (result.second - covs.begin()) << ", value is "<< covs[(result.second - covs.begin())] <<'\n';
            
            Contig * source = node->contigs[4];
            ContigNode * far_node = node->contigs[inds[0]]->otherEndNode(node);
            Contig * target = far_node->contigs[4];
            if (!(source && target)) {    
                // printf("bubble source or target missing.\n");
                continue;
            }

            // TODO: would assigning coverage of deleted contig
            // to remaining contig help downstream (e.g., disentangle)?
            if(source->getSeq().length() >= 500 && target->getSeq().length() >= 500){
                numDeleted++;
                deleteContig(node->contigs[min_contig_index]);
            }
            else{
                double ratio = max_val / min_val;
                if (ratio >= 3.0){
                    numDeleted++;
                    deleteContig(node->contigs[min_contig_index]);
                } 
            }
        }

    }  
    printf("Done popping %d bubble contigs.\n", numDeleted);
    return numDeleted;  
}


int ContigGraph::disentangle(Bloom* pair_filter, int insertSize, bool local_juncs){
    int disentangled = 0;
    double fpRate = pow(pair_filter->weight(), pair_filter->getNumHash());
    bool operationDone = false;
    std::cout << "About to disentangle.  Starting with " << nodeMap.size() << " nodes.\n";
    //looks through all contigs adjacent to nodes
    for(auto it = nodeMap.begin(); it != nodeMap.end(); ){
        ContigNode* node = &it->second;
        //printf("Testing %s\n", print_kmer(node->getKmer()));
        kmer_type kmer = it->first;
        if(node->numPathsOut() == 2){
            Contig* contig = node->contigs[4];
            ContigNode* backNode = contig->otherEndNode(node);
            
            // test for adjacent outwards facing nodes
            if(backNode && node != backNode && backNode->numPathsOut() == 2 && backNode->indexOf(contig) == 4){
                int a,b,c,d;
                
                for (int orientation = 1; orientation < 5; orientation++){ // change orienation instead of calling disentangle with different order 
                    if (operationDone) {
                        operationDone = false;
                        break;
                    } // to avoid trying to disentangle something we already have

                    if (orientation % 2 == 1) {a = backNode->getIndicesOut()[0], b = backNode->getIndicesOut()[1];}
                    else {b = backNode->getIndicesOut()[0], a = backNode->getIndicesOut()[1];}
                    if (orientation > 2){d = node->getIndicesOut()[0], c = node->getIndicesOut()[1];}
                    else {c = node->getIndicesOut()[0], d = node->getIndicesOut()[1];}
                    
                    Contig* contig_a = backNode->contigs[a]; 
                    Contig* contig_b = backNode->contigs[b];
                    Contig* contig_c = node->contigs[c];
                    Contig* contig_d = node->contigs[d];
                    std::list<JuncResult> A,B,C,D;

                    ContigNode* nodeA = contig_a->otherEndNode(backNode);
                    ContigNode* nodeB = contig_b->otherEndNode(backNode);
                    ContigNode* nodeC = contig_c->otherEndNode(node);
                    ContigNode* nodeD = contig_d->otherEndNode(node);             

                    if (local_juncs){
                        A = contig_a->getJuncResults(contig_a->getSide(backNode, a), 0, contig_a->getSeq().length());
                        B = contig_b->getJuncResults(contig_b->getSide(backNode, b), 0, contig_b->getSeq().length());
                        C = contig_c->getJuncResults(contig_c->getSide(node, c), 0, contig_c->getSeq().length());
                        D = contig_d->getJuncResults(contig_d->getSide(node, d), 0, contig_d->getSeq().length());

                    }
                    else{
                        A = backNode->getPairCandidates(a, insertSize);
                        B = backNode->getPairCandidates(b, insertSize);
                        C = node->getPairCandidates(c, insertSize);
                        D = node->getPairCandidates(d, insertSize);
                    }
                    double scoreAC = getScore(A,C, pair_filter, fpRate, insertSize);
                    double scoreAD = getScore(A,D, pair_filter, fpRate, insertSize);
                    double scoreBC = getScore(B,C, pair_filter, fpRate, insertSize);
                    double scoreBD = getScore(B,D, pair_filter, fpRate, insertSize);
                                    
                    // std::list<int> tig_lengths = {contig_a->getSeq().length(), contig_b->getSeq().length(), contig_c->getSeq().length(), contig_d->getSeq().length()};
                    // std::list<int> scores = {scoreAD,scoreBC,scoreAC,scoreBD};
                    // if (orientation==1 && (*std::max_element(scores.begin(),scores.end())>0 || contig->getSeq().length()>2*insertSize)){
                    std::cout << contig << ", contig len " << contig->getSeq().length() << ", contig cov: " << contig->getAvgCoverage() << "\n";
                    std::cout << "lenA: " << contig_a->getSeq().length() << ", lenB: "<< contig_b->getSeq().length() << ", lenC: " << contig_c->getSeq().length() << ", lenD: "<< contig_d->getSeq().length() <<'\n';
                    std::cout << "covA: " << contig_a->getAvgCoverage() << ", covB: "<< contig_b->getAvgCoverage() << ", covC: " << contig_c->getAvgCoverage() << ", covD: "<< contig_d->getAvgCoverage() <<'\n';                
                    std::cout << "scoreAD: " << scoreAD << ", scoreBC: "<< scoreBC << ", scoreAC: " << scoreAC << ", scoreBD: "<< scoreBD <<'\n';
                    
                    
                    // all distinct --> roughly linear regions when all are distinct
                    // also treat double-bubble and bubble adjacent to junction in same way
                    if((allDistinct(std::vector<Contig*>{contig, contig_a, contig_b, contig_c, contig_d}) && (nodeA!=nodeB && nodeC!=nodeD)) ){
                        // (allDistinct(std::vector<ContigNode*> {backNode, node, nodeA, nodeC}) && (nodeA==nodeB || nodeC==nodeD) &&
                        if (orientation > 2){continue;}
                        if( (std::min(scoreAC,scoreBD) > 0 && std::max(scoreAD,scoreBC) == 0)|| 
                        (std::min(scoreAC , scoreBD) > 1 && ((scoreAD == 0 && scoreBC == 1) || (scoreAD == 1 && scoreBC == 0)) )){
                        // here we try to split by junction pairs - using roughly XOR logic
                            disentanglePair(contig, backNode, node, a, b, c, d);
                            operationDone = true;
                        }
                                             
                        // else if (std::max(scoreAD , scoreBC) == 0 && std::max(scoreAC , scoreBD) == 0 && contig->getSeq().length() >= 2*insertSize){
                        // // here we try to split by coverage ratio alone, since the contig is too long  
                        // // for there to be junction pair links
                        //    double covA = contig_a->getAvgCoverage();
                        //    double covB = contig_b->getAvgCoverage();
                        //    double covC = contig_c->getAvgCoverage();
                        //    double covD = contig_d->getAvgCoverage();
                           
                        //    double L_max = std::max(covA, covB);
                        //    double L_min = std::min(covA, covB);
                        //    int L_arg_max, L_arg_min, R_arg_max, R_arg_min;
                        //    if (L_max == covA) {
                        //         L_arg_max = a;
                        //         L_arg_min = b;
                        //     }
                        //     else{
                        //         L_arg_max = b;
                        //         L_arg_min = a;   
                        //     }
                        //     auto R_max = std::max(covC, covD);
                        //     auto R_min = std::min(covC, covD);
                        //     if (R_max == covC) {
                        //         R_arg_max = c;
                        //         R_arg_min = d;
                        //     }
                        //     else{
                        //         R_arg_max = d;
                        //         R_arg_min = c;   
                        //     }
                            
                        //     if(L_max/L_min >=3 && R_max/R_min >=3 && (std::max(L_max/R_max, R_max/L_max)<=1.15 || std::max(L_min/R_min, R_min/L_min)<=1.15)){
                        //         disentanglePair(contig, backNode, node, L_arg_max,L_arg_min,R_arg_max,R_arg_min);
                        //         operationDone = true;
                                
                        //     }  
                        // }

                        if (operationDone){
                            it = nodeMap.erase(it);
                            if(it != nodeMap.end()){
                                if(backNode->getKmer() == it->first){
                                    it++;
                                }
                            }
                            nodeMap.erase(backNode->getKmer());
                            disentangled++;
                            std::cout << "made decision\n";    
                            continue;
                        }

                    }


                    else{ // not all distinct --> usually some looping or bubble on either side
                        if (orientation == 1) {std::cout << "not all distinct, " << contig << "\n";}
                        // take care of each case separately

                        if (nodeA==node && nodeC==backNode && 
                            nodeB != node && nodeB != backNode &&
                            nodeD != node && nodeD != backNode
                            ){
                            
                            if(( (scoreAD>0 || scoreBC>0) && scoreBD>0 && (contig->getSeq().length() + contig_a->getSeq().length()) <= insertSize)||
                                (std::min(scoreAD , scoreBC) > 0 && scoreBD == 0) ){
                                // II. loop - genomic repeat                            
                                
                                Contig* contigBRCRD = contig_b->concatenate(contig, contig_b->getSide(backNode), contig->getSide(backNode));
                                contigBRCRD = contigBRCRD->concatenate(contig_c, contigBRCRD->getSide(node), contig_c->getSide(node));
                                contigBRCRD = contigBRCRD->concatenate(contig, contigBRCRD->getSide(backNode), contig->getSide(backNode));
                                contigBRCRD = contigBRCRD->concatenate(contig_d, contigBRCRD->getSide(node), contig_d->getSide(node));
                                if(nodeB){
                                    nodeB->replaceContig(contig_b, contigBRCRD);
                                }
                                if(nodeD){
                                    nodeD->replaceContig(contig_d, contigBRCRD);
                                }
                                if(!nodeB && !nodeD){
                                    isolated_contigs.push_back(*contigBRCRD);
                                }

                                operationDone = true;
                            }                        
                        }
                        if (operationDone){
                            it = nodeMap.erase(it);
                            if(it != nodeMap.end()){
                                if(backNode->getKmer() == it->first){
                                    it++;
                                }
                            }
                            nodeMap.erase(backNode->getKmer());
                            disentangled++;
                            std::cout << "made decision\n";    
                            continue;
                        }

                        // if (nodeA==node && nodeD==backNode && 
                        //     nodeB != node && nodeB != backNode &&
                        //     nodeC != node && nodeC != backNode
                        //     ){

                        //     if((scoreAC>0 && scoreAD>0 && scoreBC>0 && scoreBD>0 && (contig->getSeq().length() + contig_a->getSeq().length()) <= insertSize)||
                        //         (std::min(scoreAC , scoreBD) > 0 && scoreBC == 0) ){                        
                        //         // II. loop - genomic repeat                            
                                
                        //         Contig* contigBRDRC = contig_b->concatenate(contig, contig_b->getSide(backNode), contig->getSide(backNode));
                        //         contigBRDRC = contigBRDRC->concatenate(contig_d, contigBRDRC->getSide(node), contig_d->getSide(node));
                        //         contigBRDRC = contigBRDRC->concatenate(contig, contigBRDRC->getSide(backNode), contig->getSide(backNode));
                        //         contigBRDRC = contigBRDRC->concatenate(contig_c, contigBRDRC->getSide(node), contig_c->getSide(node));
                        //         if(nodeB){
                        //             nodeB->replaceContig(contig_b, contigBRDRC);
                        //         }
                        //         if(nodeC){
                        //             nodeC->replaceContig(contig_d, contigBRDRC);
                        //         }
                        //         if(!nodeB && !nodeC){
                        //             isolated_contigs.push_back(*contigBRDRC);
                        //         }

                        //         it = nodeMap.erase(it);
                        //         if(it != nodeMap.end()){
                        //             if(backNode->getKmer() == it->first){
                        //                 it++;
                        //             }
                        //         }
                        //         nodeMap.erase(backNode->getKmer());
                        //         disentangled++;
                        //         // disentanglementCleanup(backNode, disentangled);
                        //         std::cout << "genomic repeat loop\n";                        
                        //         continue;
                        //     }


                        // }

                        // III. single hairpin

                        // IV. 

                    }
                    if (orientation==4 && !operationDone) {std::cout << "no decision\n";}

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
void ContigGraph::disentanglePair(Contig* contig, ContigNode* backNode, ContigNode* forwardNode, 
    int a, int b, int c, int d){
    Contig* contigA = backNode->contigs[a];
    Contig* contigB = backNode->contigs[b];
    Contig* contigC = forwardNode->contigs[c];
    Contig* contigD = forwardNode->contigs[d];

    ContigNode* nodeA = contigA->otherEndNode(backNode);
    ContigNode* nodeB = contigB->otherEndNode(backNode);
    ContigNode* nodeC = contigC->otherEndNode(forwardNode);
    ContigNode* nodeD = contigD->otherEndNode(forwardNode);

    // if(!allDistinct({backNode, forwardNode, nodeA, nodeB, nodeC, nodeD})){
    //     std::cout << "not all distinct\n";
    //     // return false;
    // }
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
    // return true;
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
            collapseNode(node, kmer);
            if(!nodeMap.erase(kmer)) printf("ERROR: tried to erase node %s but there was no node.\n", print_kmer(kmer));
        }
    }

    printf("Done collapsing %d nodes.\n", numCollapsed);
    return numCollapsed;
}

void ContigGraph::collapseNode(ContigNode * node, kmer_type kmer){
    Contig* backContig = node->contigs[4];
    Contig* frontContig;
    if(!backContig){
        printf("WTF no back\n");
        return;
    }
    int fronti = 0;
    for(int i = 0; i < 4; i++){
        if(node->contigs[i]) {
            frontContig = node->contigs[i];
            fronti = i;
        }
    }
    if(!frontContig){
        printf("WTF no front\n");
    }
    if(frontContig == backContig){ //circular sequence with a redundant node
        frontContig->node1_p=nullptr;
        frontContig->node2_p=nullptr;
        addIsolatedContig(*frontContig);
        delete backContig;
    }
    else { //normal case of collapsing a node between two other nodes
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
        delete backContig;
        delete frontContig;
    }
    // if(!nodeMap.erase(kmer)) printf("ERROR: tried to erase node %s but there was no node.\n", print_kmer(kmer));
}


void ContigGraph::printGraph(string fileName){
    printf("Printing graph from contig graph to fastg with iterator.\n");
    ofstream fastgFile;
    fastgFile.open(fileName);

    ContigIterator* contigIt = new ContigIterator(this);
    int node_contig_count = 0;
    //prints contigs that are adjacent to nodes
    while(contigIt->hasNextContig()){
        Contig* contig = contigIt->getContig();
        printContigFastG(&fastgFile, contig);
        node_contig_count++;
    }
    printf("printed %d node-connected contigs\n", node_contig_count);
    //prints isolated contigs
    for(auto it = isolated_contigs.begin(); it != isolated_contigs.end(); it++){
        Contig contig = *it;
        printContigFastG(&fastgFile, &contig);

    }
    printf("printed %d isolated contigs\n", isolated_contigs.size());
    //printf("Done printing contigs from contig graph.\n");
    fastgFile.close();
    printf("Done printing graph from contig graph iterator.\n");
}

void ContigGraph::printContigFastG(ofstream* fastgFile, Contig * contig){
    *fastgFile << contig->getFastGHeader(true) << "\n";
    *fastgFile << contig->getSeq() << "\n";
    *fastgFile << contig->getFastGHeader(false) << "\n";
    *fastgFile << revcomp_string(contig->getSeq()) << "\n";
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
                    jFile << ">Contig" << lineNum << "\n";
                    lineNum++;
                    // std::cout << contig->seq << "\n";
                    jFile << canon_contig(contig->getSeq() ) << "\n";
                }
            }
        }
    }

    //prints isolated contigs
    for(auto it = isolated_contigs.begin(); it != isolated_contigs.end(); it++){
        Contig contig = *it;
        jFile << ">Contig" << lineNum << "\n";
        lineNum++;
        //printf("Printing isolated contig.\n");
        jFile << canon_contig(contig.getSeq()) << "\n";
    }

    //printf("Done printing contigs from contig graph.\n");
    jFile.close();
    printf("Done printing contigs from contig graph.\n");
}

ContigGraph::ContigGraph(){
    nodeMap = {};
    it = nodeMap.begin();
    isolated_contigs = {};
    nodeVector = {};
}

//Creates a contig node and returns it or returns the already extant one if it's already extant
//Uses coverage info from junction to create the node
ContigNode * ContigGraph::createContigNode(kmer_type kmer, Junction junc){
    return &(nodeMap.insert(std::pair<kmer_type, ContigNode>(kmer, ContigNode(junc))).first->second);
}

Contig* ContigGraph::getLongestContig(){
    ContigIterator* contigIt = new ContigIterator(this);
    Contig* result;
    int maxLength = 0;
    while(contigIt->hasNextContig()){
        Contig* contig = contigIt->getContig();
        int length = contig->length();
        if(length > maxLength){
            result = contig;
            maxLength = length;
        }
    }
    return result;
}