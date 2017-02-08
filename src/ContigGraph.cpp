#include "ContigGraph.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
// #include "../utils/sparsepp.h"
// using spp::sparse_hash_map;


unordered_map<kmer_type, ContigNode> *  ContigGraph::getNodeMap(){
// sparse_hash_map<kmer_type, ContigNode> *  ContigGraph::getNodeMap(){
    return &nodeMap;
}

std::vector<Contig> * ContigGraph::getIsolatedContigs(){
    return &isolated_contigs;
}


double get_critical_val(int df, double sig){
    // critical values for t distribution at 0.05 alpha level
    std::unordered_map <int,double> df_critical_vals;
    // sparse_hash_map <int,double> df_critical_vals;
    if (sig == 0.05){
        df_critical_vals[1]=6.314;
        df_critical_vals[2]=2.92;
        df_critical_vals[3]=2.353;
        df_critical_vals[4]=2.132;
        df_critical_vals[5]=2.015;
        df_critical_vals[6]=1.943;
        df_critical_vals[7]=1.895;
        df_critical_vals[8]=1.86;
        df_critical_vals[9]=1.833;
        df_critical_vals[10]=1.812;
        df_critical_vals[11]=1.796;
        df_critical_vals[12]=1.782;
        df_critical_vals[13]=1.771;
        df_critical_vals[14]=1.761;
        df_critical_vals[15]=1.753;
        df_critical_vals[16]=1.746;
        df_critical_vals[17]=1.740;
        df_critical_vals[18]=1.734;
        df_critical_vals[19]=1.729;
        df_critical_vals[20]=1.725;
        df_critical_vals[21]=1.721;
        df_critical_vals[22]=1.717;
        df_critical_vals[23]=1.714;
        df_critical_vals[24]=1.711;
        df_critical_vals[25]=1.708;
        df_critical_vals[26]=1.706;
        df_critical_vals[27]=1.703;
        df_critical_vals[28]=1.701;
        df_critical_vals[29]=1.699;
        df_critical_vals[30]=1.697;

        if (df <=30){
            return df_critical_vals[df];
        }
        else if(df > 30 && df <=40){
            return 1.684;
        }
        else if (df > 40 && df <=50){
            return 1.676;
        }
        else if(df > 50 && df <=60){
            return 1.671;
        }
        else if(df >60 && df <=80){
            return 1.664;
        }
        else if(df >80 && df <=100){
            return 1.660;
        }
        else{
            return 1.645;
        }
    }else if(sig==0.025){
        df_critical_vals[1]=12.706;
        df_critical_vals[2]=4.303;
        df_critical_vals[3]=3.182;
        df_critical_vals[4]=2.776;
        df_critical_vals[5]=2.571;
        df_critical_vals[6]=2.447;
        df_critical_vals[7]=2.365;
        df_critical_vals[8]=2.306;
        df_critical_vals[9]=2.262;
        df_critical_vals[10]=2.228;
        df_critical_vals[11]=2.201;
        df_critical_vals[12]=2.179;
        df_critical_vals[13]=2.160;
        df_critical_vals[14]=2.145;
        df_critical_vals[15]=2.131;
        df_critical_vals[16]=2.12;
        df_critical_vals[17]=2.11;
        df_critical_vals[18]=2.101;
        df_critical_vals[19]=2.093;
        df_critical_vals[20]=2.086;
        df_critical_vals[21]=2.08;
        df_critical_vals[22]=2.074;
        df_critical_vals[23]=2.069;
        df_critical_vals[24]=2.064;
        df_critical_vals[25]=2.06;
        df_critical_vals[26]=2.056;
        df_critical_vals[27]=2.052;
        df_critical_vals[28]=2.048;
        df_critical_vals[29]=2.045;
        df_critical_vals[30]=2.042;

        if (df <=30){
            return df_critical_vals[df];
        }
        else if(df > 30 && df <=40){
            return 2.021;
        }
        else if (df > 40 && df <=60){
            return 2.0;
        }
        else if(df > 60 && df <=80){
            return 1.99;
        }
        else if(df >80 && df <=100){
            return 1.984;
        }
        else if(df >100 && df <=1000){
            return 1.962;
        }
        else{
            return 1.96;
        }
    }else{
        throw std::logic_error("tried to get critical value for significance level that is not available");
    }

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
        ++it;
        nodeMap.erase(kmer);
    }
}

bool ContigGraph::deleteTipsAndClean(){
    bool result = false;
    if(deleteTips() > 0){
        result = true;
    } 
    if (validateNoneCollapsible()){
        result = true; 
    }
    if (removeChimericExtensions(150) > 0){
        result  = true;
    }
    if (validateNoneCollapsible()){
        result = true; 
    }
    return result;
}

bool ContigGraph::breakPathsAndClean(){
    bool result = false;
    if(collapseBulges(150) > 0){
        result = true;
    }
    if (validateNoneCollapsible()){
        result = true; 
    }
    return result;
}

bool ContigGraph::disentangleAndClean(Bloom* pair_filter, double insertSize, double std){
    bool result = false;

    if(disentangleLoopPaths(pair_filter, insertSize, std) > 0){
        result = true;
    }
    if (validateNoneCollapsible()){
        result = true; 
    }
    if(disentangleParallelPaths(pair_filter, insertSize, std) > 0){
        result = true;
    }
    if (validateNoneCollapsible()){
        result = true; 
    }
    return result;
}

bool ContigGraph::cleanGraph(Bloom* short_pair_filter, Bloom* long_pair_filter){

    bool result = false;
    deleteTipsAndClean();

    if(breakPathsAndClean()){
        result = true;
    }

    Contig* longContig = getLongestContig();
    std::pair<double, double> mean_std = longContig->getPairsMeanStd(short_pair_filter);
    double insertSize = mean_std.first; 
    double std = mean_std.second;
    std::cout << "short insert size set to " << insertSize << ", std set to " << std << std::endl;
    if(disentangleAndClean(short_pair_filter, insertSize, std)){
        result = true;
    }

    longContig = getLongestContig();
    mean_std = longContig->getPairsMeanStd(long_pair_filter);
    insertSize = mean_std.first;
    std = mean_std.second;
    std::cout << "long insert size set to " << insertSize << ", std set to " << std << std::endl;
    if(disentangleAndClean(long_pair_filter, insertSize, std)){
        result = true;
    }
    
    return result;
}

bool ContigGraph::checkGraph(){
    printf("Checking node mapped contigs.\n");
    for(auto it = nodeMap.begin(); it != nodeMap.end(); ++it){
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
            // std::cout << "node is " << print_kmer(node->getKmer()) << std::endl;         
            if(node->contigs[i]){
                if(!node->contigs[i]->checkValidity()){
                    // return false;
                }
            }
        }                
    }

    printf("Checking %d isolated contigs.\n", isolated_contigs.size());
    //prints isolated contigs
    for(auto it = isolated_contigs.begin(); it != isolated_contigs.end(); ++it){
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

std::pair <Contig*,Contig*> ContigGraph::getMinMaxForwardExtensions(ContigNode * node, std::string trait){
    std::vector<int> inds = node->getIndicesOut();
    std::vector<double> covs;
    std::vector<int> lengths;
    int min_index, max_index;
    for(int i = 0; i < inds.size(); i++) {
        Contig * tig = node->contigs[inds[i]];
        if (node == node->contigs[inds[i]]->otherEndNode(node)){
            // std::cout << "is an inverted repeat at " << inds[i] << std::endl;
        }
        covs.push_back(tig->getAvgCoverage());
        lengths.push_back(tig->getSeq().length());
    }
    if (trait == "coverage"){
        auto result = std::minmax_element(covs.begin(), covs.end());
        min_index = inds[(result.first - covs.begin())]; // lowest coverage 
        max_index = inds[(result.second - covs.begin())]; // highest coverage 
    }else if(trait == "length"){
        auto result = std::minmax_element(lengths.begin(), lengths.end());
        min_index = inds[(result.first - lengths.begin())]; // shortest 
        max_index = inds[(result.second - lengths.begin())]; // longest 
    }else{
        throw std::logic_error("unknown trait");
    }
    return std::make_pair( node->contigs[min_index], node->contigs[max_index] );
}



bool ContigGraph::isTip(ContigNode* node, int index){
    Contig* contig = node->contigs[index];
    if(contig->getSeq().length() < read_length && !contig->otherEndNode(node)){//} && isLowCovContig(contig)){
        return true;
    }
    return false;
}

bool ContigGraph::isLowCovContig(Contig* contig){
    if(contig->getAvgCoverage() < 3){ //* contig->getSeq().length() < 250){
        return true;
    }  
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
    if(!node->contigs[index]){
        // printf("ERROR: tried to cut nonexistant path.");
        return;
    }
    Contig* contig = node->contigs[index];
    int side = contig->getSide(node, index);
    int otherSide = 3 - side;
    if(contig->node1_p == contig->node2_p){ //to handle inverted repeats and loops
        // std::cout << "410\n";
        int otherIndex = contig->getIndex(otherSide);
        contig->setSide(side, nullptr); //set to point to null instead of the node
        contig->setSide(otherSide, nullptr);
        node->breakPath(index);
        node->breakPath(otherIndex);
    }
    else{    
        contig->setSide(side, nullptr);
        node->breakPath(index);
    }
}


int ContigGraph::deleteIsolatedContigs(){
    int numDeleted = 0;
    printf("Deleting isolated low mass contigs. Starting with %d.\n", isolated_contigs.size());
    for(auto it = isolated_contigs.begin(); it != isolated_contigs.end();){
        Contig* contig = &*it;
        if(isLowMassContig(contig)){
            numDeleted++;
            it = isolated_contigs.erase(it);
        }
        else ++it;
    }
    printf("Deleted %d isolated low mass contigs. %d nodes remain\n", numDeleted, isolated_contigs.size());
    return numDeleted;
}
 

bool ContigGraph::testAndCutIfDegenerate(ContigNode* node){
    
    if(node->numPathsOut() == 0 && !node->contigs[4]){
        return true;
    }
    if(node->numPathsOut() == 0){
        if(node->contigs[4]){
            cutPath(node, 4);
        }
        return true;
    }
    else if(!node->contigs[4]){
        // std::cout << "found no back node\n";
        for(int j = 0; j < 4; j++){
            if(node->contigs[j]){
                // std::cout << "tried to remove contig from degenerate\n";
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
            ++it;
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
    /*
        used to test if extensions out of a node create a simple bulge
        returns list of Contig ptrs exiting Q (extension not to be collapsed) when a simple bulge exists
        otherwise returns an empty list as the path
    */
    std::list<Contig*>  path = {};
    std::list<Contig*>  cand_path;
    std::list<Contig*> alt_path;
    if (node->numPathsOut() == 2){
        // std::cout << "535\n";

        std::vector<int> inds = node->getIndicesOut();
        // target is far end of longer extension
        ContigNode* target = node->contigs[inds[1]]->otherEndNode(node);
        if (!target || !node->contigs[inds[0]]->otherEndNode(node)){ return path; }

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
        if (node->contigs[inds[1]]==node->contigs[inds[0]]){ return path; } // inverted repeat
        // std::cout << "553\n";
        // try to reach target starting from other contig
        Contig * start = node->contigs[inds[min_ind]];
        if (start->otherEndNode(node) == target){ 
            // std::cout << "557\n";

            path.push_back(start); //NodeQueueEntry(node, min_ind, 0));
            return path;
        }
        // BFS from start up to d or max_dist
        else{
            // std::cout << "564\n";
            cand_path = node->doPathsConvergeNearby(inds[max_ind], inds[min_ind], max_dist);
            if (cand_path.size()==0){
                alt_path = node->doPathsConvergeNearby(inds[min_ind], inds[max_ind], max_dist);
                if (alt_path.size()==0) {
                    return path;
                }else{
                    cand_path = alt_path;
                }
            } 
            else{
                // std::cout << "575\n";

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

int ContigGraph::deleteTips(){
    printf("Deleting tips. Starting with %d nodes. \n", nodeMap.size());
    int numDeleted = 0;
    it = nodeMap.begin();
    while(it != nodeMap.end()){  
        ContigNode* node = &it->second;
        kmer_type kmer = it->first;
        Contig * contig;
        for(int i = 0; i < 4; i++){ 
            contig = node->contigs[i];
            if(contig){
                if(isTip(node, i) && node->numPathsOut() > 1){ // just means it's short and has no node at other end
                    cutPath(node,i);   // sets node cov/ptr to 0/null, sets contig's node ptr to null on that side          
                    deleteContig(contig); // sets both node's cov/ptr to 0/null on (when not already set to null on contig), deletes contig object, sets ptr to null
                    numDeleted++;
                }
            }
        }
        if(node->contigs[4]){
            if(isTip(node,4)){ // i = 4
                contig = node->contigs[4];
                cutPath(node,4);
                deleteContig(contig);
            }
        }
        if (isCollapsible(node)){ // left with one extension on each end - redundant node
            collapseNode(node, kmer);
            it = nodeMap.erase(it); 
        }
        else if(testAndCutIfDegenerate(node)){  // one end has no extension - expired 'degenerate' node
            // calls cutpath on opposite end -- 4 when no front, all fronts when no back
            it = nodeMap.erase(it); 
        }
        else{
            ++it;
        }
    }
    printf("Deleted %d tips. %d nodes remain\n", numDeleted, nodeMap.size());
}

bool ContigGraph::isCollapsible(ContigNode * node){
    // collapsible means there are contigs at both ends, 
    // and end nodes are not the same for both contigs
    if(node->numPathsOut() != 1) {return false;}
    if(!node->contigs[4]) {return false;}

    Contig * frontContig;
    for (int i = 0; i < 4; i++){ // find the lone remaining contig
        frontContig = node->contigs[i]; 
        if (frontContig){ break;}
    }

    Contig * backContig = node->contigs[4];
    if (!frontContig || !backContig){ return false; }

    if (frontContig->node1_p && frontContig->node2_p){
        if (frontContig->node1_p == frontContig->node2_p) {return false;}    
    }
    if (backContig->node1_p && backContig->node2_p){
        if (backContig->node1_p == backContig->node2_p) {return false;}
    }
    return true;
    
}

int ContigGraph::removeChimericExtensions(int insertSize){
    printf("Removing chimeric extensions. Starting with %d nodes. \n", nodeMap.size());
    int numDeleted = 0; 
    std::unordered_set<kmer_type> seenKmers = {};

    it = nodeMap.begin();
    while(it!=nodeMap.end()){
        ContigNode* node = &it->second;
        kmer_type kmer = it->first;
        Contig* contig = node->contigs[4];

        // try to collapse highest coverage extension Q onto lowest coverage ext. P
        if (node->numPathsOut() > 1 ){ 
            std::pair <Contig *, Contig *> Pair = getMinMaxForwardExtensions(node,"coverage");
            Contig * P = Pair.first;
            Contig * Q = Pair.second;   
            ContigNode * far_node = P->otherEndNode(node);

            if (!far_node || far_node == node){
                ++it;
            }else if ( P->getSeq().length() < 150 && 
                        Q->getAvgCoverage()/P->getAvgCoverage() >= 3 && 
                        far_node->indexOf(P)!=4 && far_node->numPathsOut()>1
                    ){         

                kmer_type far_kmer = far_node->getKmer();     
                cutPath(node, node->indexOf(P));    
                cutPath(far_node, far_node->indexOf(P));                                  
                deleteContig(P);

                if(isCollapsible(node)){
                    collapseNode(node, kmer);                                                
                    it = nodeMap.erase(it);                           
                }else{
                    ++it;
                }                                                       
                numDeleted++;
            }else{
                ++it;
            }      
        }
        else if (isCollapsible(node)){ 
            collapseNode(node, kmer);
            it = nodeMap.erase(it); 
        }
        else{
            ++it;
        }
    }
    printf("removed %d chimeric contigs.\n", numDeleted);
    printf("%d nodes left in node map \n", nodeMap.size());
    return numDeleted; 
}

int ContigGraph::validateNoneCollapsible(){
    printf("Collapsing any remaining collapsible nodes. Starting with %d nodes. \n", nodeMap.size());
    int numDeleted = 0; 
    it = nodeMap.begin();
    while(it!=nodeMap.end()){
        kmer_type kmer = it->first;
        ContigNode* node = &it->second;
        if(node->numPathsOut() == 0 && !node->contigs[4]){
            it = nodeMap.erase(it); 
            numDeleted++;
        }
        else if (isCollapsible(node)){ 
            collapseNode(node, kmer);
            it = nodeMap.erase(it); 
            numDeleted++;
        }else if(testAndCutIfDegenerate(node)){
            it = nodeMap.erase(it); 
            numDeleted++;
        }
        else{
            ++it;
        }
    }
    printf("removed %d collapsible nodes.\n", numDeleted);
    printf("%d nodes left in node map \n", nodeMap.size());
    return numDeleted;     
}

int ContigGraph::collapseBulges(int max_dist){
    printf("Collapsing simple bulges. Starting with %d nodes. \n", nodeMap.size());
    int numDeleted = 0;
    std::unordered_set<kmer_type> seenKmers = {};

    it = nodeMap.begin();
    while(it!=nodeMap.end()){
        // std::cout << "749\n";
        ContigNode* node = &it->second;
        ContigNode* far_node;
        kmer_type kmer = it->first;
        std::list<Contig*> path;
        if (node->numPathsOut() == 2){
            path = getPathIfSimpleBulge(node, max_dist);
        }else{
            path = {};
        }

        // path size is > 0 only if simple bulge found
        if (path.size() > 0 && seenKmers.find(kmer) == seenKmers.end()){
            // std::cout << "756\n";

            std::pair <Contig *, Contig *> Pair = getMinMaxForwardExtensions(node,"coverage");
            Contig * P = Pair.first;
            Contig * Q = Pair.second;    
            Contig * temp;        
            double P_cov = P->getAvgCoverage(); 
            double Q_cov = Q->getAvgCoverage(); 

            if (Q_cov/P_cov < 1.5) {
                // std::cout << "764\n";
                ++it;
                continue;
            }

            if (Q->getSeq().length() == P->getSeq().length()){
                if (Q_cov==P_cov){
                    // std::cout << "771\n";
                    ++it;
                    continue; 
                }
            }else if(Q != *path.begin()) {
                // std::cout << "776\n";
                temp = P;
                P = Q;
                Q = temp;
                // std::cout << "weird - lower coverage path is one left by bulge removal\n";
            }

            // From here on we break stuff...
            // Contig* P = node->contigs[P_index];
            // Contig* Q = node->contigs[Q_index];
            // printf("P cov %f, length %d : Q cov %f, length %d\n", P_cov, P->getSeq().length(), Q_cov, Q->getSeq().length());            
            far_node = P->otherEndNode(node);
            kmer_type far_kmer;
            far_kmer = far_node->getKmer(); // always have two ends for both extensions - guaranteed by getPath method
            
            ContigJuncList  origJuncs, newJuncs;
            // std::cout << "800\n";

            for(auto it = path.begin(); it != path.end(); ++it){
                // std::cout << "796\n";

                Contig* contig = *it;
                origJuncs = contig->contigJuncs;
                // std::cout << "P cov " << P_cov << ", original juncs\n";
                // origJuncs.printJuncValues();

                newJuncs = origJuncs.getShiftedCoverageContigJuncs(P_cov);   
                contig->setContigJuncs(newJuncs);
                // std::cout <<  "updated juncs\n";
                // contig->contigJuncs.printJuncValues();
            }
            deleteContig(P);
            numDeleted++;
            if(isCollapsible(node)){
                // std::cout << "810\n";

                collapseNode(node, kmer);                                                
                it = nodeMap.erase(it);                           
            }else{
                // std::cout << "814\n";
                ++it;
            }        

        }
        else if (isCollapsible(node) ){ 
            // std::cout << "821\n";
            collapseNode(node, kmer);
            it = nodeMap.erase(it); 
            
        }
        else{
            // std::cout << "827\n";
            ++it;
        }
    }

    printf("Collapsed %d bulge contigs.\n", numDeleted);
    printf("%d nodes left in node map \n", nodeMap.size());

    return numDeleted; 
}

int ContigGraph::disentangleParallelPaths(Bloom* pair_filter, double insertSize, double std){
    int disentangled = 0;
    bool operationDone = false;
    double fpRate = pow(pair_filter->weight(), pair_filter->getNumHash());

    if (insertSize <= 0 || std <= 0) return 0;
    std::cout << "About to disentangle parallel paths.  Starting with " << nodeMap.size() << " nodes.\n";

    //looks through all contigs adjacent to nodes
    for(auto it = nodeMap.begin(); it != nodeMap.end(); ){
        ContigNode* node = &it->second;
        kmer_type kmer = it->first;
        Contig* contig = node->contigs[4];
        // std::cout << "887\n";        
        if (isCollapsible(node)){ // left with one extension on each end - redundant node
            collapseNode(node, kmer);
            // std::cout << "890\n";
            it = nodeMap.erase(it);
            continue;
            // std::cout << "872\n";
        }
        else if(testAndCutIfDegenerate(node)){  // either one has or both ends have no extension - expired 'degenerate' node
            // calls cutpath on opposite end -- 4 when no front, all fronts when no back
            // std::cout << "896\n";
            it = nodeMap.erase(it); 
            continue;
        }
        else if (node->numPathsOut() > 2){
            ++it;
            continue;
        }
        else if (!contig){
            ++it;
            continue;   
        }
        else if(node->numPathsOut()==2 && contig->node2_p && contig->node1_p){
            ContigNode* backNode = contig->otherEndNode(node);
            if (!backNode || isCollapsible(backNode)){ 
                // std::cout << "910\n";
        		++it;
                continue;
            }
            if(testAndCutIfDegenerate(backNode)){  
                // std::cout << "915\n";                
                ++it;
                continue;
            }
            int a,b,c,d;
            if (node != backNode && backNode->numPathsOut() == 2 && backNode->indexOf(contig) == 4 &&
                !node->isInvertedRepeatNode() && !backNode->isInvertedRepeatNode()){
                kmer_type back_kmer = backNode->getKmer();
                // std::cout << "923\n";
                for (int orientation = 1; orientation < 3; orientation++){
                    if (orientation % 2 == 1) {a = backNode->getIndicesOut()[0], b = backNode->getIndicesOut()[1];}
                    else {b = backNode->getIndicesOut()[0], a = backNode->getIndicesOut()[1];}
                    c = node->getIndicesOut()[0], d = node->getIndicesOut()[1];
                    
                    Contig* contig_a = backNode->contigs[a]; 
                    Contig* contig_b = backNode->contigs[b];
                    Contig* contig_c = node->contigs[c];
                    Contig* contig_d = node->contigs[d];

                    std::list<JuncResult> A,B,C,D;

                    ContigNode* Nodea = contig_a->otherEndNode(backNode);
                    ContigNode* Nodeb = contig_b->otherEndNode(backNode);
                    ContigNode* Nodec = contig_c->otherEndNode(node);
                    ContigNode* Noded = contig_d->otherEndNode(node);             
                    double scoreAC, scoreAD, scoreBC, scoreBD;

                    int len_a = contig_a->getSeq().length();
                    int len_b = contig_b->getSeq().length();
                    int len_c = contig_c->getSeq().length();
                    int len_d = contig_d->getSeq().length();

                    double cov_a = contig_a->getAvgCoverage();
                    double cov_b = contig_b->getAvgCoverage();
                    double cov_c = contig_c->getAvgCoverage();
                    double cov_d = contig_d->getAvgCoverage();

                    // add 1 to always get at least a flanking junction
                    int stepSize = (int) insertSize + 3*std;
                    A = backNode->getPairCandidates(a, std::min(len_a, stepSize));
                    B = backNode->getPairCandidates(b, std::min(len_b, stepSize));
                    C = node->getPairCandidates(c, std::min(len_c, stepSize));
                    D = node->getPairCandidates(d, std::min(len_d, stepSize));
                
                    scoreAC = getScore(A,C, pair_filter, fpRate, stepSize);
                    scoreAD = getScore(A,D, pair_filter, fpRate, stepSize);
                    scoreBC = getScore(B,C, pair_filter, fpRate, stepSize);
                    scoreBD = getScore(B,D, pair_filter, fpRate, stepSize);

                    
                    if(allDistinct(std::vector<Contig*>{contig, contig_a, contig_b, contig_c, contig_d})){
                       
                        if ((std::min(scoreAC,scoreBD) > 0 && std::max(scoreAD,scoreBC) == 0)||
                            (scoreAC >= 10 && std::max(scoreAD,scoreBC) == 0 && len_a >= 500 && len_c >= 500)){

                            // std::cout << "desired split found\n";
                            if(allDistinct(std::vector<ContigNode*>{node,backNode,Nodea,Nodeb,Nodec,Noded}) ||
                            (Nodea==Nodec && Nodea!=Nodeb && Nodec!=Noded && allDistinct(std::vector<ContigNode*>{node,backNode,Nodeb,Noded})) ||
                            (Nodeb==Noded && Nodea!=Nodeb && Nodec!=Noded && allDistinct(std::vector<ContigNode*>{node,backNode,Nodea,Nodec})) ){
                                operationDone = true;  // everything distinct
                            }
                            else if (Nodea!=Noded && Nodea!=Nodec && Nodeb!=Noded && Nodeb!=Nodec ){
                                if (Nodea && Nodec){
                                    if (Nodea==Nodeb && Nodec==Noded && 
                                        Nodea->indexOf(contig_a) != 4 && Nodea->indexOf(contig_b) != 4 &&
                                        Nodec->indexOf(contig_c) != 4 && Nodec->indexOf(contig_d) != 4){
                                        operationDone = true; // double bubble
                                    }
                                }
                                if (Nodea){
                                    if(Nodea==Nodeb && 
                                        Nodea->indexOf(contig_a) != 4 && Nodea->indexOf(contig_b) != 4){
                                        operationDone = true; // single bubble on back side
                                    }
                                }
                                if (Nodec){
                                    if(Nodec==Noded && 
                                        Nodec->indexOf(contig_c) != 4 && Nodec->indexOf(contig_d) != 4){
                                        operationDone = true; // single bubble on front side
                                    }
                                }                                   
                            }
                        }

                        else{
                            if( (areEquivalentContigCoverages(contig_a->contigJuncs, contig_c->contigJuncs, 0.10) ||
                                areEquivalentContigCoverages(contig_b->contigJuncs, contig_d->contigJuncs, 0.10) ) &&
                                areDifferentialContigCoverages(contig_a->contigJuncs, contig_d->contigJuncs) && 
                                areDifferentialContigCoverages(contig_b->contigJuncs, contig_c->contigJuncs) &&
                                (std::max(cov_a/cov_b, cov_b/cov_a) >= 1.5 || std::max(cov_c/cov_d, cov_d/cov_c) >= 1.5)
                                ){
                                // contig->getSeq().length() > insertSize
                                // std::cout << "split found by coverage\n";
                                operationDone = true;
                            }
                        }    
                    }

                    if (operationDone){
                        disentanglePair(contig, backNode, node, a, b, c, d);
                        disentangled++;
                        for (int i = 0; i< 4; i++){
                            cutPath(backNode, i);
                            cutPath(node, i);
                        }
                        // backNode->clearNode();
                        // node->clearNode();
                        it = nodeMap.erase(it);
                        operationDone = false;
                        break;
                    }
                    else{
                        if (orientation==2){ // made it through inner for without disentanglement
                            ++it;
                        }
                    }
                }
            }
            else{
                ++it;
            }
        }
        else{            
            ++it;            
        }
    }
    printf("Done disentangling %d parallel paths.\n", disentangled);
    return disentangled;
}

int ContigGraph::disentangleLoopPaths(Bloom* pair_filter, double insertSize, double std){
    int disentangled = 0;
    bool operationDone = false;
    double fpRate = pow(pair_filter->weight(), pair_filter->getNumHash());

    if (insertSize <= 0 || std <= 0) return 0;
    std::cout << "About to disentangle loop paths.  Starting with " << nodeMap.size() << " nodes.\n";

    //looks through all contigs adjacent to nodes
    for(auto it = nodeMap.begin(); it != nodeMap.end(); ){
        ContigNode* node = &it->second;
        kmer_type kmer = it->first;
        Contig* contig = node->contigs[4];
        // std::cout << "1055\n";
        if (isCollapsible(node)){ // left with one extension on each end - redundant node
            collapseNode(node, kmer);
            // std::cout << "1058\n";
            it = nodeMap.erase(it);
            continue; 
            // std::cout << "872\n";
        }
        else if(testAndCutIfDegenerate(node)){  // either one has or both ends have no extension - expired 'degenerate' node
            // calls cutpath on opposite end -- 4 when no front, all fronts when no back
            // std::cout << "1064\n";
            it = nodeMap.erase(it); 
            continue;
        }
        else if (node->numPathsOut() > 2){
            ++it;
            continue;
        }
        else if (!contig){
            ++it;
            continue;   
        }
        else if(node->numPathsOut()==2 && contig->node2_p && contig->node1_p){
            ContigNode* backNode = contig->otherEndNode(node);
            if (!backNode || isCollapsible(backNode)){ 
                // std::cout << "1078\n";
                ++it;
                continue;
            }
            if(testAndCutIfDegenerate(backNode)){  
                // std::cout << "1083\n";
                ++it;
                continue;
            }
            int a,b,c,d;
            if (node != backNode && backNode->numPathsOut() == 2 && backNode->indexOf(contig) == 4 &&
                !node->isInvertedRepeatNode() && !backNode->isInvertedRepeatNode()){
                kmer_type back_kmer = backNode->getKmer();
                // std::cout << "1091\n";
                for (int orientation = 1; orientation < 5; orientation++){
                    if (orientation % 2 == 1) {a = backNode->getIndicesOut()[0], b = backNode->getIndicesOut()[1];}
                    else {b = backNode->getIndicesOut()[0], a = backNode->getIndicesOut()[1];}
                    if (orientation > 2){d = node->getIndicesOut()[0], c = node->getIndicesOut()[1];}
                    else {c = node->getIndicesOut()[0], d = node->getIndicesOut()[1];}
                    
                    Contig* contig_a = backNode->contigs[a]; 
                    Contig* contig_b = backNode->contigs[b];
                    Contig* contig_c = node->contigs[c];
                    Contig* contig_d = node->contigs[d];

                    std::list<JuncResult> A,B,C,D;

                    ContigNode* Nodea = contig_a->otherEndNode(backNode);
                    ContigNode* Nodeb = contig_b->otherEndNode(backNode);
                    ContigNode* Nodec = contig_c->otherEndNode(node);
                    ContigNode* Noded = contig_d->otherEndNode(node);             
                    double scoreAC, scoreAD, scoreBC, scoreBD;

                    int len_a = contig_a->getSeq().length();
                    int len_b = contig_b->getSeq().length();
                    int len_c = contig_c->getSeq().length();
                    int len_d = contig_d->getSeq().length();

                    double cov_a = contig_a->getAvgCoverage();
                    double cov_b = contig_b->getAvgCoverage();
                    double cov_c = contig_c->getAvgCoverage();
                    double cov_d = contig_d->getAvgCoverage();

                    // add 1 to always get at least a flanking junction
                    int stepSize = (int) insertSize + 3*std;
                    A = backNode->getPairCandidates(a, std::min(len_a, stepSize));
                    B = backNode->getPairCandidates(b, std::min(len_b, stepSize));
                    C = node->getPairCandidates(c, std::min(len_c, stepSize));
                    D = node->getPairCandidates(d, std::min(len_d, stepSize));
                
                    scoreAC = getScore(A,C, pair_filter, fpRate, stepSize);
                    scoreAD = getScore(A,D, pair_filter, fpRate, stepSize);
                    scoreBC = getScore(B,C, pair_filter, fpRate, stepSize);
                    scoreBD = getScore(B,D, pair_filter, fpRate, stepSize);

                    // std::cout << contig << ", contig len " << contig->getSeq().length() << ", contig cov: " << contig->getAvgCoverage() << ", insert size is " << insertSize << "\n";
                    // std::cout << "lenA: " << len_a << ", lenB: "<< len_b << ", lenC: " << len_c << ", lenD: "<< len_d <<'\n';
                    // std::cout << "covA: " << cov_a << ", covB: "<< cov_b << ", covC: " << cov_c << ", covD: "<< cov_d <<'\n';                
                    // std::cout << "scoreAD: " << scoreAD << ", scoreBC: "<< scoreBC << ", scoreAC: " << scoreAC << ", scoreBD: "<< scoreBD <<'\n';
                    // std::cout << "size A: " << A.size() << ", size B: "<< B.size() << ", size C: " << C.size() << ", size D: "<< D.size() <<'\n';                

                    if (Nodea==node && Nodec==backNode && Nodeb != node && Nodeb != backNode && Noded != node && Noded != backNode){
                           
                        if((scoreAD>0 || scoreBC>0) ){ // && scoreAC == 0 && scoreBD == 0
                            // loop - genomic repeat                            

                            ContigJuncList origJuncs = contig->contigJuncs;   
                            ContigJuncList newJuncs;                             
                            double BC_weight = (cov_b*len_b + cov_c*len_c) / (len_b + len_c);
                            double CD_weight = (cov_c*len_c + cov_d*len_d) / (len_c + len_d);
                            double scale_factor_BC = BC_weight  / (BC_weight + CD_weight);
                            double scale_factor_CD = 1 - scale_factor_BC; 
                            newJuncs = origJuncs.getScaledContigJuncs(scale_factor_BC);   
                            contig->setContigJuncs(newJuncs);
                            Contig* contigBRCRD = contig_b->concatenate(contig, contig_b->getSide(backNode), contig->getSide(backNode));
                            contigBRCRD = contigBRCRD->concatenate(contig_c, contigBRCRD->getSide(node), contig_c->getSide(node));
                            newJuncs = origJuncs.getScaledContigJuncs(scale_factor_CD);   
                            contig->setContigJuncs(newJuncs);
                            contigBRCRD = contigBRCRD->concatenate(contig, contigBRCRD->getSide(backNode), contig->getSide(backNode));
                            contigBRCRD = contigBRCRD->concatenate(contig_d, contigBRCRD->getSide(node), contig_d->getSide(node));
                            if(Nodeb){
                                Nodeb->replaceContig(contig_b, contigBRCRD);
                            }
                            if(Noded){
                                Noded->replaceContig(contig_d, contigBRCRD);
                            }
                            if(!Nodeb && !Noded){
                                isolated_contigs.push_back(*contigBRCRD);
                            }
                            // std::cout << "split found for loop\n";
                            operationDone = true;
                        }      
                    }
                    if (operationDone){
                        // std::cout << "970\n";
                        disentanglePair(contig, backNode, node, a, b, c, d);
                        disentangled++;
                        for (int i = 0; i< 4; i++){
                            cutPath(backNode, i);
                            cutPath(node, i);
                        }
                        // backNode->clearNode();
                        // node->clearNode();
                        it = nodeMap.erase(it);
                        operationDone = false;
                        break;
                    }
                    else{
                        // std::cout << "985\n";
                        if (orientation==4){ // made it through inner for without disentanglement
                            ++it;
                        }
                    }
                }
            }
            else{
                ++it;
            }
        }
        else{            
            ++it;            
        }
    }
    printf("Done disentangling %d loop paths.\n", disentangled);
    return disentangled;
}


bool ContigGraph::areEquivalentContigCoverages(ContigJuncList A, ContigJuncList B, double frac){
    // two one sided T-tests: frac is portion of a's mean allowed to vary
    // 0.05 significance level 
    double ma = A.getAvgCoverage();
    double mb = B.getAvgCoverage();
    double sa = A.getCoverageSampleVariance();
    double sb = B.getCoverageSampleVariance();
    int na = A.size();
    int nb = B.size();
    int df = na + nb - 2;
    if (!((sa > 0 || sb > 0) && (na > 2 && nb > 2))){ return false; }
    double diff = ma - mb;
    double thresh_hi = frac*ma;
    double thresh_lo = -frac*ma;
    if (nb > na){
        thresh_hi = frac*mb;
        thresh_lo = -frac*mb;
    }
    double two_samp_var = pow(pow(sa,2)/na + pow(sb,2)/nb , 0.5);
    double c_t = get_critical_val(df, 0.05);
    double t_lo = (diff - thresh_lo) / two_samp_var;
    double t_hi = (diff - thresh_hi) / two_samp_var; 
    if (t_lo >= c_t && t_hi <= -c_t){
        // std::cout << "returned equivalent\n";
        return true;
    }
    else {
        // std::cout << "returned not equivalent\n"; 
        return false;
    }
}

bool ContigGraph::areDifferentialContigCoverages(ContigJuncList A, ContigJuncList B){
    // two tailed T-test: 
    // 0.05 significance (alpha) level (halved b/c two-tailed) 
    double ma = A.getAvgCoverage();
    double mb = B.getAvgCoverage();
    double sa = A.getCoverageSampleVariance();
    double sb = B.getCoverageSampleVariance();
    int na = A.size();
    int nb = B.size();
    if (!((sa > 0 || sb > 0) && (na > 2 && nb > 2))){ return false; }
    int df = std::round(pow(sa/na + sb/nb,  2) / (pow(sa/na, 2)/(na-1) + pow(sb/nb, 2)/(nb-1)));
    double two_samp_var = pow(pow(sa,2)/na + pow(sb,2)/nb , 0.5);
    double c_t = get_critical_val(df, 0.025);
    double T = std::abs((ma - mb)/two_samp_var);
    
    if (T > c_t){
        // std::cout << "returned differential\n";
        return true;
    }
    else {
        // std::cout << "returned not differential\n"; 
        return false;
    }
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

    double covA = contigA->getAvgCoverage();
    double covB = contigB->getAvgCoverage();
    double covC = contigC->getAvgCoverage();
    double covD = contigD->getAvgCoverage();

    int lenA = contigA->getSeq().length();
    int lenB = contigB->getSeq().length();
    int lenC = contigC->getSeq().length();
    int lenD = contigD->getSeq().length();

    ContigJuncList origJuncs = contig->contigJuncs;
    ContigJuncList newJuncs;
    
    //new ContigJuncList(contig->contigJuncs->seq, contig->contigJuncs->distances, contig->contigJuncs->coverages);
    double AC_weight = (covA*lenA + covC*lenC) / (lenA + lenC);
    double BD_weight = (covB*lenB + covD*lenD) / (lenB + lenD);
    double scale_factor_AC = AC_weight  / (AC_weight + BD_weight);
    double scale_factor_BD = 1 - scale_factor_AC; 
    newJuncs = origJuncs.getScaledContigJuncs(scale_factor_AC);   
    // std::cout << "AC factor " << scale_factor_AC << ", BD factor " << scale_factor_BD <<  ", original juncs\n";
    // origJuncs.printJuncValues();
      
    contig->setContigJuncs(newJuncs);
    // std::cout << "after first scaling\n";
    // newJuncs.printJuncValues();

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

    // clear coverages in newJuncs, set to original values scaled second way
    newJuncs = origJuncs.getScaledContigJuncs(scale_factor_BD);     
    contig->setContigJuncs(newJuncs);
    // std::cout << "after second scaling\n";
    // newJuncs.printJuncValues();

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
        ++it;
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
    Contig* frontContig = nullptr;
    if(!backContig){
        // printf("WTF no back\n");
        // std::cout << "1402\n";
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
        // printf("WTF no front\n");
        // std::cout << "1414\n";
        return;
    }

    if(frontContig == backContig){ //circular sequence with a redundant node
        // std::cout << "1419\n";
        // frontContig->node1_p=nullptr;
        // frontContig->node2_p=nullptr;
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
            // std::cout << "1439\n";
            addIsolatedContig(*newContig);
            delete newContig;
        } 
        delete backContig;
        delete frontContig;
    }
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
    for(auto it = isolated_contigs.begin(); it != isolated_contigs.end(); ++it){
        Contig* contig = &*it;
        if (contig->getSeq().length() >= 200){
           printContigFastG(&fastgFile, contig);
        }

    }
    printf("printed %d isolated contigs\n", isolated_contigs.size());
    //printf("Done printing contigs from contig graph.\n");
    fastgFile.close();
    printf("Done printing graph from contig graph iterator.\n");
}

void ContigGraph::printContigFastG(std::ostream* fastgFile, Contig * contig){
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
    for(auto it = nodeMap.begin(); it != nodeMap.end(); ++it){
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
    for(auto it = isolated_contigs.begin(); it != isolated_contigs.end(); ++it){
        Contig contig = *it;
        if (contig.getSeq().length() >= 200){

            jFile << ">Contig" << lineNum << "\n";
            lineNum++;
            //printf("Printing isolated contig.\n");
            jFile << canon_contig(contig.getSeq()) << "\n";
        }
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
