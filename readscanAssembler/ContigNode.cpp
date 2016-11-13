#include <fstream>
#include "ContigNode.h"
#include <time.h>
using std::ofstream;
using std::stringstream;
#include <sstream> //for std::stringstream 
#include <string>  //for std::string
#include <unordered_set>




ContigNode::ContigNode(Junction junction){
    for(int i  = 0; i < 4; i++){
        cov[i] = junction.getCoverage(i);
        contigs[i] = nullptr;
    }
    contigs[4] = nullptr;
}

ContigNode::ContigNode(){
    for(int i  = 0; i < 5; i++){
        cov[i] = 0;
        contigs[i] = nullptr;
    }   
}
    
bool ContigNode::isInvertedRepeatNode(){
    std::vector<int> inds = this->getIndicesOut();
    std::unordered_set<Contig *> seenContigs = {};
    for (auto i : inds){
        if(seenContigs.find(contigs[i]) == seenContigs.end()){
            seenContigs.insert(contigs[i]);
        }
        else{
            return true;
        }
    }
    return false;
}

std::list<JuncResult> ContigNode::getPairCandidates(int index, int maxDist) {
    //std::cout << "Getting candidate pairs.\n";
    // clock_t t = clock();
    std::unordered_set<kmer_type> seenKmers = {};
    std::vector<NodeQueueEntry> queue;
    // if (maxDist <= contigs[index]->getSeq().length()){
    //     Contig * contig = contigs[index];
    //     return contig->getJuncResults(contig->getSide(this, index),0, maxDist);
    // }  
    queue.push_back(NodeQueueEntry(this, index, 0));
    std::list<JuncResult> results = {};

    while (!queue.empty()){
        NodeQueueEntry entry = queue.front();
        // queue.pop_front();
        queue.erase(queue.begin());
        kmer_type unique_kmer;
        if (!entry.node->contigs[entry.index]){
            continue; // don't advance if at dead end
        }else {
            unique_kmer = entry.node->getUniqueKmer(entry.index);
        }
        if(seenKmers.find(unique_kmer) == seenKmers.end()){
            seenKmers.insert(unique_kmer);
            if(entry.startDist <= maxDist){
                std::list<JuncResult> newResults = entry.getJuncResults(maxDist);
                results.insert(results.end(), newResults.begin(), newResults.end());
                entry.addNeighbors(queue); //, true);
            }
        }
   }
   //std::cout << "Seen kmers: " << seenKmers.size() << "\n";
   //std::cout << "Results: " << results.size() << "\n";
    //std::cout << "Time: "<< clock()-t << " tics.\n";
    results.sort();

    return results;
}



std::list<Contig*> ContigNode::doPathsConvergeNearby(int max_ind, int min_ind, int max_dist){
    /*
        BFSearches up to max_dist away from the node to verify extensions out of node
        converge to the same node. Returns list of Contig ptrs on Q when paths do converge,
        otherwise returns empty list. 
    */
    ContigNode* target = contigs[max_ind]->otherEndNode(this);
    std::unordered_set<kmer_type> seenKmers = {};
    std::vector<NodeQueueEntry> queue;
    queue.reserve(20);

    // std::vector<NodeQueueEntry> parents = {}; // same nodes are inserted as in queue, but never popped off
    // parents.reserve(20);
    std::list<Contig*> path = {};
    
    queue.push_back(NodeQueueEntry(this, min_ind, 0));
    int pos = 0;
    while (pos < queue.size()){//!queue.empty()){
        NodeQueueEntry entry = queue.at(pos);//front();
        pos++;
        // queue.pop_front();
        // queue.erase(queue.begin());
        kmer_type unique_kmer;
        if (!entry.node->contigs[entry.index]){
            continue; // don't advance if at dead end
        }else {
            unique_kmer = entry.node->getUniqueKmer(entry.index);
        }
        if(seenKmers.find(unique_kmer) == seenKmers.end()){
        
            seenKmers.insert(unique_kmer);

            if (entry.startDist > max_dist){
                continue;
            }
            else if (entry.node->contigs[entry.index]->otherEndNode(entry.node)==target){
                // reconstruct path from parents
                path = entry.reconstructPathFromParents(queue);
                return path; //entry.startDist;
            }
            else{
                // entry.recordParents(parents);
                // entry.addNeighbors(parents);
                entry.addNeighbors(queue); 
            }
            
        }

   }
   return path;
}


bool ContigNode::checkValidity(){
    for(int i = 0; i < 5; i++){
        if(contigs[i]){
            Contig* contig = contigs[i];
            int side = contig->getSide(this, i);
            if(side == 1){
                if(contig->ind1 != i){
                    printf("GRAPHERROR: contig has wrong index.\n");
                    return false;
                }
                if(contig->node1_p != this){
                    printf("GRAPHERROR: contig points to wrong node.\n");
                    return false;
                }
            }
            if(side == 2){
                if(contig->ind2 != i){
                    printf("GRAPHERROR: contig has wrong index.\n");
                    return false;
                }
                if(contig->node2_p != this){
                    printf("GRAPHERROR: contig points to wrong node.\n");
                    return false;
                }
            }
        }
    }
    return true;
}

std::vector<std::pair<Contig*, bool>> ContigNode::getFastGNeighbors(int contigIndex){
    std::vector<std::pair<Contig*, bool>> result = {};
    if(contigIndex == 4){
        for(int i = 0; i < 4; i++){
            if(contigs[i]){
                bool RC = false;
                if(contigs[i]->getSide(this,i) == 2) {
                    RC = true;
                }
                result.push_back(std::pair<Contig*, bool>(contigs[i], RC));
            }
        }
    }
    else{
        if(contigs[4]){
            bool RC = false;
            if(contigs[4]->getSide(this,4) == 2) {
                RC = true;
            }
            result.push_back(std::pair<Contig*, bool>(contigs[4], RC));
        }
    }
    return result;
}

kmer_type ContigNode::getForwardExtension(int index){
    return next_kmer(getKmer(), index, FORWARD);
}

kmer_type ContigNode::getUniqueKmer(int index){
    if(index != 4){
        return getForwardExtension(index);
    }
    else{
        return getKmer();
    }
}

int ContigNode::numPathsOut(){
    int numPaths = 0;
    for(int i = 0; i < 4; i++){
        if(cov[i] > 0){
            numPaths++;
        }
    }
    return numPaths;
}

std::vector<int> ContigNode::getIndicesOut(){
    std::vector<int> paths = {};
    for(int i = 0; i < 4; i++){
        if(cov[i] > 0){
            paths.push_back(i);
        }
    }
    return paths;
}

int ContigNode::getTotalCoverage(){
    return getCoverage(4);
}

int ContigNode::getCoverage(int nucExt){
    if(nucExt < 4){
        return (int)cov[nucExt];
    }
    return (int)cov[0] + (int)cov[1] + (int)cov[2] + (int)cov[3];
}

void ContigNode::setCoverage(Junction junc){
    for(int i = 0; i < 4; i++){
        cov[i] = junc.getCoverage(i);
    }
}

void ContigNode::setCoverage(int nucExt, int coverage){
    cov[nucExt] = coverage;
}

void ContigNode::replaceContig(Contig* oldContig, Contig* newContig){
     for(int i = 0; i < 5; i++){
        if(contigs[i] == oldContig){
            contigs[i] = newContig;
        }
    }
}

int ContigNode::indexOf(Contig* contig){
    for(int i = 0; i < 5; i++){
        if(contigs[i] == contig){
            return i;
        }
    }
    throw std::logic_error("ERROR: tried to find index of contig that's not present.");
    // return 5;
}

void ContigNode::update(int nucExt, Contig* contig){
    contigs[nucExt] = contig;
}

void ContigNode::breakPath(int nucExt){
    cov[nucExt] = 0;
    contigs[nucExt] = nullptr;
}

void ContigNode::clearNode(){
    for (int i=0; i<5; i++){
        this->breakPath(i);
    }
}

kmer_type ContigNode::getKmer(){
    for(int i = 4; i >= 0; i--){
        if(contigs[i]){
            return contigs[i]->getNodeKmer(this);
        }
    }
   // intentionally don't return 0 here because that could be a valid kmer value
    throw std::logic_error("No valid contigs from which to getKmer()");
}

ContigNode* ContigNode::getNeighbor(int index){
    if(contigs[index]){
        return contigs[index]->otherEndNode(this);
    }
    return nullptr;
}

std::string ContigNode::getString(){
    std::stringstream result;
    for(int i = 0; i < 5; i++){
        result <<  (int)getCoverage(i) << " ";
       result << contigs[i] << " ";
    }
    return result.str();
}


NodeQueueEntry::NodeQueueEntry(ContigNode* n, int i, int s){
    node = n;
    index = i;
    startDist = s;
}  

NodeQueueEntry::NodeQueueEntry(){
    node = nullptr;
    index = -1;
    startDist = -1;
}

std::list<JuncResult> NodeQueueEntry::getJuncResults(int maxDist){
    Contig* contig = node->contigs[index];
     return contig->getJuncResults(contig->getSide(node, index),startDist, maxDist);
}

void NodeQueueEntry::addNeighbors(std::vector<NodeQueueEntry>& queue){
    Contig* contig = node->contigs[index];
    // if (node->contigs[index]){
    //     printf("no contig at this index!\n");
    // }
    int otherSide = 3 - contig->getSide(node,index);    
    ContigNode* nextNode = contig->getNode(otherSide);
    int nextIndex = contig->getIndex(otherSide);
    


    if(nextNode){
        if(nextIndex != 4){
            if(nextNode->contigs[4]){
                queue.push_back(NodeQueueEntry(nextNode, 4, startDist + contig->getTotalDistance()));                  
            }
        }
        else{
            for (int i = 0; i < 4; i++){
                if(nextNode->contigs[i]){
                    queue.push_back(NodeQueueEntry(nextNode, i, startDist + contig->getTotalDistance()));
                    
                }
            }
        }
    }
    
}

// use stack of parents to reconstruct path: start from target, get other end node
std::list<Contig*> NodeQueueEntry::reconstructPathFromParents(std::vector<NodeQueueEntry>& parents){
    std::list<Contig*> path = {};
    path.push_front(node->contigs[index]); // this is the target
    NodeQueueEntry *currEntry = this;

    // move along parents vector in reverse order
    // query for other end node using entry's contig index
    // when other end node is current entry's node, 
    // make its entry the current entry, add contig to front of path
    std::cout << "in reconstructPathFromParents\n";
    for (auto it = parents.rbegin(); it != parents.rend(); ++it){
        if (it->node->contigs[it->index]->otherEndNode(it->node) == currEntry->node){
            path.push_front(it->node->contigs[it->index]);
            currEntry = &(*it);
        }
    }


    // path.push_front(node->contigs[index]);
    // ContigNode * currNode;
    // Contig * currContig;
    // std::unordered_map<NodeQueueEntry, NodeQueueEntry>::iterator currentEntry = parents.find(*this);
    // while(currentEntry != parents.end()){
    //     currNode = currentEntry->second.node;
    //     currContig = currNode->contigs[currentEntry->second.index];
    //     path.push_front(currContig);
    //     currentEntry = parents.find(currentEntry->second);
    // }
    return path;
}




// void NodeQueueEntry::recordParents(std::unordered_map<NodeQueueEntry, NodeQueueEntry>& parents){
//     Contig* contig = node->contigs[index];
//     // if (node->contigs[index]){
//     //     printf("no contig at this index!\n");
//     // }
//     int otherSide = 3 - contig->getSide(node,index);    
//     ContigNode* nextNode = contig->getNode(otherSide);
//     int nextIndex = contig->getIndex(otherSide);
    


//     if(nextNode){
//         if(nextIndex != 4){
//             if(nextNode->contigs[4]){
//                 parents[NodeQueueEntry(nextNode, 4, startDist + contig->getTotalDistance())] = *this;                  
//             }
//         }
//         else{
//             for (int i = 0; i < 4; i++){
//                 if(nextNode->contigs[i]){
//                     parents[NodeQueueEntry(nextNode, i, startDist + contig->getTotalDistance())] = *this;
                    
//                 }
//             }
//         }
//     }
    
// }

