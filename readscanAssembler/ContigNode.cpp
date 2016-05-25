#include <fstream>
#include "ContigNode.h"
#include <time.h>
using std::ofstream;
using std::stringstream;
#include <sstream> //for std::stringstream 
#include <string>  //for std::string

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
    
std::list<JuncResult> ContigNode::getPairCandidates(int index, int maxDist) {
    //std::cout << "Getting candidate pairs.\n";
    clock_t t = clock();
    std::set<kmer_type> seenKmers = {};
    std::deque<NodeQueueEntry> queue= {};
    queue.push_back(NodeQueueEntry(this, index, 0));
    std::list<JuncResult> results = {};

    while (!queue.empty()){
        NodeQueueEntry entry = queue.front();
        queue.pop_front();
        kmer_type unique_kmer = entry.node->getUniqueKmer(entry.index);
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



bool ContigNode::doPathsConvergeNearby(int max_ind, int min_ind, int max_dist){
    ContigNode* target = contigs[max_ind]->otherEndNode(this);
    std::set<kmer_type> seenKmers = {};
    std::deque<NodeQueueEntry> queue= {};
    // start from shorter branch
    // int start_dist = contigs[min_ind]->getSeq().length();
    queue.push_back(NodeQueueEntry(this, min_ind, 0));

    // printf("Max, min: %d %d\n", max_ind, min_ind);

    while (!queue.empty()){
        // printf("queue size entering while: %d\n", queue.size());
        NodeQueueEntry entry = queue.front();
        queue.pop_front();
        // printf("Just popped: %li, %d, %d\n", entry.node, entry.index, entry.startDist);
        
        // if (!entry.node){
        //     printf("no node!\n");
        // }
        // if( !entry.node->contigs[entry.index]){
        //     printf("contig missing at extension\n");   
        // }
        if (!entry.node->contigs[entry.index]->node1_p || 
            !entry.node->contigs[entry.index]->node2_p){
            // printf("contig node pointer not present\n");
            continue;
        }
        if (entry.node->contigs[entry.index]->isDegenerateLoop()){
            // printf("degenerate loop\n");
            continue;
        }
        // if(!entry.node->contigs[4]){
        //     printf("Has no backcontig; won't be able to getKmer()\n");
        // }
        // printf("Before seg fault\n");
        //seg fault definitely on this line
        kmer_type unique_kmer = entry.node->getUniqueKmer(entry.index);
        // printf("After seg fault\n");

        if(seenKmers.find(unique_kmer) == seenKmers.end()){
        
            seenKmers.insert(unique_kmer);

            // printf("unseen kmer\n");
            if (entry.startDist > max_dist){
                // printf("too far\n");
                continue;
            }
            else if (entry.node->contigs[entry.index]->otherEndNode(entry.node)==target){ 
                // printf("found target\n");
                return true;
            }
            else{
                // printf("added neighbors\n");
                entry.addNeighbors(queue); //, true); // pushes to front 
                // printf("queue size after adding neighbors: %d\n", queue.size());
            }
            
        }

   }
   return false;
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
    printf("ERROR: tried to find index of contig that's not present.");
    return 5;
}

void ContigNode::update(int nucExt, Contig* contig){
    contigs[nucExt] = contig;
}

void ContigNode::breakPath(int nucExt){
    cov[nucExt] = 0;
    contigs[nucExt] = nullptr;
}

kmer_type ContigNode::getKmer(){
    for(int i = 4; i >= 0; i--){
        if(contigs[i]){
            return contigs[i]->getNodeKmer(this);
        }
    }
    printf("ERROR: No valid contigs from which to getKmer()\n");
    return 0;
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

std::list<JuncResult> NodeQueueEntry::getJuncResults(int maxDist){
    Contig* contig = node->contigs[index];
     return contig->getJuncResults(contig->getSide(node, index),startDist, maxDist);
}

void NodeQueueEntry::addNeighbors(std::deque<NodeQueueEntry>& queue){
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
                // if (to_back){
                    queue.push_back(NodeQueueEntry(nextNode, 4, startDist + contig->getTotalDistance()));
                    // printf("should add one extension\n");
                // }else{
                //     queue.push_front(NodeQueueEntry(nextNode, 4, startDist + contig->getTotalDistance()));                    
                // }
            }
        }
        else{
            for (int i = 0; i < 4; i++){
                if(nextNode->contigs[i]){
                    // if (to_back){
                        queue.push_back(NodeQueueEntry(nextNode, i, startDist + contig->getTotalDistance()));
                        // printf("should add multiple extensions\n");

                    // }else{
                    //     queue.push_front(NodeQueueEntry(nextNode, i, startDist + contig->getTotalDistance()));                        
                    // }
                }
            }
        }
    }
    
}