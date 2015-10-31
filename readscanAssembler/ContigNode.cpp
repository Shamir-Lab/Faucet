#include <fstream>
#include "ContigNode.h"
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

ContigNode::ContigNode(Node node){
	for(int i  = 0; i < 5; i++){
		cov[i] = node.cov[i];
		contigs[i] = nullptr;
	}
}

ContigNode::ContigNode(){
	for(int i  = 0; i < 5; i++){
		cov[i] = -1;
		contigs[i] = nullptr;
	}	
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
    return contigs[4]->getNodeKmer(this);
}

ContigNode* ContigNode::getNeighbor(int index){
    if(contigs[index]){
        return contigs[index]->otherEndNode(this);
    }
    return nullptr;
}

std::string ContigNode::getString(){
    std::string result = "";
    result += print_kmer(getKmer());
    for(int i = 0; i < 5; i++){
        result += " ";
        result +=  getCoverage(i);
    }
    // for(int i = 0; i < 5; i++){
    //     result += " ";
    //     ContigNode * neighbor = getNeighbor(i);
    //     if(neighbor){
    //         result += print_kmer(neighbor->getKmer());
    //     }
    //     else{
    //         result += "X";
    //     }
    // }
    result += "\n";
}


