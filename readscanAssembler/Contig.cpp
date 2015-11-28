#include <fstream>
#include <sstream>
#include <algorithm>    // std::reverse
#include <vector>       // std::vector
#include "Contig.h"
using std::stringstream;
using std::ofstream;

Contig* Contig::concatenate(Contig* otherContig, int thisSide, int otherSide){
	if(thisSide == 1){
		reverse();
	}
	if(otherSide == 2){
		otherContig->reverse();
	}
	return concatenate(otherContig);
}

//utility for linking them if they're both facing "forward"
Contig* Contig::concatenate(Contig* otherContig){
	Contig* result = new Contig();
	result->setEnds(node1_p, ind1, otherContig->node2_p, otherContig->ind2);
	if(seq.length() < sizeKmer){
		printf("ERROR: seq less than k long in Contig::Concatenate.\n");
	}
	result->setSeq(seq.substr(0, seq.length()-sizeKmer) + otherContig->seq);
	result->addJuncDistances(juncDistances);
	result->addJuncDistances(otherContig->juncDistances);
	result->setCoverage(coverageSum + otherContig->coverageSum);
	return result;
}

void Contig::reverse(){
	{ContigNode * temp = node1_p;
		node1_p = node2_p;
		node2_p = temp;}

	{int temp = ind1;
		ind1 = ind2;
		ind2 = temp;}

	seq = revcomp_string(seq);

	std::reverse(juncDistances.begin(), juncDistances.end());
}

void Contig::setEnds( ContigNode* n1, int i1, ContigNode* n2, int i2){
	node1_p = n1;
	node2_p = n2;
	setIndices(i1, i2);
	if(node1_p){
		node1_p->contigs[i1] = this;
	}
	if(node2_p){
		node2_p->contigs[i2] = this;
	}
}

void Contig::setSeq(std::string cont){
	seq = cont;
}

void Contig::setCoverage(int cov){
	coverageSum = cov;
}

void Contig::addCoverage(int cov){
	coverageSum += cov;
}

long int Contig::getAvgCoverage(){
	//printf("Cov: %d\n", coverageSum);
	//printf("Junc dist size: %u\n", juncDistances.size());
	if( juncDistances.size() < 0){
		return -1;
	}
	return (long int) coverageSum / (long int)(juncDistances.size()+1);
}

int Contig::getMinAdjacentCoverage(){
	if(isIsolated()) return 0;
	int min = 1000000;
	if(node1_p){
		min = node1_p->getTotalCoverage();
	}
	if(node2_p){
		min = std::min(min, node2_p->getTotalCoverage());
	}
	return min;
}

float Contig::getMass(){
	return getAvgCoverage()*seq.length();
}

void Contig::addJuncDistances(std::vector<unsigned char> distances){
	for(auto it = distances.begin(); it != distances.end(); it++){
		unsigned char dist = *it;
		juncDistances.push_back(dist);
	}
}

void Contig::addJuncDistance(unsigned char dist){
	juncDistances.push_back(dist);
}

void Contig::setIndices(int i1, int i2){
	ind1 = i1;
	ind2 = i2;
}

int Contig::getMinIndex(){
	return std::min(ind1, ind2);
}

ContigNode* Contig::otherEndNode(ContigNode * oneEnd){
	if(node1_p == oneEnd){
		return node2_p;
	}
	if(node2_p == oneEnd){
		return node1_p;
	}
	printf("ERROR: tried to get other end of a contig, but the given pointer didn't point to either end!.\n");
	std::cout << "node1_p: " << node1_p << " node2_p: " << node2_p << " oneEnd: " << oneEnd << "\n";
	return nullptr;
}

//Assumes the given contig node points to one end of this contig
kmer_type Contig::getNodeKmer(ContigNode * contigNode){
	if(node1_p == contigNode){
		return getSideKmer(1);
	}
	if(node2_p == contigNode){
		return getSideKmer(2);
	}
	printf("ERROR: tried to get the kmer corresponding to a node not adjacent to this contig from this contig.\n");
}

//Gets kmer for node1_p if side == 1, node2_p if side == 2
kmer_type Contig::getSideKmer(int side){
	if(side == 1){
		kmer_type kmer = getKmerFromRead(seq, 0);
		if(ind1 == 4) return revcomp(kmer);
		return kmer;
	}
	if(side == 2){
		kmer_type kmer = getKmerFromRead(seq, seq.length()-sizeKmer);
		if(ind2 == 4) return kmer;
		return revcomp(kmer);
	}
	printf("ERROR: tried to get a kmer corresponding to a side other than one or two from a contig.\n");
}

int Contig::getSide(ContigNode* node){
	if(node1_p == node){
		return 1;
	}
	if(node2_p == node){
		return 2;
	}
	printf("ERROR: tried to get the side of a contig node not adjacent to the contig.\n");
	std::cout << "Node1: " << node1_p << ", Node2: " << node2_p << " Input: " << node << "\n";
	return -1;
}

int Contig::getSide(ContigNode* node, int index){
	if((node1_p == node) && (ind1 == index)){
		return 1;
	}
	if((node2_p == node) && (ind2 == index)){
		return 2;
	}
	printf("ERROR: tried to get the side of a contig node,index pair, but didn't find it on either side.\n");
	std::cout << "Node1: " << node1_p << ", Node2: " << node2_p << " Input: " << node << "\n";
	return -1;
}

void Contig::setSide(int side, ContigNode* node){
	if(side == 1){
		node1_p = node;
	}
	else if(side == 2){
		node2_p = node;
	}
	else printf("ERROR: tried to set side for side other than 1,2.\n");	
}

bool Contig::isIsolated(){
	return node1_p == nullptr && node2_p == nullptr;
}

std::vector<std::pair<Contig*, bool>> Contig::getNeighbors(bool RC){
	if(!RC){ //forward node continuations 
	    if(node2_p){ //if node exists in forward direction 
	    	return node2_p->getFastGNeighbors(ind2);
		}
	}
	else{ //backward node continuations
		if(node1_p){ //if node exists in backward direction
			return node1_p->getFastGNeighbors(ind1);
		}
	}
	return {};
}

bool Contig::checkValidity(){
	if(node1_p){
		if(node1_p->contigs[ind1] != this){
			printf("CONTIG_ERROR: adjacent node at specified index doesn't point back to this contig.\n");
			std::cout << "At " << getFastGName(true) << "\n";			
			return false;
		}
	}
	if(node2_p){
		if(node2_p->contigs[ind2] != this){
			printf("CONTIG_ERROR: adjacent node at specified index doesn't point back to this contig.\n");
			std::cout << "At " << getFastGName(true) << "\n";
			return false;
		}
	}
	return true;
}

string Contig::getFastGName(bool RC){
	stringstream stream;
    stream << "NODE_" << this << "_length_" << seq.length() << "_cov_" << getAvgCoverage();
    if(RC){
    	stream << "'";
    }
    return stream.str();
}

string Contig::getFastGHeader(bool RC){
	stringstream stream;
	stream << ">";
    stream << getFastGName(RC);

    //get neighbors in direction corresponding to RC value
    std::vector<std::pair<Contig*, bool>> neighbors = getNeighbors(RC);

    //if empty return now
    if(neighbors.empty()){
    	stream << ";" ;
    	return stream.str();
    }

    //not empty, add neighbors to line
    stream << ":";
    for(auto it = neighbors.begin(); it != neighbors.end(); it++){
    	Contig* neighbor = it->first;
    	bool RC = it->second;
    	stream << neighbor->getFastGName(RC) << ",";
    }
    string result = stream.str();
    result[result.length()-1] = ';';
    return result;
}

string Contig::getStringRep(){
	stringstream stream;
    stream << seq << "\n";
    stream << node1_p << "," << ind1 << " " << node2_p << "," << ind2 << "\n";
    for(auto it = juncDistances.begin(); it != juncDistances.end(); it++){
        stream << (int)*it << " ";
    }
    stream << "\n";
    stream << coverageSum << "\n";
    return stream.str();
}

Contig::Contig(){
	seq = "";
	node1_p = nullptr;
	node2_p = nullptr;
	ind1 = -1;
	ind2 = -1;
	coverageSum = 0;
	juncDistances = {};
}

Contig::~Contig(){
	juncDistances.clear();
}