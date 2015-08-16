#include <fstream>
#include <sstream>
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
	result->setSeq(seq.substr(0, seq.length()-sizeKmer) + otherContig->seq);
	std::list<unsigned char> newDistances(juncDistances);
	newDistances.insert(newDistances.end(), otherContig->juncDistances.begin(), otherContig->juncDistances.end());
	result->setJuncDistances(newDistances);
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

	juncDistances.reverse();
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

float Contig::getAvgCoverage(){
	return (float) coverageSum / (float) (juncDistances.size()+1);
}

void Contig::setJuncDistances(std::list<unsigned char> juncDists){
	juncDistances = juncDists;
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

