#include <fstream>
#include "Contig.h"
using std::ofstream;


void Contig::setEnds( ContigNode* n1, int i1, ContigNode* n2, int i2){
	node1_p = n1;
	ind1 = i1;
	node2_p = n2;
	ind2 = i2;
}

void Contig::setSeq(const std::string& cont){
	*seq_p += cont;
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

Contig::Contig(){
	seq_p = new std::string("");
	node1_p = nullptr;
	node2_p = nullptr;
	ind1 = -1;
	ind2 = -1;
}

