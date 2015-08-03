#include <fstream>
#include "Contig.h"
using std::ofstream;


void Contig::setEnds(const ContigNode* n1, int i1, const ContigNode* n2, int i2){
	node1_p = n1;
	ind1 = ind1;
	node2_p = n2;
	ind2 = i2;
}

void Contig::setSeq(const std::string& cont){
	*seq_p += cont;
}

Contig::Contig(){
	seq_p = new std::string;
	*seq_p = "";
	node1_p = nullptr;
	node2_p = nullptr;
	ind1 = -1;
	ind2 = -1;
}

