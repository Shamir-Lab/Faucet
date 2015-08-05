#include <fstream>
#include "ContigNode.h"
using std::ofstream;



ContigNode::ContigNode(Node node){
	for(int i  = 0; i < 5; i++){
		cov[i] = node.cov[i];
		contigs[i] = new Contig();
	}
}

ContigNode::ContigNode(){
	for(int i  = 0; i < 5; i++){
		cov[i] = -1;
		contigs[i] = new Contig();
	}	
}

void ContigNode::update(int nucExt, Contig* contig){
      *contigs[nucExt] = *contig;
}

