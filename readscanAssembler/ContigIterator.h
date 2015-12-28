#ifndef CONTIG_ITERATOR
#define CONTIG_ITERATOR

#include "ContigGraph.h"
#include "Contig.h"
#include "ContigNode.h"
#include "../utils/Kmer.h"
#include <iterator>

class ContigIterator{ 
private:
    ContigGraph* graph;
    std::unordered_map<kmer_type,ContigNode>::iterator nodeIt;
    int index;
    Contig* findNextContig(); //gets the contig but doesn't increment it and index
    void increment(); //increments index and nodeIt to point to next possible contig
public:
    ContigIterator(ContigGraph* graph);
    Contig* getContig(); //gets contig, moves pointer to next contig or end
    bool hasNextContig();
};

#endif