#ifndef CONTIG_GRAPH
#define CONTIG_GRAPH

class ContigGraph; //forward declare

#include <fstream>
#include <iostream> 
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <vector>
#include "../utils/Kmer.h"
#include "../utils/Junction.h"
#include "../utils/Bloom.h"
#include "../utils/JunctionMap.h"
#include "../utils/JChecker.h"
#include "Contig.h"
#include "ContigNode.h"
#include "BfSearchResult.h"
using std::unordered_set;
using std::unordered_map;
using std::string;

class ContigGraph 
{
unordered_map<kmer_type,ContigNode> nodeMap;
std::vector<Contig> isolated_contigs;

public: 
    void addIsolatedContig(Contig contig); 
    ContigNode * getContigNode(kmer_type kmer);//returns a pointer to the ContigNode if it exists or NULL otherwise
    bool isContigNode(kmer_type kmer); //true if a contig node exists for that kmer
    bool isErrorContig(Contig* contig);
    void deleteContig(Contig* contig);

    void deleteErrorContigs();   //remove tips, chimeras, and bubbles. Return number of deleted contigs.
    int collapseDummyNodes(); //removes nodes with only one real extension, merges forward and back contigs
    void printContigs(string filename); //prints the contigs raw
    void printGraph(string fileName); //prints graph : TBD print format- fastg?
    ContigGraph();

    //Creates a contig node if it doesn't already exist
    //If it exists, does nothing and returns the existing one.
    //Otherwise, returns the new one
    ContigNode * createContigNode(kmer_type kmer, Junction junction); 
private:
    unordered_map<kmer_type, ContigNode> contigNodeMap; // maps kmers to ContigNodes after contigs constructed
    void collapseNode(ContigNode * node);
    void cutPath(ContigNode* node, int index); //used on nodes with no backward contig
};

#endif
