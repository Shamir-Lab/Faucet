#ifndef BF_SEARCH_RESULT
#define BF_SEARCH_RESULT

//Stores all the info involved in moving from one node to another in the graph by doing a bloom scan
//Also used for scanning junction to junction
struct BfSearchResult{
    BfSearchResult() : kmer(-1){} //use this to tell whether a BfSearchResult has been set or if it was just declared
    BfSearchResult(kmer_type km, bool node, int i, int d, string cont) : kmer(km), isNode(node), index(i), distance(d), contig(cont){} 
    kmer_type kmer;
    bool isNode; //either node/junction or sink
    int index; //index from it that points back to the start point
    int distance; //how far away it was
    string contig; //the contig!
};


#endif