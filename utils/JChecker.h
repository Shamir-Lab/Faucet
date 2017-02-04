#ifndef JCHECKER 
#define JCHECKER

#include "BloomFilter.hpp"
#include "Bloom.h"
#include "Kmer.h"

class JChecker 
{
    private:
        bloom_filter* bloom;

        //for incremental hashing
        uint64_t ** lastHashes;
        uint64_t ** nextHashes;
        uint64_t ** tempor;
        uint64_t nextHash0, nextHash1;

        //to store the kmers in the BFS
        kmer_type* lastKmers;
        kmer_type* nextKmers;
        kmer_type* temp;

    public:
        int j; //value of j!
        bool jcheck(char* kmerSeq, uint64_t nextH0, uint64_t nextH1);//incremental version
        bool jcheck(kmer_type kmer);//normal version
        JChecker(int jVal, bloom_filter* bloo);
};
#endif