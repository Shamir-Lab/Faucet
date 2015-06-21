#ifndef JCHECKER 
#define JCHECKER

#include "Bloom.h"
#include "Kmer.h"

class JChecker 
{
    private:
        int j;
        Bloom* bloom;

        uint64_t ** lastHashes;
        uint64_t ** nextHashes;
        uint64_t ** tempor;
        uint64_t nextHash0, nextHash1;

        kmer_type* lastKmers;
        kmer_type* nextKmers;
        kmer_type* temp;

    public:
        bool jcheck(char* kmerSeq, uint64_t nextH0, uint64_t nextH1);//incremental version
        bool jcheck(kmer_type kmer);//normal version
        JChecker(int jVal, Bloom* bloo);
};
#endif