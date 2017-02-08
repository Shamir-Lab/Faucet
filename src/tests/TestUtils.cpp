#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <inttypes.h>
#include <stdint.h>
#include <algorithm> // for max/min
#include <vector> // for sorting_kmers
#include <sys/time.h>
#include <string>
#include "TestUtils.h"
using std::string;

using namespace std;
int64_t nb_reads;
kmer_type test_kmer;
//Note: nuc order is ACTG = 0123

bool kmer_matches_readseq(char* read, kmer_type kmer, int i){
    char* kmerSeq = new char[sizeKmer];
    code2seq(kmer, kmerSeq);
    for(int pos = 0; pos < sizeKmer; pos++){
        if(read[i+pos] != kmerSeq[pos]) return false;
    }
    return true;
}

bool kmer_matches_kmer(kmer_type kmer1, int i1, kmer_type kmer2, int i2){
    char* kmerSeq1 = new char[sizeKmer];
    code2seq(kmer1, kmerSeq1);
    char* kmerSeq2 = new char[sizeKmer];
    code2seq(kmer2, kmerSeq2);
    int length = min(sizeKmer - i1, sizeKmer - i2);
    for(int pos = 0; pos < length; pos++){
        if(kmerSeq1[i1+pos] != kmerSeq2[i2+pos]) return false;
    }
    return true;
}

kmer_type getKmerFromString(string kmerString){
    kmer_type kmer;
    getFirstKmerFromRead(&kmer, &(kmerString[0]));
    return kmer;
}

Bloom* loadBloom(string list[], int numKmers, int k){
    Bloom* fakeBloom = new Bloom((uint64_t)10000, k);

    std::set<kmer_type> valids;

    kmer_type kmer;
    for(int i = 0; i < numKmers; i++){
        valids.insert(getKmerFromString(list[i]));
    }
    fakeBloom->fakify(valids);
    return fakeBloom;
}

void fail(char* testName, char* errorMessage){
    printf("%s: %s \n", testName, errorMessage);
}


void fail(char* testName){
    printf("%s: fail. \n", testName);
}

void succeed(char* testName){
    printf("%s: success! \n", testName);
}
