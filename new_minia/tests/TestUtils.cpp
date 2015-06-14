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
ReadScanner* scanner;
kmer_type test_kmer;
//Note: nuc order is ACTG = 0123

void setSizeKmer(int k){
    sizeKmer = k;
    //kmer mask is correct!
    if (sizeKmer == (int)(sizeof(kmer_type)*4))
        kmerMask = -1;
    else
        kmerMask=(((kmer_type)1)<<(sizeKmer*2))-(kmer_type)1;
    //test kmer is a reasonable nontrivial sequence.
    test_kmer = ((uint64_t)3091834598457230) & kmerMask;
}

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

void loadBloom(Bloom* fakeBloom, string list[], int numKmers){
    fakeBloom = new Bloom((uint64_t)0, 25);

    std::set<kmer_type> valids;

    kmer_type kmer;
    for(int i = 0; i < numKmers; i++){
        valids.insert(get_canon(getKmerFromString(list[i])));
    }
    fakeBloom->fakify(valids);
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
