#include "JChecker.h"
#include <stdio.h>

//incremental version
//j = 0 always returns true
//j > 0 checks extensions up to j deep from kmer, and returns true if there is a sequence of j extensions
//which returns all positive in the bloom filter
bool JChecker::jcheck(char* kmerSeq, uint64_t nextH0, uint64_t nextH1){
  if(j == 0){
    return true;
  }
  uint64_t workingHash0, workingHash1;
  int lastCount, nextCount;
  lastCount = 1;
  lastHashes[0][0] = nextH0;
  lastHashes[0][1] = nextH1;

  for(int i = 0; i < j; i++){//for each level up to j
    nextCount = 0; //have found no extensions yet
    for(int k = 0; k < lastCount; k++){ //for each kmer in the last level
      workingHash0 = lastHashes[k][0];
      workingHash1 = lastHashes[k][1];
      for(int nt = 0; nt < 4; nt++){ //for each possible extension
        nextHash0 = bloom->roll_hash(workingHash0, NT2int(kmerSeq[i]), nt, 0);
        nextHash1 = bloom->roll_hash(workingHash1, NT2int(kmerSeq[i]), nt, 1);
        if(bloom->contains(nextHash0, nextHash1)){ //add to next level if it's in the bloom filter
          if(i == (j-1)){
            return true;//if this is the last level return true after the first check
          }
          nextHashes[nextCount][0] = nextHash0;
          nextHashes[nextCount][1] = nextHash1;
          nextCount++;
        }
      }
    }

    if(nextCount == 0){ //if there are no kmers in the list now, return false
      return false;
    }
    //reset counts and lists for next level of th search
    lastCount = nextCount;

    tempor = lastHashes;
    lastHashes = nextHashes;
    nextHashes = tempor;
  }
}
    
//Normal version of jchecking, without rolling hash.  
//Old hash! use only for old hash!  For kpomerscanner
bool JChecker::jcheck(kmer_type kmer){
  kmer_type this_kmer, nextKmer;
  int lastCount, nextCount;

  lastCount = 1;
  lastKmers[0] = kmer;

  for(int i = 0; i < j; i++){ //for up to j levels
    nextCount = 0;
    for(int k = 0; k < lastCount; k++){ //for every kmer in the last level
      this_kmer = lastKmers[k];
      for(int nt = 0; nt < 4; nt++){ //for every possible extension
        nextKmer = next_kmer(this_kmer, nt, FORWARD);
        if(bloom->oldContains(get_canon(nextKmer))){//add any positive extensions to the next level
          nextKmers[nextCount] = nextKmer;
          nextCount++;
        }
      }
    }
    if(nextCount == 0){
      return false; //if there are ever no valid kmers in the next level, the kmer does not jcheck
    }
    lastCount = nextCount;
    //switch the pointers to the "last" and "next" arrays so we can use the current "next" one as the next "last" one
    temp = lastKmers;
    lastKmers = nextKmers;
    nextKmers = temp;
  }
  return true;
}

JChecker::JChecker(int jVal, Bloom* bloo){
    j = jVal;
    bloom = bloo;

    //all this is for rolling hash function.. not relevant now
    lastHashes = new uint64_t*[20000];
    nextHashes = new uint64_t*[20000];
    for(int i = 0; i < 20000; i++){
        lastHashes[i] = new uint64_t[2];
        nextHashes[i] = new uint64_t[2];
    }
    lastKmers = new kmer_type[1000];
    nextKmers = new kmer_type[1000];
}