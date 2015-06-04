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

#include "Kmer.h"
int64_t nb_reads;


void preliminaries(){
    sizeKmer = 32;
    if (sizeKmer == (int)(sizeof(kmer_type)*4))
        kmerMask = -1;
    else
        kmerMask=(((kmer_type)1)<<(sizeKmer*2))-1;
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

void next_kmer_forward_testOverlap(){
    char* testName = (char*)"next_kmer_forward_testOverlap";
    kmer_type kmer = (uint64_t)123235983256209835;
    char* originalString = new char[sizeKmer];
    code2seq(kmer, originalString);
    char* nextString = new char[sizeKmer];
    code2seq(next_kmer(kmer,0,0), nextString);
    bool match = true;
    for(int i = 1; i < sizeKmer; i++){
        if(originalString[i] != nextString[i-1])
            match = false;
    }
    if(!match){
        fail(testName, (char*)"overlapping portion of strings don't match.");
        return;
    }

    succeed(testName);
}

void next_kmer_forward_testLastChar(){
    char* testName = (char*)"ext_kmer_forward_testLastChar";
    kmer_type kmer = (uint64_t)2457368451437897;
    char* nextString = new char[sizeKmer];
    code2seq(next_kmer(kmer,1,0), nextString);
    
    if(nextString[sizeKmer-1] != 'C'){  
        fail(testName);
        return;
    }

    succeed(testName);
}

void next_kmer_backward_testOverlap(){
    char* testName = (char*)"next_kmer_fbackward_testOverlap";
    kmer_type kmer = (uint64_t)123235983256209835;
    char* originalString = new char[sizeKmer];
    code2seq(kmer, originalString);
    char* nextString = new char[sizeKmer];
    code2seq(next_kmer(kmer,2,1), nextString);
    bool match = true;
    for(int i = 0; i < sizeKmer-1; i++){
        if(originalString[i] != nextString[i+1])
            match = false;
    }
    if(!match){
        fail(testName, (char*)"overlapping portion of strings don't match.");
        return;
    }

    succeed(testName);
}

void next_kmer_backward_testFirstChar(){
    char* testName = (char*)"ext_kmer_backward_testFirstChar";
    kmer_type kmer = (uint64_t)2457368451437897;
    char* nextString = new char[sizeKmer];
    code2seq(next_kmer(kmer,1,1), nextString);

    if(nextString[0] != 'C'){  
        fail(testName);
        return;
    }

    succeed(testName);
}

void getFirstKmerFromRead_test(){
    char* testName = (char*) "getFirstKmerFromRead_test";
    char* read = (char*)"ACGGGGGTCAAAATCGGGAATCCGGGGGGAGGCCCTAGT";
    kmer_type kmer;
    getFirstKmerFromRead(&kmer, read);
    char* kmerSeq = new char[sizeKmer];
    code2seq(kmer, kmerSeq);
    for(int i = 0; i < sizeKmer; i++){
        if(kmerSeq[i] != read[i]){
            fail(testName);
            return;
        }
    }

    succeed(testName);
}

void nextKmerInRead_test_forward(){
    char* testName = (char*) "nextKmerInRead_test_forward";
    kmer_type kmer = (uint64_t) 2394867203721228976;
    char* kmerSeq = new char [sizeKmer];
    code2seq(kmer, kmerSeq);

    char* read = new char [100];
    read[0] = 'C';
    strcpy(&read[1], kmerSeq);
    read[1+sizeKmer] = 'T';

    kmer_type next_kmer = next_kmer_in_read(kmer, 1, read, 0);
    char* nextKmerSeq = new char[sizeKmer];
    code2seq(next_kmer, nextKmerSeq);

    for(int i = 0; i < sizeKmer; i++){
        if(nextKmerSeq[i] != read[i+2]){
            fail(testName);
            return;
        }
    }

    succeed(testName);
}


void nextKmerInRead_test_backward(){
    char* testName = (char*) "nextKmerInRead_test_backward";
    kmer_type kmer = (uint64_t) 2394867203721228976;
    char* kmerSeq = new char [sizeKmer];
    code2seq(kmer, kmerSeq);

    char* read = new char [100];
    read[0] = 'C';
    strcpy(&read[1], kmerSeq);
    read[1+sizeKmer] = 'T';

    kmer_type next_kmer = next_kmer_in_read(kmer, 1, read, 1);
    char* nextKmerSeq = new char[sizeKmer];
    code2seq(next_kmer, nextKmerSeq);

    for(int i = 0; i < sizeKmer; i++){
        if(nextKmerSeq[i] != read[i]){
            fail(testName);
            return;
        }
    }

    succeed(testName);
}

int main(int argc, char *argv[]){
    preliminaries();

    next_kmer_forward_testOverlap();
    next_kmer_forward_testLastChar();
    next_kmer_backward_testOverlap();
    next_kmer_backward_testFirstChar();

    getFirstKmerFromRead_test();
    
    nextKmerInRead_test_forward();
    nextKmerInRead_test_backward();


    return 0;
}

