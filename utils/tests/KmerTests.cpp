#include "KmerTests.h"
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
#include "TestUtils.h"

using namespace std;

namespace kmerTests{

void shift_kmer_forward_testOverlap(){
    char* testName = (char*)"shift_kmer_forward_testOverlap";
    kmer_type kmer = test_kmer;
    kmer_type originalKmer = kmer;

    shift_kmer(&kmer, 1,0);
    
    if(!kmer_matches_kmer(kmer, 0, originalKmer,1)){
        fail(testName, (char*)"overlapping portion of strings don't match.");
        return;
    }
    succeed(testName);
}

void shift_kmer_forward_testLastChar(){
    char* testName = (char*)"shift_kmer_forward_testLastChar";
    kmer_type kmer = test_kmer;
    char* kmerSeq = new char[sizeKmer];

    shift_kmer(&kmer, 1,0);
    code2seq(kmer,kmerSeq);

    if(kmerSeq[sizeKmer-1] != 'C'){ 
        fail(testName);
        return;
    }
    succeed(testName);
}

void shift_kmer_backward_testOverlap(){
     char* testName = (char*)"shift_kmer_forward_testOverlap";
    kmer_type kmer = test_kmer;
    kmer_type originalKmer = kmer;

    shift_kmer(&kmer, 1,1);
    
    if(!kmer_matches_kmer(kmer, 1, originalKmer,0)){
        fail(testName, (char*)"overlapping portion of strings don't match.");
        return;
    }
    succeed(testName);
}

void shift_kmer_backward_testFirstChar(){
    char* testName = (char*)"shift_kmer_forward_testLastChar";
    kmer_type kmer = test_kmer;
    char* kmerSeq = new char[sizeKmer];
    shift_kmer(&kmer, 2,1);
    code2seq(kmer,kmerSeq);
    
    if(kmerSeq[0] != 'T'){  
        fail(testName);
        return;
    }
    succeed(testName);
}

void next_kmer_forward_testOverlap(){
    char* testName = (char*)"next_kmer_forward_testOverlap";
    kmer_type kmer = test_kmer;
    
    kmer_type nextKmer= next_kmer(kmer,0,0);
    
    if(!kmer_matches_kmer(kmer, 1, nextKmer,0)){
        fail(testName, (char*)"overlapping portion of strings don't match.");
        return;
    }
    succeed(testName);
}

void next_kmer_forward_testLastChar(){
    char* testName = (char*)"ext_kmer_forward_testLastChar";
    kmer_type kmer = test_kmer;
    char* nextString = new char[sizeKmer];
    
    code2seq(next_kmer(kmer,1,0), nextString);
    
    if(nextString[sizeKmer-1] != 'C'){  
        fail(testName);
        return;
    }
    succeed(testName);
}

void next_kmer_backward_testOverlap(){
     char* testName = (char*)"next_kmer_forward_testOverlap";
    kmer_type kmer = test_kmer;
    
    kmer_type nextKmer= next_kmer(kmer,0,1);
    
    if(!kmer_matches_kmer(kmer, 0, nextKmer,1)){
        fail(testName, (char*)"overlapping portion of strings don't match.");
        return;
    }
    succeed(testName);
}

void next_kmer_backward_testFirstChar(){
    char* testName = (char*)"ext_kmer_backward_testFirstChar";
    kmer_type kmer = test_kmer;
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
    
    if(!kmer_matches_readseq(read, kmer, 0)){
        fail(testName);
        return;
    }
    succeed(testName);
}

void nextKmerInRead_test_forward(){
    char* testName = (char*) "nextKmerInRead_test_forward";
    kmer_type kmer = test_kmer;
    char* kmerSeq = new char [sizeKmer];
    code2seq(kmer, kmerSeq);

    char* read = new char [100];
    read[0] = 'C';
    strcpy(&read[1], kmerSeq);
    read[1+sizeKmer] = 'T';

    kmer_type next_kmer = next_kmer_in_read(kmer, 1, read, 0);

    if(!kmer_matches_readseq(read, next_kmer,2)){
        fail(testName);
        return;
    }
    succeed(testName);
}


void nextKmerInRead_test_backward(){
    char* testName = (char*) "nextKmerInRead_test_backward";
    kmer_type kmer = test_kmer;
    char* kmerSeq = new char [sizeKmer];
    code2seq(kmer, kmerSeq);

    char* read = new char [100];
    read[0] = 'C';
    strcpy(&read[1], kmerSeq);
    read[1+sizeKmer] = 'T';

    kmer_type next_kmer = next_kmer_in_read(kmer, 1, read, 1);

    if(!kmer_matches_readseq(read, next_kmer,0)){
        fail(testName);
        return;
    }

    succeed(testName);
}

void advanceKmer_test(){
    char* testName = (char*) "advanceKmer_test";
    kmer_type kmer = test_kmer;
    char* kmerSeq = new char [sizeKmer];
    code2seq(kmer, kmerSeq);
    char* read = new char [100];
    strcpy(&read[0], (char*) "CGGT"); 
    strcpy(&read[4], kmerSeq);
    strcpy(&read[4+sizeKmer], (char*)"ACCCGTTTAAACGTTTAGCCTCTCTGAGAGAAAA");
    
    advance_kmer(&read[0], &kmer, 4,15);

    if(!kmer_matches_readseq(read, kmer, 15)){
        fail(testName);
        return;
    }
    succeed(testName);
}

void runKmerTests(){
    setSizeKmer(27);


    shift_kmer_forward_testOverlap();
    shift_kmer_forward_testLastChar();
    shift_kmer_backward_testOverlap();
    shift_kmer_backward_testFirstChar();

    next_kmer_forward_testOverlap();
    next_kmer_forward_testLastChar();
    next_kmer_backward_testOverlap();
    next_kmer_backward_testFirstChar();

    getFirstKmerFromRead_test();
    
    nextKmerInRead_test_forward();
    nextKmerInRead_test_backward();

    advanceKmer_test();
}

}