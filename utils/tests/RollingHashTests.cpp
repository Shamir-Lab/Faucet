#include "RollingHashTests.h"
#include <stdio.h>

namespace rollingHashTests
{
    
Bloom* bloom;
int kVal = 27;

void rotate_right_test_moveBitRight(){
    char* testName = (char*)"rotate_right_test_move1BitRight";
    bloom = new Bloom((uint64_t)1000, kVal); //hashSize should be 10
    uint64_t hash = (uint64_t)(1 << 7);

    uint64_t rotated = bloom->rotate_right(hash, 7);

    if(rotated != 1){
        fail(testName);
    }
    succeed(testName);
}


void rotate_right_test_wrapBit(){
    char* testName = (char*)"rotate_right_test_wrapBit";
    bloom = new Bloom((uint64_t)1000, kVal); //hashSize should be 10
    uint64_t hash = (uint64_t)1;

    uint64_t rotated = bloom->rotate_right(hash, 7);

    if(rotated != (1 << 3)){
        fail(testName);
    }
    succeed(testName);
}



void rotate_right_test_noOverFlow(){
    char* testName = (char*)"rotate_right_test_noOverFlow";
    bloom = new Bloom((uint64_t)1000, kVal); //hashSize should be 10
    uint64_t hash = (uint64_t)1000;

    uint64_t rotated = bloom->rotate_right(hash, 7);

    if(rotated > 1024){
        fail(testName);
    }
    succeed(testName);
}


void rotate_left_test_moveBitLeft(){
    char* testName = (char*)"rotate_Left_test_move1BitLeft";
    bloom = new Bloom((uint64_t)1000, kVal); //hashSize should be 10
    uint64_t hash = (uint64_t)1;

    uint64_t rotated = bloom->rotate_left(hash, 7);

    if(rotated != (1 << 7)){
        fail(testName);
    }
    succeed(testName);
}

void rotate_left_test_wrapBit(){
    char* testName = (char*)"rotate_left_test_wrapBit";
    bloom = new Bloom((uint64_t)1000, kVal); //hashSize should be 10
    uint64_t hash = (uint64_t)(1 << 5);

    uint64_t rotated = bloom->rotate_left(hash, 7);

    if(rotated != (1 << 2)){
        fail(testName);
    }
    succeed(testName);
}


void rotate_left_test_noOverFlow(){
    char* testName = (char*)"rotate_left_test_noOverFlow";
    bloom = new Bloom((uint64_t)1000, kVal); //hashSize should be 10
    uint64_t hash = (uint64_t)1000;

    uint64_t rotated = bloom->rotate_left(hash, 7);

    if(rotated > 1024){
        fail(testName);
    }
    succeed(testName);
}

void roll_hash_hash_func0_bigbloom_checkSame(){
    bloom = new Bloom((uint64_t)1000,kVal);
    char* testName = (char*)"roll_hash_hash_func0_bigbloom_checkSame";
    char* kmerSeq = (char*)"ACTTACTGGGCTCTATTGCGTATCGATCGATCGATGCATCTACCCCCATCTAATTAGAGTGAATAGATCGATCGATCGCATACTCAGCATAGCTATA";
    kmer_type firstKmer, secondKmer, kmer;

    getFirstKmerFromRead(&firstKmer, &(kmerSeq[0]));
    getFirstKmerFromRead(&secondKmer, &(kmerSeq[sizeKmer]));

    uint64_t rolledHash =  bloom->get_rolling_hash(firstKmer, 0);
    uint64_t calculatedHash = bloom->get_rolling_hash(secondKmer, 0);

    for(int i = 0; i < sizeKmer; i++){
        rolledHash = bloom->roll_hash(rolledHash, NT2int(kmerSeq[i]), NT2int(kmerSeq[i+sizeKmer]),0);
    }

    if(rolledHash != calculatedHash){
        fail(testName);
        return;
    }
    succeed(testName);
}

void roll_hash_hash_func1_bigbloom_checkSame(){
    bloom = new Bloom((uint64_t)10000, kVal);
    char* testName = (char*)"roll_hash_hash_func1_bigbloom_checkSame";
    char* kmerSeq = (char*)"ACTTACTGGGCTCTATTGCGTATCGATCGATCGATGCATCTACCCCCATCTAATTAGAGTGAATAGATCGATCGATCGCATACTCAGCATAGCTATA";
    kmer_type firstKmer, secondKmer, kmer;

    getFirstKmerFromRead(&firstKmer, &(kmerSeq[0]));
    getFirstKmerFromRead(&secondKmer, &(kmerSeq[sizeKmer]));

    uint64_t rolledHash =  bloom->get_rolling_hash(firstKmer, 1);
    uint64_t calculatedHash = bloom->get_rolling_hash(secondKmer, 1);

   
    for(int i = 0; i < sizeKmer; i++){
         rolledHash = bloom->roll_hash(rolledHash, NT2int(kmerSeq[i]), NT2int(kmerSeq[i+sizeKmer]),1);
    }

    if(rolledHash != calculatedHash){
        fail(testName);
        return;
    }
    succeed(testName);
}

void roll_hash_hash_func0_smallbloom_checkSame(){
    bloom = new Bloom((uint64_t)100000, kVal);
    char* testName = (char*)"roll_hash_hash_func0_smallbloom_checkSame";
    char* kmerSeq = (char*)"ACTTACTGGGCTCTATTGCGTATCGATCGATCGATGCATCTACCCCCATCTAATTAGAGTGAATAGATCGATCGATCGCATACTCAGCATAGCTATA";
    kmer_type firstKmer, secondKmer, kmer;

    getFirstKmerFromRead(&firstKmer, &(kmerSeq[0]));
    getFirstKmerFromRead(&secondKmer, &(kmerSeq[sizeKmer]));

    uint64_t rolledHash =  bloom->get_rolling_hash(firstKmer, 0);
    uint64_t calculatedHash = bloom->get_rolling_hash(secondKmer, 0);

    for(int i = 0; i < sizeKmer; i++){
         rolledHash = bloom->roll_hash(rolledHash, NT2int(kmerSeq[i]), NT2int(kmerSeq[i+sizeKmer]),0);
    }

    if(rolledHash != calculatedHash){
        fail(testName);
        return;
    }
    succeed(testName);
}

void roll_hash_hash_func1_smallbloom_checkSame(){
    bloom = new Bloom((uint64_t)1000000, kVal);
    char* testName = (char*)"roll_hash_hash_func1_smallbloom_checkSame";
    char* kmerSeq = (char*)"ACTTACTGGGCTCTATTGCGTATCGATCGATCGATGCATCTACCCCCATCTAATTAGAGTGAATAGATCGATCGATCGCATACTCAGCATAGCTATA";
    kmer_type firstKmer, secondKmer, kmer;

    getFirstKmerFromRead(&firstKmer, &(kmerSeq[0]));
    getFirstKmerFromRead(&secondKmer, &(kmerSeq[sizeKmer]));

    uint64_t rolledHash =  bloom->get_rolling_hash(firstKmer, 1);
    uint64_t calculatedHash = bloom->get_rolling_hash(secondKmer, 1);

   
    for(int i = 0; i < sizeKmer; i++){
         rolledHash = bloom->roll_hash(rolledHash, NT2int(kmerSeq[i]), NT2int(kmerSeq[i+sizeKmer]),1);
    }

    if(rolledHash != calculatedHash){
        fail(testName);
        return;
    }
    succeed(testName);
}

void advance_hash_test_checkSame(){
    bloom = new Bloom((uint64_t)10000, kVal);
    char* testName = (char*) "advance_hash_test_checkSame";
    char* kmerSeq = (char*)"ACTTACTGGGCTCTATTGCGTATCGATCGATCGATGCATCTACCCCCATCTAATTAGAGTGAATAGATCGATCGATCGCATACTCAGCATAGCTATA";
    kmer_type firstKmer, secondKmer;

    getFirstKmerFromRead(&firstKmer, &(kmerSeq[0]));
    getFirstKmerFromRead(&secondKmer, &(kmerSeq[kVal]));

    uint64_t advancedHash0 =  bloom->get_rolling_hash(firstKmer, 0);
    uint64_t advancedHash1 =  bloom->get_rolling_hash(firstKmer, 1);
    uint64_t calculatedHash0 = bloom->get_rolling_hash(secondKmer, 0);
    uint64_t calculatedHash1 = bloom->get_rolling_hash(secondKmer, 1);

    bloom->advance_hash(&kmerSeq[0], &advancedHash0, &advancedHash1,0,kVal);

    if(advancedHash0 != calculatedHash0){
        fail(testName, (char*)"hash 0 doesn't match.");
        return;
    }
    if(advancedHash1 != calculatedHash1){
        fail(testName, (char*)"hash 1 doesn't match.");
        return;
    }
    succeed(testName);

}


void runRollingHashTests(){
    setSizeKmer(kVal);

    rotate_right_test_moveBitRight();
    rotate_right_test_wrapBit();
    rotate_right_test_noOverFlow();

    rotate_left_test_moveBitLeft();
    rotate_left_test_wrapBit();
    rotate_left_test_noOverFlow();

    roll_hash_hash_func0_bigbloom_checkSame();
     roll_hash_hash_func1_bigbloom_checkSame();
     roll_hash_hash_func0_smallbloom_checkSame();
     roll_hash_hash_func1_smallbloom_checkSame();

    advance_hash_test_checkSame();  
}

}