#include "BloomTests.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

namespace bloomTests{

Bloom* bloom;

char get_random_nuc(){
    return rand() % 4;
}

kmer_type get_random_kmer(){
    uint64_t rand_kmer = 0;
    for(int i = 0; i < sizeKmer; i++){
        rand_kmer ^= (uint64_t)get_random_nuc();
        rand_kmer <<= 2;
    }
    return rand_kmer;
}

void test_false_positive_rate(uint64_t bloomSize, int sampleSize, float fpRate){
    
    setSizeKmer(25);
    bloom = bloom->create_bloom_filter_optimal(bloomSize, fpRate);
    for(int i = 0; i < bloomSize; i++){
        bloom->add(get_random_kmer());
    }
    int fpCount = 0;
    for(int i = 0; i < sampleSize; i++){
        if(bloom->contains(get_random_kmer())){
            fpCount += 1;
        }
    }
    printf("Size %lli, desired rate %f: %f \n", 
        bloomSize, fpRate, (float)fpCount/(float)sampleSize);
}

void test_speed_raw(uint64_t bloomSize, int sampleSize, float fpRate){
    
    setSizeKmer(25);
    bloom = bloom->create_bloom_filter_optimal(bloomSize, fpRate);
    for(int i = 0; i < bloomSize; i++){
        bloom->add(get_random_kmer());
    }
    time_t start, stop;
    time(&start);
    int fpCount = 0;
    kmer_type kmer = (uint64_t)0;
    for(int i = 0; i < sampleSize; i++){
        bloom->contains(kmer);
        kmer++;
    }
    time(&stop);
    printf("Raw queries per second: %f \n", sampleSize / difftime(stop,start));
}

void test_speed_incremental(uint64_t bloomSize, int sampleSize, float fpRate){
    
    setSizeKmer(25);
    bloom = bloom->create_bloom_filter_optimal(bloomSize, fpRate);
    kmer_type kmer = get_random_kmer();
    uint64_t hash0 = bloom->get_rolling_hash(kmer,0);
    uint64_t hash1 = bloom->get_rolling_hash(kmer,1);
    for(int i = 0; i < bloomSize; i++){
        bloom->add(get_random_kmer());
    }
    time_t start, stop;
    time(&start);
    int fpCount = 0;
    int newNuc = 0;
    int oldNuc = 3;
    for(int i = 0; i < sampleSize; i++){
        bloom->contains(hash0, hash1);
        hash0 = bloom->roll_hash(hash0, oldNuc, newNuc, 0);
        hash1 = bloom->roll_hash(hash1, oldNuc, newNuc, 1);
        newNuc = (newNuc + 1) %4;
        oldNuc = (oldNuc +1)%4;
    }
    time(&stop);
    printf("Incremental Queries per second: %f \n", sampleSize / difftime(stop,start));
}


void test_speed_readscan(uint64_t bloomSize, int sampleSize, float fpRate){
    
    setSizeKmer(25);
    bloom = bloom->create_bloom_filter_optimal(bloomSize, fpRate);
    kmer_type kmer = (uint64_t)0;
    uint64_t hash0 = bloom->get_rolling_hash(kmer,0);
    uint64_t hash1 = bloom->get_rolling_hash(kmer,1);
    for(int i = 0; i < bloomSize; i++){
        bloom->add(get_random_kmer());
    }
    time_t start, stop;
    time(&start);
    int fpCount = 0;
    int newNuc = 0;
    int oldNuc = 3;
    for(int i = 0; i < sampleSize; i++){
        if(i % 100 == 0){
            bloom->contains(kmer);
            kmer++;
        }
        else{
            bloom->contains(hash0, hash1);
            hash0 = bloom->roll_hash(hash0, oldNuc, newNuc, 0);
            hash1 = bloom->roll_hash(hash1, oldNuc, newNuc, 1);
            newNuc = (newNuc + 1) %4;
            oldNuc = (oldNuc +1)%4;
        }
    }
    time(&stop);
    printf("Read Queries per second: %f \n", sampleSize / difftime(stop,start));
}

void test_speed_old(uint64_t bloomSize, int sampleSize, float fpRate){
    
    setSizeKmer(25);
    bloom = bloom->create_bloom_filter_optimal(bloomSize, fpRate);
    for(int i = 0; i < bloomSize; i++){
        bloom->add(get_random_kmer());
    }
    time_t start, stop;
    time(&start);
    int fpCount = 0;
    kmer_type kmer = (uint64_t)0;
    for(int i = 0; i < sampleSize; i++){    
        bloom->oldContains(kmer);
        kmer++;
    }
    time(&stop);
    printf("Old queries per second: %f \n", sampleSize / difftime(stop,start));
}

void runBloomTests(){
    test_false_positive_rate(100000, 100000, .001);
    test_false_positive_rate(100000, 100000, .01);
    test_false_positive_rate(100000, 100000, .1);
    test_speed_raw(100000, 2000000, .01);
    test_speed_incremental(100000, 20000000, .1);
    test_speed_readscan(100000, 20000000, .1);
    test_speed_old(100000, 5000000, .1);
}

}