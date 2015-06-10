#include "TraverseReadsTests.h"
#include <stdio.h>
#include <map>
using std::map;

string fake_readX = "ACGGGCGAACTTTCATAGGA";
string fake_readY = "GGCGAACTAGTCCAT";
string fake_readZ  = "AACTTTCATACGATT";
Bloom* bloom_fake;

//this is all the kmers from the reads plus two error kmers that cause a 
//TACGA --> ACGATT, ACGAAA branch (fake of length 2)
string inserted_5mers[] = {"ACGGG","CGGGC","GGGCG","GGCGA","GCGAA","CGAAC","GAACT"
    ,"AACTT","ACTTT","CTTTC","TTTCA","TTCAT","TCATA","CATAG","ATAGG","TAGGA", "AACTA","ACTAG"
    , "CTAGT", "TAGTC", "AGTCC","GTCCA", "TCCAT" ,"CATAC", "ATACG", "TACGA", "ACGAT","CGATT", "ACGAC", "CGACA"};


void traverseReads(int j){
    smart_traverse_read(fake_readX, bloom_fake, j);
    smart_traverse_read(fake_readY, bloom_fake, j);
    smart_traverse_read(fake_readZ, bloom_fake, j);
}

void testTraverseReads_J0(){
    char* testName = (char*)"testTraverseReads_J0";
    junctionMap = {};

    traverseReads(0);

    if(junctionMap.size() != 4){
        fail(testName, (char*)"junction map size was wrong.");
        printf("Size: %d \n", junctionMap.size());
        for (auto& kv : junctionMap){
            printf("%s \n", print_kmer(kv.first));
        }
        return;
    }
    succeed(testName);
}


void testTraverseReads_J1(){
    char* testName = (char*)"testTraverseReads_J1";
    junctionMap = {};

    traverseReads(1);

    if(junctionMap.size() != 4){
        fail(testName, (char*)"junction map size was wrong.");
        printf("Size: %d \n", junctionMap.size());
        for (auto& kv : junctionMap){
            printf("%s \n", print_kmer(kv.first));
        }
        return;
    }
    succeed(testName);
}


void testTraverseReads_J2(){
    char* testName = (char*)"testTraverseReads_J2";
    junctionMap = {};

    traverseReads(2);

    if(junctionMap.size() != 3){
        fail(testName, (char*)"junction map size was wrong.");
        printf("Size: %d \n", junctionMap.size());
        for (auto& kv : junctionMap){
            printf("%s \n", print_kmer(kv.first));
        }
        return;
    }
    succeed(testName);
}

void testTraverseReadTwice_SameJuncs(){
    char* testName = (char*)"testTraverseReadsTwice_SameJuncs";
    junctionMap = {};

    traverseReads(1);
    traverseReads(1);

    if(junctionMap.size() != 4){
        fail(testName, (char*)"junction map size was wrong.");
        printf("Size: %d \n", junctionMap.size());
        for (auto& kv : junctionMap){
            printf("%s \n", print_kmer(kv.first));
        }
        return;
    }
    succeed(testName);
}

void runTraverseReadsTests(){
    setSizeKmer(5);
    allocateJunctionMap(1000);
    loadBloom(bloom_fake, inserted_5mers,30);

    testTraverseReads_J0();
    testTraverseReads_J1();
    testTraverseReads_J2();

    testTraverseReadTwice_SameJuncs();
}