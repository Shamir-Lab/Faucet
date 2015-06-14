#include "FindNextJunctionTests.h"
#include <stdio.h>

string fake_readA = "ACGGGCGAACTTTCATAGGA";
string fake_readB = "GGCGAACTAGTCCAT";
string fake_readC  = "AACTTTCATACGATT";
Bloom* bloomFake;

//this is all the kmers from the reads plus two error kmers that cause a 
//TACGA --> ACGATT, ACGAAA branch (fake of length 2)
string real_5mers[] = {"ACGGG","CGGGC","GGGCG","GGCGA","GCGAA","CGAAC","GAACT"
    ,"AACTT","ACTTT","CTTTC","TTTCA","TTCAT","TCATA","CATAG","ATAGG","TAGGA", "AACTA","ACTAG"
    , "CTAGT", "TAGTC", "AGTCC","GTCCA", "TCCAT" ,"CATAC", "ATACG", "TACGA", "ACGAT","CGATT", "ACGAC", "CGACA"};

void findNextJunction_J1_testFromStart(){
    //int* pos, kmer_type * kmer, string read, int j, Bloom* bloo1
    char* testName = (char*)"findNextJunction_J1_testFromStart";
    int pos = 0;
    kmer_type kmer = getKmerFromString("ACGGG");
    scanner->setJ(1);

    bool foundJunc = scanner->find_next_junction(&pos, &kmer, fake_readA);

    if(!foundJunc){
        fail(testName, (char*)"didn't find a junction.");
        return;
    }
    if(pos != 6){
        fail(testName, (char*)"pos was incorrect.");
        return;
    }    
    if(kmer != getKmerFromString("GAACT")){
        fail(testName, (char*)"kmer was incorrect.");
        return;
    }
    succeed(testName);
}

void findNextJunction_J2_testOffEnd(){
    //int* pos, kmer_type * kmer, string read, int j, Bloom* bloo1
    char* testName = (char*)"findNextJunction_J2_testOffEnd";
    scanner->setJ(2);

    int pos = 13;
    kmer_type kmer = getKmerFromString("CATAG");
    
    bool foundJunc = scanner->find_next_junction(&pos, &kmer, fake_readA);

    if(foundJunc){
        fail(testName, (char*)"found a junction.");
        return;
    }
    succeed(testName);
}

void findNextJunction_J2_Ignore2Branch(){
    //int* pos, kmer_type * kmer, string read, int j, Bloom* bloo1
    char* testName = (char*)"findNextJunction_J2_Ignore2Branch";
    scanner->setJ(2);

    int pos = 6;
    kmer_type kmer = getKmerFromString("CATAC");
    
    bool foundJunc = scanner->find_next_junction(&pos, &kmer, fake_readC);

    if(foundJunc){
        fail(testName, (char*)"found a junction.");
        return;
    }
    succeed(testName);
}

void findNextJunction_J1_Find2Branch(){
    //int* pos, kmer_type * kmer, string read, int j, Bloom* bloo1
    char* testName = (char*)"findNextJunction_J1_Find2Branch";
    scanner->setJ(1);

    int pos = 6;
    kmer_type kmer = getKmerFromString("CATAC");
    
    bool foundJunc = scanner->find_next_junction(&pos, &kmer, fake_readC);

    if(!foundJunc){
        fail(testName, (char*)"didn't find a junction.");
        return;
    }
    if(pos != 8){
        fail(testName, (char*)"pos was incorrect");
        return;
    }
    if(kmer != getKmerFromString("TACGA")){
        fail(testName, (char*)"kmer was incorrect.");
        return;
    }
    succeed(testName);
}

void runFindNextJunctionTests(){
    loadBloom(bloomFake, real_5mers,30);
    scanner = new ReadScanner("mockFile", bloomFake);
    
   findNextJunction_J1_testFromStart();
   findNextJunction_J2_testOffEnd();
   findNextJunction_J2_Ignore2Branch();
   findNextJunction_J1_Find2Branch();
}
