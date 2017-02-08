#include "FindNextJunctionTests.h"
#include <stdio.h>

namespace findNextJunctionTests
{

string fake_read1 = "ACGGGCGAACTTTCATAGGA";
string fake_read2 = "GGCGAACTAGTCCAT";
string fake_read3  = "AACTTTCATACGATT";
Bloom* bloom;
ReadScanner* scanner;

//this is all the kmers from the reads plus two error kmers that cause a 
//TACGA --> ACGATT, ACGAAA branch (fake of length 2)
string valid_5mers[] = {"ACGGG","CGGGC","GGGCG","GGCGA","GCGAA","CGAAC","GAACT"
    ,"AACTT","ACTTT","CTTTC","TTTCA","TTCAT","TCATA","CATAG","ATAGG","TAGGA", "AACTA","ACTAG"
    , "CTAGT", "TAGTC", "AGTCC","GTCCA", "TCCAT" ,"CATAC", "ATACG", "TACGA", "ACGAT","CGATT", "ACGAC", "CGACA"};

void findNextJunction_J1_testFromStart(){
    //int* pos, kmer_type * kmer, string read, int j, Bloom* bloo1
    char* testName = (char*)"findNextJunction_J1_testFromStart";
    int pos = 0;
    kmer_type kmer = getKmerFromString("ACGGG");
     scanner = new ReadScanner("mockFile", bloom, new JChecker(1, bloom));
    scanner->resetHashes(kmer);

    Junction* junc = scanner->find_next_junction(&pos, &kmer, fake_read1);

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
    int pos = 13;
    kmer_type kmer = getKmerFromString("CATAG");  
    scanner = new ReadScanner("mockFile", bloom, new JChecker(2, bloom));
    scanner->resetHashes(kmer);
    
    Junction* junc = scanner->find_next_junction(&pos, &kmer, fake_read1);

    if(junc){
        fail(testName, (char*)"Returned a junction.");
        return;
    }
    if(pos != 15 ){
        printf("Position %d\n", pos);
        fail(testName, (char*)"Incorrect position.");
        return;
    }
    succeed(testName);
}
void findNextJunction_J1_testAtJunction(){
    //int* pos, kmer_type * kmer, string read, int j, Bloom* bloo1
    char* testName = (char*)"findNextJunction_J1_testAtJunction";
    int pos = 13;
    kmer_type kmer = getKmerFromString("CATAG");  
    scanner = new ReadScanner("mockFile", bloom, new JChecker(2, bloom));
    scanner->getJunctionMap()->createJunction(kmer);
    scanner->resetHashes(kmer);
    
    Junction* junc = scanner->find_next_junction(&pos, &kmer, fake_read1);

    if(pos != 13 ){
        fail(testName, (char*)"Did not return initial position.");
        return;
    }
    succeed(testName);
}

void runFindNextJunctionTests(){
    setSizeKmer(5);
    bloom = loadBloom(valid_5mers,30,5);

   findNextJunction_J1_testFromStart();
   findNextJunction_J2_testOffEnd();
   findNextJunction_J1_testAtJunction();
}

}