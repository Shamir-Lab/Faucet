#include "TraverseReadsTests.h"
#include <stdio.h>
#include <map>
using std::map;

namespace traverseReadsTests {

string fake_read1 = "ACGGGCGAACTTTCATAGGA";
string fake_read2 = "GGCGAACTAGTCCAT";
string fake_read3  = "AACTTTCATACGATT";
Bloom* bloom;

//this is all the kmers from the reads plus two error kmers that cause a 
//TACGA --> ACGATT, ACGAAA branch (fake of length 2)
string valid_5mers[] = {"ACGGG","CGGGC","GGGCG","GGCGA","GCGAA","CGAAC","GAACT"
    ,"AACTT","ACTTT","CTTTC","TTTCA","TTCAT","TCATA","CATAG","ATAGG","TAGGA", "AACTA","ACTAG"
    , "CTAGT", "TAGTC", "AGTCC","GTCCA", "TCCAT" ,"CATAC", "ATACG", "TACGA", "ACGAT","CGATT", "ACGAC", "CGACA"};


void printJunctionMap(){
    printf("Size: %d \n", scanner->getJunctionMap().size());
    for (auto& kv : scanner->getJunctionMap()){
        printf("%s \n", print_kmer(kv.first));
    }
}
void traverseReads(int j){
    scanner->setJ(j);
    scanner->smart_traverse_read(fake_read1);
    scanner->smart_traverse_read(fake_read2);
    scanner->smart_traverse_read(fake_read3);
}

void testTraverseReads_J0(){
    char* testName = (char*)"testTraverseReads_J0";
    scanner = new ReadScanner("mockfile", bloom);

    traverseReads(0);

    if(scanner->getJunctionMap().size() != 4){
        fail(testName, (char*)"junction map size was wrong.");
        printJunctionMap();
        return;
    }
    succeed(testName);
}


void testTraverseReads_J1(){
    char* testName = (char*)"testTraverseReads_J1";
    scanner = new ReadScanner("mockfile", bloom);

    traverseReads(1);

    if(scanner->getJunctionMap().size() != 4){
        fail(testName, (char*)"junction map size was wrong.");
        printJunctionMap();
        return;
    }
    succeed(testName);
}


void testTraverseReads_J2(){
    char* testName = (char*)"testTraverseReads_J2";
    scanner = new ReadScanner("mockfile", bloom);

    traverseReads(2);

    if(scanner->getJunctionMap().size() != 3){
        fail(testName, (char*)"junction map size was wrong.");
        printJunctionMap();
        return;
    }
    succeed(testName);
}

void testTraverseReadTwice_SameJuncs(){
    char* testName = (char*)"testTraverseReadsTwice_SameJuncs";
    scanner = new ReadScanner("mockfile", bloom);

    traverseReads(1);
    traverseReads(1);

    if(scanner->getJunctionMap().size() != 4){
        fail(testName, (char*)"junction map size was wrong.");
        printJunctionMap();
        return;
    }
    succeed(testName);
}

void runTraverseReadsTests(){
    setSizeKmer(5);
    bloom = loadBloom(valid_5mers,30,5);

    testTraverseReads_J0();
    testTraverseReads_J1();
    testTraverseReads_J2();

    testTraverseReadTwice_SameJuncs();
}

}