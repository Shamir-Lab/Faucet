#include "GetReadJunctionsTests.h"
#include <stdio.h>
#include <map>
using std::map;

namespace findReadJunctionsTests {

string fake_read1 = "ACGGGCGAACTTTCATAGGA";
string fake_read2 = "GGCGAACTAGTCCAT";
string fake_read3  = "AACTTTCATACGATT";
string fake_read4 = "CGCGCGCGCGC";
Bloom* bloom;
ReadScanner* scanner;

//this is all the kmers from the reads plus two error kmers that cause a 
//TACGA --> ACGATT, ACGAAA branch (fake of length 2)
string valid_5mers[] = {"ACGGG","CGGGC","GGGCG","GGCGA","GCGAA","CGAAC","GAACT"
    ,"AACTT","ACTTT","CTTTC","TTTCA","TTCAT","TCATA","CATAG","ATAGG","TAGGA", "AACTA","ACTAG"
    , "CTAGT", "TAGTC", "AGTCC","GTCCA", "TCCAT" ,"CATAC", "ATACG", "TACGA", "ACGAT","CGATT", "ACGAC", "CGACA"
    , "CGCGC", "GCGCG"};

Junction* getReadJunction(int pos, int dir,  Junction** readJunctions){
    return readJunctions[2*pos + dir];
}

void testFindReadJunctions_NoJunctions(){
    char* testName = (char*)"testFindReadJunctions_NoJunctions";

     scanner->find_read_junctions(fake_read4);

     Junction** readJunctions = scanner->get_read_junctions();

    for(int i = 0; i < fake_read4.length(); i++){
        if(getReadJunction(i,0,readJunctions) || getReadJunction(i,1,readJunctions) ){
            fail(testName);
            return;
        }
    }  
     succeed(testName);
}

void testFindReadJunctions_2Junctions(){
    char* testName = (char*)"testFindReadJunctions_2Junctions";

     scanner->find_read_junctions(fake_read1);//expect junction at 6 and 12

     Junction** readJunctions = scanner->get_read_junctions();

     if(!getReadJunction(6,0, readJunctions)){
        fail(testName, (char*)"Didn't find the forward junction at 6.");
        return;
     } 
     if(!getReadJunction(12,0, readJunctions)){
        fail(testName, (char*)"Didn't find the forward junction at 12");
     }
     succeed(testName);
}

void testFindReadJunctions_backwardJunction(){
    char* testName = (char*)"testFindReadJunctions_BackwardJunction";
    string backwardRead = fake_read1;
    revcomp_sequence(&backwardRead[0], backwardRead.length());

     scanner->find_read_junctions(backwardRead);//expect junction at 6 and 12

     Junction** readJunctions = scanner->get_read_junctions();

     if(!getReadJunction(3,1, readJunctions)){
        fail(testName, (char*)"Didn't find the backward junction at 3.");
        return;
     }  
     if(!getReadJunction(9,1, readJunctions)){
        fail(testName, (char*)"Didn't find the backward junction at 9");
     }
     succeed(testName);
}

void runFindReadJunctionsTests(){
    setSizeKmer(5);
    bloom = loadBloom(valid_5mers,31,5);
    scanner = new ReadScanner("mockfile", bloom, new JChecker(0, bloom));
    
    testFindReadJunctions_NoJunctions();
    testFindReadJunctions_2Junctions();
    testFindReadJunctions_backwardJunction();
}

}