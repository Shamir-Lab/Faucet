#include "JCheckTests.h"
#include <string>
#include <set>
using std::string;

string fake_read1 = "ACGGGCGAACTTTCATAGGA";
string fake_read2 = "GGCGAACTAGTCCAT";
string fake_read3  = "AACTTTCATACGATT";
Bloom* fakeBloom;
string valid_5mers[] = {"ACGGG","CGGGC","GGGCG","GGCGA","GCGAA","CGAAC","GAACT"
    ,"AACTT","ACTTT","CTTTC","TTTCA","TTCAT","TCATA","CATAG","ATAGG","TAGGA", "AACTA","ACTAG"
    , "CTAGT", "TAGTC", "AGTCC","GTCCA", "TCCAT" ,"CATAC", "ATACG", "TACGA", "ACGAT", "CGATT"};

void jcheck_testForwardJ1NoExtension(){
    char* testName = (char*)"jcheck_testForwardJ1WithNoExtension";
    scanner->setJ(1);

    bool jchecked = scanner->jcheck(getKmerFromString("TCCAT"));

    if(jchecked){
        fail(testName,(char*)"It thought there was an extension.");
        return;
    }
    succeed(testName);
}

void jcheck_testForwardJ1WithExtension(){
    char* testName = (char*)"jcheck_testForwardJ1WithExtension";
    scanner->setJ(1);

    bool jchecked = scanner->jcheck(getKmerFromString("GTCCA"));

    if(!jchecked){
        fail(testName, (char*)"It thought there was no extension.");
        return;
    }
    succeed(testName);
}

void jcheck_testForwardJ2WithExtension(){
    char* testName = (char*)"jcheck_testForwardJ2WithExtension";
    scanner->setJ(2);

    bool jchecked = scanner->jcheck(getKmerFromString("GAACT"));
    
    if(!jchecked){
        fail(testName, (char*)"It thought there was no extension.");
        return;
    }
    succeed(testName);
}

void jcheck_testForwardJ2WithNoExtension(){
    char* testName = (char*)"jcheck_testForwardJ2WithNoExtension";
    scanner->setJ(2);

    bool jchecked = scanner->jcheck(getKmerFromString("ACGAT"));

    if(jchecked){
        fail(testName, (char*)"It thought there was an extension.");
        return;
    }
    succeed(testName);
}

void jcheck_testForwardJ17WithExtension(){
    char* testName = (char*)"jcheck_testForwardJ17WithExtension";
    scanner->setJ(17);

    bool jchecked = scanner->jcheck(getKmerFromString("ACGGG"));

    if(!jchecked){
        fail(testName, (char*)"It thought there was no extension.");
        return;
    }
    succeed(testName);
}

void jcheck_testForwardJ18WithNoExtension(){
    scanner->setJ(18);
    char* testName = (char*)"jcheck_testForwardJ18WithNoExtension";

    bool jchecked = scanner->jcheck(getKmerFromString("ACGGG"));

    if(jchecked){
        fail(testName, (char*)"It thought there was an extension.");
        return;
    }
    succeed(testName);
}



void runJCheckTests(){
    setSizeKmer(5);
    loadBloom(fakeBloom,valid_5mers, 28);
    scanner = new ReadScanner("mockFileName", fakeBloom);

    jcheck_testForwardJ1WithExtension();
    jcheck_testForwardJ1NoExtension();
    jcheck_testForwardJ2WithExtension();
    jcheck_testForwardJ2WithNoExtension();
    jcheck_testForwardJ17WithExtension();
    jcheck_testForwardJ18WithNoExtension();
}

