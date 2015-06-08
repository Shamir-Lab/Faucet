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

    bool jchecked = jcheck(getKmerFromString("TCCAT"),1,0,fakeBloom);

    if(jchecked){
        fail(testName,(char*)"It thought there was an extension.");
        return;
    }
    succeed(testName);
}

void jcheck_testForwardJ1WithExtension(){
    char* testName = (char*)"jcheck_testForwardJ1WithExtension";

    bool jchecked = jcheck(getKmerFromString("GTCCA"),1,0,fakeBloom);

    if(!jchecked){
        fail(testName, (char*)"It thought there was no extension.");
        return;
    }
    succeed(testName);
}

void jcheck_testForwardJ2WithExtension(){
    char* testName = (char*)"jcheck_testForwardJ2WithExtension";

    bool jchecked = jcheck(getKmerFromString("GAACT"),2,0,fakeBloom);
    
    if(!jchecked){
        fail(testName, (char*)"It thought there was no extension.");
        return;
    }
    succeed(testName);
}

void jcheck_testForwardJ2WithNoExtension(){
    char* testName = (char*)"jcheck_testForwardJ2WithNoExtension";

    bool jchecked = jcheck(getKmerFromString("ACGAT"),2,0,fakeBloom);

    if(jchecked){
        fail(testName, (char*)"It thought there was an extension.");
        return;
    }
    succeed(testName);
}

void jcheck_testForwardJ17WithExtension(){
    char* testName = (char*)"jcheck_testForwardJ17WithExtension";

    bool jchecked = jcheck(getKmerFromString("ACGGG"),17,0,fakeBloom);

    if(!jchecked){
        fail(testName, (char*)"It thought there was no extension.");
        return;
    }
    succeed(testName);
}

void jcheck_testForwardJ18WithNoExtension(){
    char* testName = (char*)"jcheck_testForwardJ18WithNoExtension";

    bool jchecked = jcheck(getKmerFromString("ACGGG"),18,0,fakeBloom);

    if(jchecked){
        fail(testName, (char*)"It thought there was an extension.");
        return;
    }
    succeed(testName);
}



void runJCheckTests(){
    setSizeKmer(5);
    
    loadBloom(fakeBloom,valid_5mers, 28);

    jcheck_testForwardJ1WithExtension();
    jcheck_testForwardJ1NoExtension();
    jcheck_testForwardJ2WithExtension();
    jcheck_testForwardJ2WithNoExtension();
    jcheck_testForwardJ17WithExtension();
    jcheck_testForwardJ18WithNoExtension();
}

