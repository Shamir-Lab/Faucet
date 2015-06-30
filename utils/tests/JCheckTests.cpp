#include "JCheckTests.h"
#include <string>
#include <set>
using std::string;

namespace jCheckTests{

string fake_read1 = "ACGGGCGAACTTTCATAGGA";
string fake_read2 = "GGCGAACTAGTCCAT";
string fake_read3  = "AACTTTCATACGATT";
Bloom* bloom;
string valid_5mers[] = {"ACGGG","CGGGC","GGGCG","GGCGA","GCGAA","CGAAC","GAACT"
    ,"AACTT","ACTTT","CTTTC","TTTCA","TTCAT","TCATA","CATAG","ATAGG","TAGGA", "AACTA","ACTAG"
    , "CTAGT", "TAGTC", "AGTCC","GTCCA", "TCCAT" ,"CATAC", "ATACG", "TACGA", "ACGAT", "CGATT"};
JChecker* jchecker;

bool jcheck(string kmer, int j){
    jchecker = new JChecker(j, bloom);
    uint64_t hash0 =  bloom->get_rolling_hash(getKmerFromString(kmer),0);
    uint64_t hash1 =  bloom->get_rolling_hash(getKmerFromString(kmer),1);
    return jchecker->jcheck(&kmer[0],hash0,hash1);
}

void jcheck_testForwardJ1NoExtension(){
    char* testName = (char*)"jcheck_testForwardJ1WithNoExtension";

    bool jchecked = jcheck("TCCAT", 1);

    if(jchecked){
        fail(testName,(char*)"It thought there was an extension.");
        return;
    }
    succeed(testName);
}

void jcheck_testForwardJ1WithExtension(){
    char* testName = (char*)"jcheck_testForwardJ1WithExtension";

    bool jchecked = jcheck("GTCCA",1);

    if(!jchecked){
        fail(testName, (char*)"It thought there was no extension.");
        return;
    }
    succeed(testName);
}

void jcheck_testForwardJ2WithExtension(){
    char* testName = (char*)"jcheck_testForwardJ2WithExtension";

    bool jchecked = jcheck("GAACT",2);
    
    if(!jchecked){
        fail(testName, (char*)"It thought there was no extension.");
        return;
    }
    succeed(testName);
}

void jcheck_testForwardJ2WithNoExtension(){
    char* testName = (char*)"jcheck_testForwardJ2WithNoExtension";

    bool jchecked = jcheck("ACGAT", 2);

    if(jchecked){
        fail(testName, (char*)"It thought there was an extension.");
        return;
    }
    succeed(testName);
}

void runJCheckTests(){
    setSizeKmer(5);
    bloom = loadBloom(valid_5mers, 28,5);

    jcheck_testForwardJ1WithExtension();
    jcheck_testForwardJ1NoExtension();
    jcheck_testForwardJ2WithExtension();
    jcheck_testForwardJ2WithNoExtension();
}

}