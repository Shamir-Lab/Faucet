#include <stdio.h>
#include <map>
#include "gtest/gtest.h"
#include "../ReadScanner.h"
#include "../../utils/Bloom.h"
#include "../../utils/JunctionMap.h"
using std::map;

string fake_read1 = "ACGGGCGAACTTTCATAGGA";
string fake_read2 = "GGCGAACTAGTCCAT";
string fake_read3  = "AACTTTCATACGATT";
Bloom* bloom;
ReadScanner* scanner;
int j = 0;
int read_length = 30;
int estimated_kmers = 35;
int maxSpacerDist = 8;
double fpRate = .1;

//this is all the kmers from the reads plus two error kmers that cause a 
//TACGA --> ACGATT, ACGAAA branch (fake of length 2)
string valid_5mers[] = {"ACGGG","CGGGC","GGGCG","GGCGA","GCGAA","CGAAC","GAACT"
    ,"AACTT","ACTTT","CTTTC","TTTCA","TTCAT","TCATA","CATAG","ATAGG","TAGGA", "AACTA","ACTAG"
    , "CTAGT", "TAGTC", "AGTCC","GTCCA", "TCCAT" ,"CATAC", "ATACG", "TACGA", "ACGAT","CGATT", "ACGAC", "CGACA"};

kmer_type getKmerFromString(string kmerString){
    kmer_type kmer;
    getFirstKmerFromRead(&kmer, &(kmerString[0]));
    return kmer;
}

Bloom* loadBloom(string list[], int numKmers, int k){
    Bloom* fakeBloom = new Bloom((uint64_t)10000, k);

    std::set<kmer_type> valids;

    kmer_type kmer;
    for(int i = 0; i < numKmers; i++){
        valids.insert(getKmerFromString(list[i]));
    }
    fakeBloom->fakify(valids);
    return fakeBloom;
}

void printJunctionMap(ReadScanner scanner) {
    auto map = scanner.getJunctionMap()->junctionMap;
    printf("Size: %d \n", map.size());
    for (auto& kv : map){
        printf("%s \n", print_kmer(kv.first));
    }
}

TEST(readscan_test, build_map) {
    setSizeKmer(5);

    JChecker* jchecker = new JChecker(j, bloom);
    JunctionMap* junctionMap = new JunctionMap(bloom, jchecker, read_length);

    string read_scan_file = "mock_file";

    bloom = loadBloom(valid_5mers,30,5);

    Bloom* short_pair_filter = short_pair_filter->create_bloom_filter_optimal(estimated_kmers/9, fpRate);
    Bloom* long_pair_filter = long_pair_filter->create_bloom_filter_optimal(estimated_kmers/6, fpRate);

    scanner = new ReadScanner(junctionMap, read_scan_file, bloom, short_pair_filter, long_pair_filter, jchecker, maxSpacerDist);

    printf("Scanning reads!\n");
    scanner->scanInputRead(fake_read1, true);
    scanner->scanInputRead(fake_read2, true);
    scanner->scanInputRead(fake_read3, true);

    printJunctionMap(*scanner);
}

