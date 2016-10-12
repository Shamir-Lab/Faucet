#include <stdio.h>
#include <map>
#include "gtest/gtest.h"
#include "../ReadScanner.h"
#include "../../utils/Bloom.h"
#include "../../utils/JunctionMap.h"
using std::map;

class readScan : public ::testing::Test {

protected:
    std::vector<string> reads;
    std::vector<string> kmers;

    Bloom* bloom;
    ReadScanner* scanner;
    int j;
    int read_length;
    int estimated_kmers;
    int maxSpacerDist;
    double fpRate;

    JChecker* jchecker;
    JunctionMap* junctionMap;
    Bloom* short_pair_filter;
    Bloom* long_pair_filter;

    // Build a kmer out of a string input
    kmer_type getKmerFromString(string kmerString){
        kmer_type kmer;
        getFirstKmerFromRead(&kmer, &(kmerString[0]));
        return kmer;
    }

    // Add the kmers from a vector of strings to a fake bloom filter
    std::set<kmer_type> addKmers(Bloom* bloom, std::vector<string> kmers) {
        std::set<kmer_type> valids {};
        for (string kmer : kmers) {
            valids.insert(get_canon(getKmerFromString(kmer)));
        }

        bloom->addFakeKmers(valids);
        valids.clear();
    }

    // Create a bloom filter, but make it a fake one
    Bloom* createBloom(){
        Bloom* fakeBloom =  fakeBloom->create_bloom_filter_optimal(estimated_kmers, fpRate);
        fakeBloom->fakify();
        return fakeBloom;
    }

    // This method should be used and modified to print whatever we want ot check about the resulting junction map
    void printJunctionMap(ReadScanner scanner) {
        auto map = scanner.getJunctionMap()->junctionMap;
        printf("Size: %d \n", map.size());
        for (auto& kv : map){
            printf("%s \n", print_kmer(kv.first));
            printf("%d %d %d %d %d \n",
                kv.second.dist[0], kv.second.dist[1], kv.second.dist[2], kv.second.dist[3], kv.second.dist[4]);
        }
    }

    // set up blooms, junction map, jchecker, readscanner for testing
    readScan() {
        j = 0;
        read_length = 30;
        estimated_kmers = 35;
        maxSpacerDist = 8;
        fpRate = .1;
        kmers = {};
        reads = {};

        bloom = createBloom(); 
        jchecker = new JChecker(j, bloom);

        junctionMap = new JunctionMap(bloom, jchecker, read_length);
        string read_scan_file = "mock_file";

        short_pair_filter = short_pair_filter->create_bloom_filter_optimal(estimated_kmers/9, fpRate);
        long_pair_filter = long_pair_filter->create_bloom_filter_optimal(estimated_kmers/6, fpRate);

        scanner = new ReadScanner(junctionMap, read_scan_file, bloom, short_pair_filter, long_pair_filter, jchecker, maxSpacerDist);
        printf("Done initializing!\n");
    }
    ~readScan(){
        reads.clear();
        kmers.clear();
        delete jchecker;
        delete short_pair_filter;
        delete long_pair_filter;
        delete bloom;
        delete junctionMap;
        delete scanner;
    }
};

// This test adds one read, and adds the reads kmers to the bloom filter, scans and prints the junction map 
TEST_F(readScan, singleReadNoJunctions) {
    setSizeKmer(5);

    reads = {"ACGGGCGAACTTTCATAGGA"};
    kmers = {"ACGGG","CGGGC","GGGCG","GGCGA","GCGAA","CGAAC","GAACT","AACTT",
        "ACTTT","CTTTC","TTTCA","TTCAT","TCATA","CATAG","ATAGG","TAGGA"};

    addKmers(bloom, kmers);

    scanner->scanInputRead(reads[0], true);
    std::unordered_map<kmer_type, Junction> map = scanner->getJunctionMap()->junctionMap;
    // Expected junctions & distances
    // TCCTA 
    // 0 0 15 0 1 
    // AACTT 
    // 0 0 15 0 15 
    //assert junction k-mers in map
    ASSERT_EQ(map.count(getKmerFromRead("TCCTA", 0)),1);
    ASSERT_EQ(map.count(getKmerFromRead("AACTT", 0)),1);
    // only these junction in map
    ASSERT_EQ(map.size(),2);

    for (auto& kv : map){
        // assert distances are correct
        if (print_kmer(kv.first)=="TCCTA"){
            ASSERT_EQ(kv.second.dist[2],15);
            ASSERT_EQ(kv.second.dist[4],1);            
        }
        if (print_kmer(kv.first)=="AACTT"){
            ASSERT_EQ(kv.second.dist[2],15);
            ASSERT_EQ(kv.second.dist[4],15);    
        }
    }
}

TEST_F(readScan, singleReadOneFakeJunction) {
    setSizeKmer(5);

    // added k-mers in BF "AACTC", "ACTCC" create fake junction and branch of length 2
    reads = {"ACGGGCGAACTTTCATAGGA"};
    kmers = {"ACGGG","CGGGC","GGGCG","GGCGA","GCGAA","CGAAC","GAACT","AACTT","AACTC","ACTCC",
        "ACTTT","CTTTC","TTTCA","TTCAT","TCATA","CATAG","ATAGG","TAGGA"};

    addKmers(bloom, kmers);

    scanner->scanInputRead(reads[0], true);
    std::unordered_map<kmer_type, Junction> map = scanner->getJunctionMap()->junctionMap;
    // Expected junctions & distances
    // CCTAT 
    // 0 0 0 15 3 
    // GAACT 
    // 0 0 15 0 13 
    //assert junction k-mers in map
    ASSERT_EQ(map.count(getKmerFromRead("CCTAT", 0)),1);
    ASSERT_EQ(map.count(getKmerFromRead("GAACT", 0)),1);
    // only these junction in map
    ASSERT_EQ(map.size(),2);
    for (auto& kv : map){
        // assert distances are correct
        if (print_kmer(kv.first)=="CCTAT"){
            ASSERT_EQ(kv.second.dist[3],15);
            ASSERT_EQ(kv.second.dist[4],3);            
        }
        if (print_kmer(kv.first)=="GAACT"){
            ASSERT_EQ(kv.second.dist[2],15);
            ASSERT_EQ(kv.second.dist[4],13);    
        }
    }
}

// Long read, no junctions
TEST_F(readScan, LongReadNoJunctions) {
    setSizeKmer(5);
    reads = {"ACGGGCGAACTTTCATAGGATCGCACTCAC"};
    kmers = {"ACGGG","CGGGC","GGGCG","GGCGA","GCGAA","CGAAC","GAACT","AACTT",
        "ACTTT","CTTTC","TTTCA","TTCAT","TCATA","CATAG","ATAGG","TAGGA",
        "AGGAT", "GGATC", "GATCG", "ATCGC", "TCGCA", "CGCAC", "GCACT", "GCACT", 
        "CACTC", "ACTCA", "CTCAC"};
        
    addKmers(bloom, kmers);

    scanner->scanInputRead(reads[0], true);
    std::unordered_map<kmer_type, Junction> map = scanner->getJunctionMap()->junctionMap;
    // Expected junctions & distances
    // ATCGC 
    // 1 0 0 0 3 
    // TGCGA 
    // 0 0 1 0 11 
    // CGATC 
    // 0 1 0 0 3 
    // GGATC 
    // 0 0 0 1 12 
    // TTCAT 
    // 12 0 0 0 15 
    // TTCGC 
    // 0 1 0 0 15 
    // GGCGA 
    // 1 0 0 0 7 

    //assert some of junction k-mers in map
    ASSERT_EQ(map.count(getKmerFromRead("ATCGC", 0)),1);
    ASSERT_EQ(map.count(getKmerFromRead("GGATC", 0)),1);
    // only map is correct size
    ASSERT_EQ(map.size(),7);
    for (auto& kv : map){
        // assert distances are correct
        if (print_kmer(kv.first)=="GGATC"){
            ASSERT_EQ(kv.second.dist[3],1);
            ASSERT_EQ(kv.second.dist[4],12);            
        }
        if (print_kmer(kv.first)=="TTCGC"){
            ASSERT_EQ(kv.second.dist[1],1);
            ASSERT_EQ(kv.second.dist[4],15);    
        }
    }
}


// additional tests wanted:
// edge case: read that's a tandem repeat, no junctions
// read with k-mer missing - see getValidReads mechanism works


// Same thing but with three reads
TEST_F(readScan, buildFullMap) {
    setSizeKmer(5);
    reads = {"ACGGGCGAACTTTCATAGGA", "GGCGAACTAGTCCAT", "AACTTTCATACGATT"};
    kmers = {"ACGGG","CGGGC","GGGCG","GGCGA","GCGAA","CGAAC","GAACT","AACTT","ACTTT",
        "CTTTC","TTTCA","TTCAT","TCATA","CATAG","ATAGG","TAGGA","GGCGA", "GCGAA", "CGAAC", 
        "GAACT", "AACTA","ACTAG", "CTAGT", "TAGTC", "AGTCC","GTCCA", "TCCAT","AACTT", 
        "ACTTT", "CTTTC", "TTTCA", "TTCAT", "TCATA", "CATAC", "ATACG", "TACGA", "ACGAT","CGATT"};
    addKmers(bloom, kmers);

    scanner->scanInputRead(reads[0], true);
    scanner->scanInputRead(reads[1], true);
    scanner->scanInputRead(reads[2], true);
    std::unordered_map<kmer_type, Junction> map = scanner->getJunctionMap()->junctionMap;

    // Expected junctions & distances
    // CTAGT 
    // 0 8 3 0 3 
    // TCATA 
    // 0 10 0 6 12 
    // GAACT 
    // 3 0 12 0 13 
    //assert some of junction k-mers in map
    ASSERT_EQ(map.count(getKmerFromRead("CTAGT", 0)),1);
    ASSERT_EQ(map.count(getKmerFromRead("GAACT", 0)),1);
    // only map is correct size
    ASSERT_EQ(map.size(),3);
    for (auto& kv : map){
        // assert distances are correct
        if (print_kmer(kv.first)=="GAACT"){
            ASSERT_EQ(kv.second.dist[0],3);
            ASSERT_EQ(kv.second.dist[4],13);            
        }
        if (print_kmer(kv.first)=="CTAGT"){
            ASSERT_EQ(kv.second.dist[1],8);
            ASSERT_EQ(kv.second.dist[4],3);    
        }
    }
    // printJunctionMap(*scanner);
}

TEST_F(readScan, buildFSmallMap) {
    setSizeKmer(3);
    reads = {"CATTG", "GATTC"};
    kmers = {"CAT", "GAT", "ATT", "TTG" ,"TTC"};
    addKmers(bloom, kmers);

    scanner->scanInputRead(reads[0], true);
    scanner->scanInputRead(reads[1], true);
    std::unordered_map<kmer_type, Junction> map = scanner->getJunctionMap()->junctionMap;

    
    printJunctionMap(*scanner);
}

// add separate to JunctionMapTest
// test building of map, then removal of complex junctions - closer to 
// then

// int main(int ac, char* av[])
// {
//   testing::InitGoogleTest(&ac, av);
//   return RUN_ALL_TESTS();
// }
