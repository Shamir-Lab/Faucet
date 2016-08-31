#include <stdio.h>
#include <map>
#include "gtest/gtest.h"
#include "../ReadScanner.h"
#include "../../utils/Bloom.h"
#include "../../utils/JunctionMap.h"
using std::map;

class readscan : public ::testing::Test {

protected:
    string fake_read1 = "ACGGGCGAACTTTCATAGGA";
    string fake_read2 = "GGCGAACTAGTCCAT";
    string fake_read3 = "AACTTTCATACGATT";
    Bloom* bloom;
    ReadScanner* scanner;
    int j = 0;
    int read_length = 30;
    int estimated_kmers = 35;
    int maxSpacerDist = 8;
    double fpRate = .1;

    //this is all the kmers from the reads plus two error kmers that cause a 
    //TACGA --> ACGATT, ACGAAA branch (fake of length 2)
    std::vector<string> kmers_1 {"ACGGG","CGGGC","GGGCG","GGCGA","GCGAA","CGAAC","GAACT"
        ,"AACTT","ACTTT","CTTTC","TTTCA","TTCAT","TCATA","CATAG","ATAGG","TAGGA"};
    std::vector<string> kmers_2 {"GGCGA", "GCGAA", "CGAAC", "GAACT", "AACTA","ACTAG"
        , "CTAGT", "TAGTC", "AGTCC","GTCCA", "TCCAT"};
    std::vector<string> kmers_3 {"AACTT", "ACTTT", "CTTTC", "TTTCA", "TTCAT", "TCATA", "CATAC", "ATACG", "TACGA", "ACGAT","CGATT"};

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
    readscan() {
        setSizeKmer(5);

        bloom = createBloom(); 
        jchecker = new JChecker(j, bloom);

        junctionMap = new JunctionMap(bloom, jchecker, read_length);
        string read_scan_file = "mock_file";

        short_pair_filter = short_pair_filter->create_bloom_filter_optimal(estimated_kmers/9, fpRate);
        long_pair_filter = long_pair_filter->create_bloom_filter_optimal(estimated_kmers/6, fpRate);

        scanner = new ReadScanner(junctionMap, read_scan_file, bloom, short_pair_filter, long_pair_filter, jchecker, maxSpacerDist);
        printf("Done initializing!\n");
    }
};

/**
* Right now only one of these tests can work- if you comment either out, the other will run!
* I suspect the solution is something to do with the test construction and destruction- right now we're not using a destructor, 
* which could be the problem.
*/

// Two example test cases
// These aren't engineered to test anything in particular yet, but they give an idea of how to write tests.

// This test adds one read, and adds the reads kmers to the bloom filter, scans and prints the junction map 
TEST_F(readscan, scan_one_read) {
    addKmers(bloom, kmers_1);

    scanner->scanInputRead(fake_read1, true);

    printJunctionMap(*scanner);
}

// Same thing but with three reads
TEST_F(readscan, build_full_map) {
    addKmers(bloom, kmers_1);
    addKmers(bloom, kmers_2);
    addKmers(bloom, kmers_3);

    scanner->scanInputRead(fake_read1, true);
    scanner->scanInputRead(fake_read2, true);
    scanner->scanInputRead(fake_read3, true);

    printJunctionMap(*scanner);
}

