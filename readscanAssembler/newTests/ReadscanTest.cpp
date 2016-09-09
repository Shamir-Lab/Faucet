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


// Two example test cases
// These aren't engineered to test anything in particular yet, but they give an idea of how to write tests.

// This test adds one read, and adds the reads kmers to the bloom filter, scans and prints the junction map 
TEST_F(readScan, singleReadNoJunctions) {
    reads = {"ACGGGCGAACTTTCATAGGA"};
    kmers = {"ACGGG","CGGGC","GGGCG","GGCGA","GCGAA","CGAAC","GAACT","AACTT",
        "ACTTT","CTTTC","TTTCA","TTCAT","TCATA","CATAG","ATAGG","TAGGA"};

    addKmers(bloom, kmers);

    scanner->scanInputRead(reads[0], true);

    printJunctionMap(*scanner);
}

TEST_F(readScan, singleReadOneFakeJunction) {
    // added k-mers in BF "AACTC", "ACTCC" create fake junction and branch of length 2
    reads = {"ACGGGCGAACTTTCATAGGA"};
    kmers = {"ACGGG","CGGGC","GGGCG","GGCGA","GCGAA","CGAAC","GAACT","AACTT","AACTC","ACTCC",
        "ACTTT","CTTTC","TTTCA","TTCAT","TCATA","CATAG","ATAGG","TAGGA"};

    addKmers(bloom, kmers);

    scanner->scanInputRead(reads[0], true);

    printJunctionMap(*scanner);
}

// Long read, no junctions
TEST_F(readScan, LongReadNoJunctions) {
    reads = {"ACGGGCGAACTTTCATAGGATCGCACTCAC"}; //CCTTAAACGAGAG
    kmers = {"ACGGG","CGGGC","GGGCG","GGCGA","GCGAA","CGAAC","GAACT","AACTT",
        "ACTTT","CTTTC","TTTCA","TTCAT","TCATA","CATAG","ATAGG","TAGGA",
        "AGGAT", "GGATC", "GATCG", "ATCGC", "TCGCA", "CGCAC", "GCACT", "GCACT", 
        "CACTC", "ACTCA", "CTCAC"};
        
    addKmers(bloom, kmers);

    scanner->scanInputRead(reads[0], true);

    printJunctionMap(*scanner);
}

// additional tests wanted:
// edge case: read that's a tandem repeat, no junctions
// read with k-mer missing - see getValidReads mechanism works



// // Same thing but with three reads
// TEST_F(readScan, buildFullMap) {
//     reads = {"ACGGGCGAACTTTCATAGGA", "GGCGAACTAGTCCAT", "AACTTTCATACGATT"};
//     kmers = {"ACGGG","CGGGC","GGGCG","GGCGA","GCGAA","CGAAC","GAACT","AACTT","ACTTT",
//         "CTTTC","TTTCA","TTCAT","TCATA","CATAG","ATAGG","TAGGA","GGCGA", "GCGAA", "CGAAC", 
//         "GAACT", "AACTA","ACTAG", "CTAGT", "TAGTC", "AGTCC","GTCCA", "TCCAT","AACTT", 
//         "ACTTT", "CTTTC", "TTTCA", "TTCAT", "TCATA", "CATAC", "ATACG", "TACGA", "ACGAT","CGATT"};
//     addKmers(bloom, kmers);

//     scanner->scanInputRead(reads[0], true);
//     scanner->scanInputRead(reads[1], true);
//     scanner->scanInputRead(reads[2], true);

//     printJunctionMap(*scanner);
// }

int main(int ac, char* av[])
{
  testing::InitGoogleTest(&ac, av);
  return RUN_ALL_TESTS();
}
