#include <stdio.h>
#include <map>
#include "gtest/gtest.h"
#include "../ReadScanner.h"
#include "../../utils/Bloom.h"
#include "../../utils/JunctionMap.h"
#include "../ContigGraph.h"
using std::map;

class juncMapData : public ::testing::Test {

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
    ContigGraph* contigGraph = new ContigGraph();

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

    void printContigGraph(ContigGraph* graph){
        ContigIterator* contigIt = new ContigIterator(graph);
        //prints contigs that are adjacent to nodes
        while(contigIt->hasNextContig()){
            Contig* contig = contigIt->getContig();
            graph->printContigFastG(&std::cout, contig);
        }
        //prints isolated contigs
        std::vector<Contig> * isolated_contigs = graph->getIsolatedContigs();
        for(auto it = isolated_contigs->begin(); it != isolated_contigs->end(); ++it){
            Contig* contig = &*it;       
            graph->printContigFastG(&std::cout, contig);
        }
    }

    // set up blooms, junction map, jchecker, readscanner for testing
    juncMapData() {
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
        contigGraph = new ContigGraph();
    }
    ~juncMapData(){
        reads.clear();
        kmers.clear();
        delete jchecker;
        delete short_pair_filter;
        delete long_pair_filter;
        delete bloom;
        delete junctionMap;
        delete scanner;
        delete contigGraph;
    }
};


// build junction map of three reads
// TEST_F(juncMapData, buildBranchingPaths) {
//     setSizeKmer(5);
        // j = 0;
//     reads = {"ACGGGCGAACTTTCATAGGA", "GGCGAACTAGTCCAT", "AACTTTCATACGATT"};
//     kmers = {"ACGGG","CGGGC","GGGCG","GGCGA","GCGAA","CGAAC","GAACT","AACTT","ACTTT",
//         "CTTTC","TTTCA","TTCAT","TCATA","CATAG","ATAGG","TAGGA","GGCGA", "GCGAA", "CGAAC", 
//         "GAACT", "AACTA","ACTAG", "CTAGT", "TAGTC", "AGTCC","GTCCA", "TCCAT","AACTT", 
//         "ACTTT", "CTTTC", "TTTCA", "TTCAT", "TCATA", "CATAC", "ATACG", "TACGA", "ACGAT","CGATT"};
//     addKmers(bloom, kmers);

//     scanner->scanInputRead(reads[0], true);
//     scanner->scanInputRead(reads[1], true);
//     scanner->scanInputRead(reads[2], true);
//     std::unordered_map<kmer_type, Junction> map = scanner->getJunctionMap()->junctionMap;

//     // Expected junctions & distances before changes
//     // CTAGT 
//     // 0 8 3 0 3 
//     // TCATA 
//     // 0 10 0 6 12 
//     // GAACT 
//     // 3 0 12 0 13 
//     std::cout << "map before changes\n";
//     printJunctionMap(*scanner);
//     junctionMap->buildBranchingPaths(contigGraph);
//     std::cout << "built branching paths\n";
//     printContigGraph(contigGraph);

//     printf("Destroying complex junctions.\n");
//     junctionMap->destroyComplexJunctions();
//     std::cout << "map before changes\n";
//     printJunctionMap(*scanner);
    
//     printf("Building linear regions.\n");
//     junctionMap->buildLinearRegions(contigGraph);
//     printContigGraph(contigGraph);

//     printf("Checking graph.\n");
//     contigGraph->checkGraph();
// }

TEST_F(juncMapData, smallDoubleJuncMap) {
    setSizeKmer(7);
    j = 0;

    reads = {"AAAAACAGCGATTC", "AAAAAGAGCGATTTA"};
    kmers = {"AAAAACA", "AAAAAGA", "AAAACAG", "AAAAGAG", "AAACAGC", "AAAGAGC",
        "AACAGCG", "AAGAGCG", "ACAGCGA","AGAGCGA","CAGCGAT", "GAGCGAT", "AGCGATT", 
        "GCGATTT" ,"GCGATTC", "CGATTTA"};
    addKmers(bloom, kmers);

    scanner->scanInputRead(reads[0], true);
    scanner->scanInputRead(reads[1], true);
    std::unordered_map<kmer_type, Junction> map = scanner->getJunctionMap()->junctionMap;

    printJunctionMap(*scanner);

    std::cout << "map before changes\n";
    printJunctionMap(*scanner);
    junctionMap->buildBranchingPaths(contigGraph);
    std::cout << "built branching paths\n";
    printContigGraph(contigGraph);

    printf("Destroying complex junctions.\n");
    junctionMap->destroyComplexJunctions();
    std::cout << "map before changes\n";
    printJunctionMap(*scanner);
    
    printf("Building linear regions.\n");
    junctionMap->buildLinearRegions(contigGraph);
    printContigGraph(contigGraph);

    printf("Checking graph.\n");
    contigGraph->checkGraph();

}

TEST_F(juncMapData, endJuncMap) {
    setSizeKmer(7);
    j = 1;

    reads = {"AAAAAACAGCGATTC", "AAAAAACTAAAAAA"}; // single read, first kmer is junction, should poinnt back one
    kmers = {"AAAAAAC", "AAAAACA", "AAAAACT", "AAAACAG", "AAACAGC", "AACAGCG", "ACAGCGA","CAGCGAT", "AGCGATT", "GCGATTC",
        "AAAACTA", "AAACTAA", "AACTAAA", "ACTAAAA", "CTAAAAA", "CTAAAAAA"};
    addKmers(bloom, kmers);

    scanner->scanInputRead(reads[0], true);
    scanner->scanInputRead(reads[1], true);

    std::unordered_map<kmer_type, Junction> map = scanner->getJunctionMap()->junctionMap;

    printJunctionMap(*scanner);

    std::cout << "map before changes\n";
    printJunctionMap(*scanner);
    junctionMap->buildBranchingPaths(contigGraph);
    std::cout << "built branching paths\n";
    printContigGraph(contigGraph);

    printf("Destroying complex junctions.\n");
    junctionMap->destroyComplexJunctions();
    std::cout << "map before changes\n";
    printJunctionMap(*scanner);
    
    printf("Building linear regions.\n");
    junctionMap->buildLinearRegions(contigGraph);
    printContigGraph(contigGraph);

    printf("Checking graph.\n");
    contigGraph->checkGraph();

}

