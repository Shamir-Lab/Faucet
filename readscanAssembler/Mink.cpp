#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <inttypes.h>
#include <stdint.h>
#include <algorithm> // for max/min
#include <vector> // for sorting_kmers
#include <sys/time.h>


using namespace std;

#define NNKS 4 // default minimal abundance for solidity
#define MIN_CONTIG_SIZE (2*sizeKmer+1)

float fpRate = .01;
int j = 1;

string reads_file;
string bloom_input_file;
string junctions_input_file;

int read_length;
uint64_t estimated_kmers;
bool two_hash = false;
bool from_bloom = false;
bool from_junctions = false;
int64_t nb_reads;

#include "../utils/Bloom.h"
#include "../utils/Kmer.h"
#include "ReadScanner.h"
#include "../utils/JChecker.h"
#include "Graph.h"
#include <iostream>
#include <fstream>
#include <string>

set<kmer_type> all_kmers;
string file_prefix;
/*
To run Mink, first type make in the directory to compile.

Then, type ./mink, followed by the following arguments:
-reads_file <>, name of reads file (current format is each line has a string of characters representing the read)
-size_kmer k
-max_read_length <>, upper bound on the size of a read
-estimaed_kmers <>, number of number of distinct kmers.  This will be directly used to size the bloom filter so try to have a good estimate.
-fp <>, false positive rate, default .01
-j j, j = 0 corresponds to taking direct extensions of the reads, j = 1 is extensions of extensions, etc. Default value 1
-file_prefix <>, used for junctions file and contigs file and graph file
--two_hash, if this option is selected a bigger bloom filter with only two hash functions is used
-bloom_file <>, used to shortcut loading the filter if you already have it on file
-junctions_file <>, used to shortcut the readscan if you have access to a junctions file

Note: cannot use junctions_file option without also using bloom_file option

This will load a bloom filter with all the kmers from the reads, then scan through them to find all of the junctions.
It will print all of the junctions to a file named file_prefix.junctions, 
and print some summary info about the run to stdout.
It will then build a graph from the junction map, clean tips from the graph, and print the contigs and the graph to files
file_prefix.contigs, file_prefix.graph.
*/

void argumentError(){
    fprintf (stderr,"usage:\n");
    fprintf (stderr,"./mink -read_file filename -size_kmer k -max_read_length length -estimated_kmers num_kmers -file_prefix prefix");
    fprintf(stderr, "Optional arguments: -fp rate -j j  --two_hash -bloom_file filename -junctions_file filename\n");
}

int handle_arguments(int argc, char *argv[]){
    if(argc <  8)
    {
        fprintf (stderr,"usage:\n");
        fprintf (stderr,"./mink -read_file filename -size_kmer k -max_read_length length -estimated_kmers num_kmers [-fp rate] [-j j] -file_prefix prefix [--two_hash]\n");
        return 1;
    }

    for(int i = 1; i < argc; i++){
        if(0 == strcmp(argv[i] , "-read_file")) //1st arg: read file name
                reads_file = string(argv[i+1]), i++;
        else if(0 == strcmp(argv[i] , "-size_kmer")) // 2rd arg: kmer size.
                setSizeKmer(atoi(argv[i+1])), i++;
        else if(0 == strcmp(argv[i], "-max_read_length")) //3rd arg: read length
                read_length = atoi(argv[i+1]), i++;
        else if(0 == strcmp(argv[i] , "-estimated_kmers")) //estimated number of distinct kmers
                estimated_kmers = atoll(argv[i+1]), i++;
        else if(0 == strcmp(argv[i] , "-fp")) //false posiive rate
                fpRate = atof(argv[i+1]), i++;
        else if(0 == strcmp(argv[i] , "-j")) //value of j for jchecking
                j = atoi(argv[i+1]), i++;
        else if(0 == strcmp(argv[i] , "-file_prefix")) //file prefix used for output files
                file_prefix = string(argv[i+1]), i++;
        else if(0 == strcmp(argv[i] , "--two_hash")) //use two hash function BF instead of optimal size BF
                two_hash = true; 
        else if(0 == strcmp(argv[i] , "-bloom_file")){
                bloom_input_file = string(argv[i+1]);
                from_bloom = true, i++;
        }  
        else if(0 == strcmp(argv[i] , "-junctions_file")){
                junctions_input_file = string(argv[i+1]);
                from_junctions = true, i++;
        }
        else {
                fprintf (stderr, "Cannot parse tag %s\n", argv[i]);
                argumentError();
                return 1; }
        
    }

    if(from_junctions && !from_bloom){
        fprintf(stderr, "Cannot start from junctions without a bloom file.\n");
        argumentError();
        return 1;
    }

    if(from_junctions && from_bloom)
        printf("Starting from after read scan based on bloom and junction files.\n");
    else if(from_bloom)
        printf("Starting from after bloom load based on bloom file.\n");
    else
        printf("Starting at the beginning: will load bloom and find junctions from the read set.");

    std::cout << "Reads file name: " << reads_file << "\n";

    printf("k: %d: \n", sizeKmer);

    printf("Maximal read length: %d\n", read_length);

    printf("Estimated number of distinct kmers, for sizing bloom filter: %lli.\n", estimated_kmers);

    printf("False positive rate: %f\n", fpRate);

    printf("j: %d \n", j);
    
    printf("File prefix: %s\n", &file_prefix[0]);
    
    if(two_hash){
        printf("Using 2 hash functions.\n");
    }
    else{
        printf("Using space-optimal hash settings.\n");
    }

    printf("Size of junction: %d\n", sizeof(Junction));
    printf("Size of node: %d\n", sizeof(Node));
}
 
//create and load bloom filter
Bloom* getBloomFilterFromFile(){
    Bloom* bloom;
        if(two_hash){
            bloom = bloom->create_bloom_filter_2_hash(estimated_kmers, fpRate);
        }
        else{
            bloom = bloom->create_bloom_filter_optimal(estimated_kmers, fpRate);
        }
        bloom->load(&bloom_input_file[0]);
        return bloom;
}

Bloom* getBloomFilterFromReads(){ //handles loading from reads
    Bloom* bloo1;
    Bloom* bloo2;

    if(two_hash){
        bloo1 = bloo1->create_bloom_filter_2_hash(estimated_kmers, fpRate);
        bloo2 = bloo2->create_bloom_filter_2_hash(estimated_kmers, fpRate);
    }
    else{
        bloo1 = bloo1->create_bloom_filter_optimal(estimated_kmers, fpRate);
        bloo2 = bloo2->create_bloom_filter_optimal(estimated_kmers, fpRate);
    }
    load_two_filters(bloo1, bloo2, reads_file);
    delete(bloo1);
    return bloo2;
}

//Builds the junction map from either a file or the readscan
void buildJunctionMapFromReads(JunctionMap* junctionMap, Bloom* bloom, JChecker* jchecker){
     Bloom* junc_bloom;
     ReadScanner* scanner = new ReadScanner(junctionMap, reads_file, bloom, junc_bloom, jchecker);
     
    //scan reads, print summary
    scanner->scanReads();
    scanner->printScanSummary();
}

int main(int argc, char *argv[])
{
    //get all parameters from arguments
    if (handle_arguments(argc, argv) == 1){
        return 1;
    }

    //Build bloom filter from reads or file, and dump to file
    Bloom* bloom;
    if(from_bloom){
        bloom = getBloomFilterFromFile();
    }
    else{
        bloom = getBloomFilterFromReads();
    }
    bloom->dump(&(file_prefix + ".bloom")[0]);

    //create JChecker
    JChecker* jchecker = new JChecker(j, bloom);

    //Build junction map from either reads or file
    JunctionMap* junctionMap = new JunctionMap(bloom, jchecker, read_length);
    if(from_junctions){
        junctionMap->buildFromFile(junctions_input_file);
    }
    else{
        buildJunctionMapFromReads(junctionMap, bloom, jchecker);
    }
    //dump junctions to file
    junctionMap->writeToFile(file_prefix + ".junctions");

    //build graph, dump graph and contigs to file
    Graph* graph = new Graph(bloom, jchecker);
    graph->buildGraph(junctionMap);
    graph->linkNodes();
    graph->cutTips(read_length);
    graph->printContigs(file_prefix + ".contigs");
    graph->printGraph(file_prefix + ".graph");

    //done!
    printf("Program reached end. \n");
    return 0;
}