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


#define NNKS 4 // default minimal abundance for solidity
#define MIN_CONTIG_SIZE (2*sizeKmer+1)

float fpRate = .01;
int j = 1;

char* solids_file = new char[100];

int read_length;
uint64_t estimated_kmers;
bool TwoHash = false;
int64_t nb_reads;
char* kmer_filename = (char*)"solid_27mers_100k";

#include "../utils/Bloom.h"
#include "../utils/Kmer.h"
#include "ReadScanner.h"
#include "../utils/JChecker.h"
#include "Graph.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

set<kmer_type> all_kmers;
string* file_prefix;
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

This will load a bloom filter with all the kmers from the reads, then scan through them to find all of the junctions.
It will print all of the junctions to a file named file_prefix.junctions, 
and print some summary info about the run to stdout.
It will then build a graph from the junction map, clean tips from the graph, and print the contigs and the graph to files
file_prefix.contigs, file_prefix.graph.
*/

inline int handle_arguments(int argc, char *argv[]){
    if(argc <  8)
    {
        fprintf (stderr,"usage:\n");
        fprintf (stderr,"./mink -read_file filename -size_kmer k -max_read_length length -estimated_kmers num_kmers [-fp rate] [-j j] -file_prefix prefix [--two_hash]\n");
        return 1;
    }

    for(int i = 1; i < argc; i++){
        if(0 == strcmp(argv[i] , "-read_file")) //1st arg: read file name
                strcpy(solids_file,argv[i+1]), i++;
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
                file_prefix = new string(argv[i+1]), i++;
        else if(0 == strcmp(argv[i] , "--two_hash")) //use two hash function BF instead of optimal size BF
                TwoHash = true; 
        else {
                fprintf (stderr, "Cannot parse tag %s\n", argv[i]);
                fprintf (stderr,"usage:\n");
                fprintf (stderr,"./mink -read_file filename -size_kmer k -max_read_length length -estimated_kmers num_kmers [-fp rate] [-j j] -file_prefix prefix [--two_hash]\n");
                return 1; }
        
    }
    
    printf("Reads file name: %s \n", solids_file);

    printf("k: %d: \n", sizeKmer);

    printf("Maximal read length: %d\n", read_length);

    printf("Estimated number of distinct kmers, for sizing bloom filter: %lli.\n", estimated_kmers);

    printf("False positive rate: %f\n", fpRate);

    printf("j: %d \n", j);
    
    printf("File prefix: %s\n", &((*file_prefix)[0]));
    
    if(TwoHash){
        printf("Using 2 hash functions.\n");
    }
    else{
        printf("Using space-optimal hash settings.\n");
    }

    printf("Size of junction: %d\n", sizeof(Junction));
    printf("Size of node: %d\n", sizeof(Node));
}
 
int main(int argc, char *argv[])
{
    //get all parameters from arguments
    if (handle_arguments(argc, argv) == 1){
        return 1;
    }

    //create and load bloom filter
    Bloom* bloo1;
    Bloom* bloo2;
    if(TwoHash){
        bloo1 = bloo1->create_bloom_filter_2_hash(estimated_kmers, fpRate);
        bloo2 = bloo2->create_bloom_filter_2_hash(estimated_kmers, fpRate);
    }
    else{
        bloo1 = bloo1->create_bloom_filter_optimal(estimated_kmers, fpRate);
        bloo2 = bloo2->create_bloom_filter_optimal(estimated_kmers, fpRate);
    }

    load_two_filters(bloo1, bloo2, solids_file);

    Bloom* bloom = bloo2;

    //create JChecker, JunctionMap, and ReadScanner
    JChecker* jchecker = new JChecker(j, bloom);
    JunctionMap* junctionMap = new JunctionMap(bloom, jchecker, read_length);
    ReadScanner* scanner = new ReadScanner(junctionMap, solids_file, bloom, jchecker);
    
    //scan reads, print summary
    scanner->scanReads();
    scanner->printScanSummary();

    //dump junctions to file
    junctionMap->writeToFile(*file_prefix + ".junctions");

    //build graph, dump graph and contigs to file
    Graph* graph = new Graph(bloom, jchecker);
    graph->buildGraph(junctionMap);
    graph->linkNodes();
    graph->cutTips(read_length);
    graph->printContigs(*file_prefix + ".contigs");
    graph->printGraph(*file_prefix + ".graph");

    //done!
    printf("Program reached end. \n");
    return 0;
}