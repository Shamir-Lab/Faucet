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
int j = 0;

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

Then, type ./mink 1 2 3 4 5 6 7 [8], where
1 = name of reads file (current format is each line has a string of characters representing the read)
2 = k
3 = read length
4 = number of distinct kmers.  This will be directly used to size the bloom filter so try to have a good estimate.
5 = false positive rate (.01 is a good base rate) 
6 = j.  j = 0 corresponds to taking direct extensions of the reads, j = 1 is extensions of extensions, etc.
7 = File prefix. Used for junctions file and contigs file.
8 = TwoHash.  True if you want to use a bigger bloom filter with only two hash functions.
    Can only use this is you also put in something for 7.
[] indicate optional arguments.

This will load a bloom filter with all the kmers from the reads, then scan through them to find all of the junctions.
It will print all of the junctions to a file named file_prefix.junctions, 
and print some summary info about the run to stdout.
*/

inline int handle_arguments(int argc, char *argv[]){
if(argc <  8)
    {
        fprintf (stderr,"usage:\n");
        fprintf (stderr,"./mink reads_file k read_length num_kmers fpRate j file_prefix two_hash\n");
        return 1;
    }

    //1st arg: read file name
    strcpy(solids_file,argv[1]);
    printf("Reads file name: %s \n", solids_file);

    // 2rd arg: kmer size.
    setSizeKmer(atoi(argv[2]));
    printf("k: %d: \n", sizeKmer);

    //3rd arg: read length
    read_length = atoi(argv[3]);

    //4th arg: number of reads
    estimated_kmers = atoll(argv[4]);
    printf("Genome size: %lli .\n", estimated_kmers);
    
    //5th arg: false posiive rate
    fpRate = atof(argv[5]);

    //6th arg: j
    j = atoi(argv[6]);
    printf("j: %d \n", j);

    //7th arg: file prefix
    file_prefix = new string(argv[7]);
    printf("File prefix: %s\n", &((*file_prefix)[0]));

    //8th arg: optional, TwoHash
    if(argc > 8){
        TwoHash = atoi(argv[8]);
    }
    if(TwoHash){
        printf("Using 2 hash functions.\n");
    }
    else{
        printf("Using space-optimal hash settings.\n");
    }

    printf("Size of junction: %d\n", sizeof(Junction));
}
 
int main(int argc, char *argv[])
{
    //get all parameters from arguments
    if (handle_arguments(argc, argv) == 1){
        return 1;
    }

    //create and load bloom filter
    Bloom* bloo1;
    if(TwoHash){
        bloo1 = bloo1->create_bloom_filter_2_hash(estimated_kmers, fpRate);
    }
    else{
        bloo1 = bloo1->create_bloom_filter_optimal(estimated_kmers, fpRate);
    }
    bloo1->load_from_reads(solids_file);

    //create ReadScanner
    JChecker* jchecker = new JChecker(j, bloo1);
    JunctionMap* junctionMap = new JunctionMap(bloo1, jchecker);
    ReadScanner* scanner = new ReadScanner(junctionMap, solids_file, bloo1, jchecker);
    
    //scan reads
    scanner->scanReads();
    scanner->printScanSummary();

    //dump junctions to file
    junctionMap->writeToFile(*file_prefix + ".junctions");

    //build and print graph
    Graph* graph = new Graph(bloo1, jchecker);
    graph->buildGraph(junctionMap);
    graph->linkNodesPrintContigs(*file_prefix + ".contigs");
    graph->printGraph(*file_prefix + ".graph");

    //done!
    printf("Program reached end. \n");
    return 0;
}