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

char* junctions_filename = new char[100];
char* solids_file = new char[100];
int read_length;
uint64_t genome_size;
bool TwoHash = false;
int64_t nb_reads;
char* kmer_filename = (char*)"solid_27mers_100k";
#include "Bloom.h"
#include "Kmer.h"
#include "ReadScanner.h"
#include "JChecker.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

set<kmer_type> all_kmers;
/*
To run the new version, type make in the directory to compile.

Then, type ./minia 1 2 3 4 5 6 [7] [8], where
1 = name of reads file (current format is each line has a string of characters representing the read)
2 = k
3 = read length
4 = number of distinct kmers.  This will be directly used to size the bloom filter so try to have a good estimate.
5 = false positive rate (.01 is a good base rate) 
6 = j.  j = 0 corresponds to taking direct extensions of the reads, j = 1 is extensions of extensions, etc.
7 = junctions file name.  Will print out to this file. If no filename is given, no file will be printed.
8 = TwoHash.  True if you want to use a bigger bloom filter with only two hash functions.
    Can only use this is you also put in something for 7.
[] indicate optional arguments.

This will load a bloom filter with all the kmers from the reads, then scan through them inserting potential false positives.
No error correction corrently, and no assembly.
*/

inline int handle_arguments(int argc, char *argv[]){
if(argc <  7)
    {
        fprintf (stderr,"usage:\n");
        fprintf (stderr," %s input_file kmer_size min_abundance estimated_nb_reads prefix\n",argv[0]);
        fprintf (stderr,"hints:\n min_abundance ~ 3\n estimated_nb_reads is in bp, does not need to be accurate, only controls memory usage\n prefix is any name you want the results to start with\n");

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
    genome_size = atoll(argv[4]);
    printf("Genome size: %lli .\n", genome_size);
    
    //5th arg: false posiive rate
    fpRate = atof(argv[5]);

    //6th arg: j
    j = atoi(argv[6]);
    printf("j: %d \n", j);

    if(argc > 6){
        strcpy(junctions_filename,argv[7]);
        printf("Junctions file name: %s \n", junctions_filename);
    }

    if(argc > 7){
        TwoHash = atoi(argv[8]);
    }
    if(TwoHash){
        printf("Using 2 hash functions.");
    }
    else{
        printf("Using space-optimal hash settings.");
    }
}
 
int main(int argc, char *argv[])
{
    if (handle_arguments(argc, argv) == 1){
        return 1;
    }

    int estimated_kmers = genome_size;

    Bloom* bloo1;

    if(TwoHash){
        bloo1 = bloo1->create_bloom_filter_2_hash(estimated_kmers, fpRate);
    }
    else{
        bloo1 = bloo1->create_bloom_filter_optimal(estimated_kmers, fpRate);
    }

    bloo1->load_from_reads(solids_file);
    JChecker* jchecker = new JChecker(j, bloo1);
    ReadScanner* scanner = new ReadScanner(solids_file, bloo1, jchecker);
    
    scanner->scanReads();
    scanner->printScanSummary();
    if(argc > 8){
        scanner->getJunctionMap()->writeToFile(junctions_filename);
    }   
    printf("Program reached end. \n");
    return 0;
}