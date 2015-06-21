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
uint64_t genome_size;bool TwoHash = false;
int64_t nb_reads;
char* kmer_filename = (char*)"solid_27mers_100k";
#include "../utils/Bloom.h"
#include "../utils/Kmer.h"
#include "KpomerScanner.h"
#include "../utils/JChecker.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

set<kmer_type> all_kmers;
/*
To run the new version, type make in the directory to compile.

Then, type ./minia 1 2 3 4 5 6 [7] [8], where
1 = name of kmers file (current format is each line has a string of characters representing the kmer)
2 = k
3 = number of distinct kmers.  This will be used to size the bloom filter so try to have a good estimate.
4 = false positive rate (.01 is a good base rate) 
5 = j.  j = 0 corresponds to taking direct extensions of the reads, j = 1 is extensions of extensions, etc.
6 = junctions file name.  Will print out to this file. If no filename is given, no file will be printed.
7 = TwoHash.  True if you want to use a bigger bloom filter with only two hash functions.
    Can only use this is you also put in something for 7.
[] indicate optional arguments.

This will load a bloom filter with all the kmers from the reads, then scan through them inserting potential false positives.
No error correction corrently, and no assembly.
*/

inline int handle_arguments(int argc, char *argv[]){
if(argc <  6)
    {
        fprintf (stderr,"usage:\n");
        fprintf (stderr," %s input_file kmer_size min_abundance estimated_nb_reads prefix\n",argv[0]);
        fprintf (stderr,"hints:\n min_abundance ~ 3\n estimated_nb_reads is in bp, does not need to be accurate, only controls memory usage\n prefix is any name you want the results to start with\n");

        return 1;
    }

    //1st arg: read file name
    strcpy(solids_file,argv[1]);
    printf("Kmers file name: %s \n", solids_file);

    // 2rd arg: kmer size.
    setSizeKmer(atoi(argv[2]));
    printf("k: %d: \n", sizeKmer);

    //4th arg: number of reads
    genome_size = atoll(argv[3]);
    printf("Genome size: %lli .\n", genome_size);
    
    //5th arg: false posiive rate
    fpRate = atof(argv[4]);
    printf("Fp rate: %f \n", fpRate);

    //6th arg: j
    j = atoi(argv[5]);
    printf("j: %d \n", j);

    if(argc > 6){
        strcpy(junctions_filename,argv[6]);
        printf("Junctions file name: %s \n", junctions_filename);
    }

    if(argc > 7){
        TwoHash = atoi(argv[7]);
    }
    if(TwoHash){
        printf("Using 2 hash functions.\n");
    }
    else{
        printf("Using space-optimal hash settings.\n");
    }
}

void write_kmers(const char* reads_filename){
    ifstream solidReads;
    solidReads.open(reads_filename);
    kmer_type kmer;
    char* kmerSeq = new char[sizeKmer];
    int readsProcessed = 0;
    printf("Building kmer set.\n");
    string read;
    while (getline(solidReads, read))
    {
        getFirstKmerFromRead(&kmer,&read[0]);

        for (int i = 0; i <= read.length() - sizeKmer ; i++, 
            shift_kmer(&kmer, NT2int(read[i+sizeKmer-1]), 0)){
            all_kmers.insert(kmer);
            readsProcessed++;
            if ((readsProcessed%10000)==0) fprintf (stderr,"%c %lld",13,(long long)readsProcessed);
        }
    }
    solidReads.close();
    printf("Done building kmer set.\n");
    ofstream solidKmers;
    solidKmers.open(kmer_filename);
    
    printf("Writing to kmer file\n");
    set<kmer_type>::iterator it;
    for(it = all_kmers.begin(); it != all_kmers.end(); it++){
        code2seq(*it, kmerSeq);
        solidKmers << kmerSeq;
        solidKmers << '\n';
    }
    printf("Done writing to kmer file.\n");
    solidKmers.close();

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

    bloo1->load_from_kmers(solids_file);
    JChecker* jchecker = new JChecker(j, bloo1);
    KpomerScanner* scanner = new KpomerScanner(solids_file, bloo1, jchecker);
    scanner->scan_all_kpomers(solids_file);
    scanner->printScanSummary();

    if(argc > 6){
        scanner->getJunctionMap()->writeToFile(junctions_filename);
    }
    printf("Program reached end. \n");
    return 0;
}