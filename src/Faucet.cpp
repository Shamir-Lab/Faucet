    
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
#include "Faucet.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

#define NNKS 4 // default minimal abundance for solidity
#define MIN_CONTIG_SIZE (2*sizeKmer+1)

/*
To run Faucet, first type make in the directory to compile.

Then, type ./faucet, followed by the following arguments:
-read_load_file <>, name of reads file  to load bloom filter from (current format is each line has a string of characters representing the read)
-read_scan_file <>, name of reads file  for read scan (current format is each line has a string of characters representing the read)
-size_kmer k
-max_read_length <>, upper bound on the size of a read
-estimated_kmers <>, number of number of distinct kmers.  This will be directly used to size the bloom filter so try to have a good estimate.
-fp <>, false positive rate, default .01
-file_prefix <>, used for junctions file and contigs file and graph file
--two_hash, if this option is selected a bigger bloom filter with only two hash functions is used
-bloom_file <>, used to shortcut loading the filter if you already have it on file
-junctions_file <>, used to shortcut the readscan if you have access to a junctions file
--just_load_bloom, if this option is selected the bloom will be loaded and dumped, then the program will terminate
--fastq, use fastq files
--paired_ends, file is given as interleaved paired end data.  Beginning of each read corresponds to end of overall fragment.

Note: cannot use junctions_file option without also using bloom_file option

This will load a bloom filter with all the kmers from the reads, then scan through them to find all of the junctions.
It will print all of the junctions to a file named file_prefix.junctions, 
and print some summary info about the run to stdout.
It will then build a graph from the junction map, clean tips from the graph, and print the contigs and the graph to files
file_prefix.contigs, file_prefix.graph.
*/

void argumentError(){
    fprintf (stderr,"Usage:\n");
    fprintf (stderr,"./faucet -read_load_file <filename> -read_scan_file <filename> -size_kmer <k> -max_read_length <length> -estimated_kmers <num_kmers> -singletons <num_kmers> -file_prefix <prefix>");
    fprintf(stderr, "\nOptional arguments: --fastq --mercy --high_cov -max_spacer_dist <dist> -fp rate <rate> -j <int> --two_hash -bloom_file <filename> -junctions_file <filename> --paired_ends --no_cleaning\n");
}


int handle_arguments(int argc, char *argv[]){
    if (argc==1){
        argumentError();
        return 1;
    }
    for(int i = 1; i < argc; i++){
        if(0 == strcmp(argv[i] , "-read_load_file")) //1st arg: read file name
                read_load_file = string(argv[i+1]), i++, load_file_flag=true;
        else if(0 == strcmp(argv[i] , "-read_scan_file")) //1st arg: read file name
                read_scan_file = string(argv[i+1]), i++, scan_file_flag=true;
        else if(0 == strcmp(argv[i] , "-size_kmer")) // 2rd arg: kmer size.
                setSizeKmer(atoi(argv[i+1])), i++, k_val_flag=true;
        else if(0 == strcmp(argv[i], "-max_read_length")) //3rd arg: read length
                read_length = atoi(argv[i+1]), i++, max_len_flag=true;
        else if(0 == strcmp(argv[i] , "-estimated_kmers")) //estimated number of distinct kmers
                estimated_kmers = atoll(argv[i+1]), i++, est_kmers_flag=true;
        else if(0 == strcmp(argv[i] , "-singletons")) //estimated number of distinct kmers
                singletons = atoll(argv[i+1]), i++, est_sing_flag=true;    
        else if(0 == strcmp(argv[i] , "-fp")) //false posiive rate
                fpRate = atof(argv[i+1]), i++;
        else if(0 == strcmp(argv[i] , "-j")) //value of j for jchecking
                j = atoi(argv[i+1]), i++;
        else if(0 == strcmp(argv[i] , "-file_prefix")) //file prefix used for output files
                file_prefix = string(argv[i+1]), i++, pref_flag=true;
        else if(0 == strcmp(argv[i] , "--two_hash")) //use two hash function BF instead of optimal size BF
                two_hash = true; 
        else if(0 == strcmp(argv[i] , "--just_load_bloom")) //stop after loading bloom
                just_load = true;
        else if(0 == strcmp(argv[i] , "--no_cleaning")) //stop after building contigmap
                no_cleaning = true;            
        else if(0 == strcmp(argv[i] , "--fastq")) //input is fastq file
                fastq = true;
        else if(0 == strcmp(argv[i] , "--mercy")) //input is fastq file
                mercy = true;            
        else if(0 == strcmp(argv[i] , "--high_cov")) //input is likely high coverage
                high_cov = true;            
        else if(0 == strcmp(argv[i], "--node_graph")) //use node graph rather than contig graph
                node_graph = true; 
        else if(0 == strcmp(argv[i], "--paired_ends")) //assume interleaved paired end reads
                paired_ends = true;
        else if(0 == strcmp(argv[i] , "-bloom_file")){
                bloom_input_file = string(argv[i+1]);
                from_bloom = true, i++;
        }  
        else if(0 == strcmp(argv[i] , "-max_spacer_dist")){
                maxSpacerDist = atoi(argv[i+1]), i++;
        }  
        else if(0 == strcmp(argv[i] , "-junctions_file")){
                junctions_input_file = string(argv[i+1]) + ".junctions";
                short_pair_filter_file = string(argv[i+1]) + ".short_pair_filter";
                long_pair_filter_file = string(argv[i+1]) + ".long_pair_filter";
                from_junctions = true, i++;
        }
        else if (0 == strcmp(argv[i] , "--help") || 0 == strcmp(argv[i] , "-h")){
            argumentError();
            return 1;
        }
        else if (0 == strcmp(argv[i] , "--help") || 0 == strcmp(argv[i] , "-h")){
            argumentError();
            return 1;
        }
        else {
            fprintf (stderr, "Cannot parse tag %s\n", argv[i]);
            argumentError();
            return 1; 
        }
        
    }
    if (! (load_file_flag && scan_file_flag && k_val_flag && max_len_flag && est_kmers_flag && est_sing_flag && pref_flag)){
        fprintf (stderr, "Some required argument is missing.\n");
        argumentError();
        return 1; 
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
        printf("Starting at the beginning: will load bloom and find junctions from the read set.\n");

    if(just_load)
        printf("Only loading bloom, dumping and termination.\n");

    if(node_graph)
        printf("Using node graph.\n");
    else 
        printf("Using contig graph.\n");

    std::cout << "Read load file name: " << read_load_file << "\n";

    std::cout << "Read scan file name: " << read_scan_file << "\n";

    printf("k: %d: \n", sizeKmer);

    printf("Maximal read length: %d\n", read_length);

    printf("Estimated number of distinct kmers, for sizing bloom filter: %lli.\n", estimated_kmers);

    printf("False positive rate: %f\n", fpRate);

    // printf("j: %d \n", j);
    
    printf("File prefix: %s\n", &file_prefix[0]);

    printf("Max spacer dist: %d\n", maxSpacerDist);
    
    if(two_hash){
        printf("Using 2 hash functions.\n");
    }
    else{
        printf("Using space-optimal hash settings.\n");
    }

    std::cout <<  "Paired ends: " << paired_ends << "\n";
    printf("Size of junction: %d\n", sizeof(Junction));
    printf("Size of contigNode: %d\n", sizeof(ContigNode));
    printf("Size of contig: %d\n", sizeof(Contig));
    printf("Size of int: %d\n", sizeof(int));
    printf("Size of long: %d\n", sizeof(long));
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

double my_func(double p1) { 
    double c = (estimated_kmers - (1-p1)*singletons) /estimated_kmers;
    return log(2)*log(fpRate) + log(p1)*log(1 - pow(2, -c));
    // return (log(fpRate)* (estimated_kmers - (1-p1)*singletons) /estimated_kmers) - log(p1);
}
     

Bloom* getBloomFilterFromReads(){ //handles loading from reads
    Bloom* bloo1;
    Bloom* bloo2;

    std::function<double (double)> f = my_func;
    double p1 = brents_fun(f, fpRate, 0.50, 0.0001, 1000);
    cout << "p2 is " << fpRate << " p1 estimated as " << p1 << endl;
    
    // if(two_hash){
    //     bloo1 = bloo1->create_bloom_filter_2_hash(estimated_kmers, fpRate);
    //     bloo2 = bloo2->create_bloom_filter_2_hash(estimated_kmers, fpRate);
    // }
    // else{
    bloo1 = bloo1->create_bloom_filter_optimal(estimated_kmers, p1);
    bloo2 = bloo2->create_bloom_filter_optimal(estimated_kmers, p1);
    // }
    load_two_filters(bloo1, bloo2, read_load_file, fastq, mercy);
    delete(bloo1);
    return bloo2;
}


Bloom* getBloomFilterFromReadsSingle(){ //handles loading from reads
    Bloom* bloo1;

    if(two_hash){
        bloo1 = bloo1->create_bloom_filter_2_hash(estimated_kmers, fpRate);
    }
    else{
        bloo1 = bloo1->create_bloom_filter_optimal(estimated_kmers, fpRate);
    }
    load_single_filter(bloo1, read_load_file, fastq);
    return bloo1;
}

//Builds the junction map from either a file or the readscan
void buildJunctionMapFromReads(JunctionMap* junctionMap, Bloom* bloom, Bloom* short_pair_filter, Bloom* long_pair_filter, JChecker* jchecker, bool no_cleaning){
    ReadScanner* scanner = new ReadScanner(junctionMap, read_scan_file, bloom, short_pair_filter, long_pair_filter, jchecker, maxSpacerDist);
     
    //scan reads, print summary
    scanner->scanReads(fastq, paired_ends, no_cleaning);
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
        bloom->dump(&(file_prefix + ".bloom")[0]);
    }
    
    Bloom* short_pair_filter;
    Bloom* long_pair_filter;
    if (!high_cov){
        short_pair_filter = short_pair_filter->create_bloom_filter_optimal(estimated_kmers/20, 0.01);
        if (paired_ends) long_pair_filter = long_pair_filter->create_bloom_filter_optimal(estimated_kmers/10, 0.01);
        else long_pair_filter = nullptr;
    }else {
        short_pair_filter = short_pair_filter->create_bloom_filter_optimal(estimated_kmers/2, 0.01);
        if (paired_ends) long_pair_filter = long_pair_filter->create_bloom_filter_optimal(estimated_kmers/2, 0.01);
        else long_pair_filter = nullptr;
    }
    if(just_load) return 0;
    
    //create JChecker
    JChecker* jchecker = new JChecker(j, bloom);

    //Build junction map from either reads or file
    JunctionMap* junctionMap = new JunctionMap(bloom, jchecker, read_length);
    if(from_junctions){
        junctionMap->buildFromFile(junctions_input_file);
        short_pair_filter->load(&short_pair_filter_file[0]);
        long_pair_filter->load(&long_pair_filter_file[0]);
    }
    else{
        buildJunctionMapFromReads(junctionMap, bloom, short_pair_filter, long_pair_filter, jchecker, no_cleaning);
        junctionMap->writeToFile(file_prefix + ".junctions");
        if(!no_cleaning){
            short_pair_filter->dump(&(file_prefix + ".short_pair_filter")[0]);
            if (paired_ends) long_pair_filter->dump(&(file_prefix + ".long_pair_filter")[0]);
        }
    }
    //dump junctions to file
    
    printf("Weight of short pair filter: %f\n", short_pair_filter->weight());
    if (paired_ends) printf("Weight of long pair filter: %f\n", long_pair_filter->weight());
    printf("Number of junctions: %d\n", junctionMap->junctionMap.size());
    ContigGraph* contigGraph = junctionMap->buildContigGraph();
    contigGraph->setReadLength(read_length);
    delete(bloom);

    // contigGraph->checkGraph();
    // contigGraph->printGraph(file_prefix + ".raw_graph.fastg");
    // contigGraph->printContigs(file_prefix + ".raw_contigs.fasta");

    if (no_cleaning){ return 0; }
    while(contigGraph->cleanGraph(short_pair_filter, long_pair_filter));

    contigGraph->checkGraph();
    contigGraph->printContigs(file_prefix + ".cleaned_contigs.fasta");
    contigGraph->printGraph(file_prefix + ".cleaned_graph.fastg");

    printf("Program reached end. \n");
    return 0;
}