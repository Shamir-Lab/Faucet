#ifndef MINK_MAIN
#define MINK_MAIN

#include "../utils/Bloom.h"
#include "../utils/Kmer.h"
#include "../utils/Junction.h"
#include "../utils/JChecker.h"
#include "ReadScanner.h"
#include "ContigNode.h"
#include "Contig.h"
#include "ReadScanner.h"
#include "ContigGraph.h"

float fpRate = .04;
int j = 1;

string read_load_file;
string read_scan_file;
string bloom_input_file;
string junctions_input_file;
string short_pair_filter_file;
string long_pair_filter_file;

int read_length;
uint64_t estimated_kmers;
uint64_t singletons;


// requred arguments:
bool load_file_flag = false;
bool scan_file_flag = false;
bool k_val_flag = false;
bool max_len_flag = false;
bool est_kmers_flag = false;
bool est_sing_flag = false;
bool pref_flag = false;

// optional arguments:
bool two_hash = false;
bool from_bloom = false;
bool from_junctions = false;
bool just_load = false;
bool fastq = false;
bool mercy = false;
bool node_graph = false;
bool paired_ends = false;
bool no_cleaning = false;
int maxSpacerDist = 100; //max is 128, smaller --> more frequent spacers, bigger --> less frequent.  Measured in base pairs
int64_t nb_reads;
bool high_cov = false;

set<kmer_type> all_kmers;
string file_prefix;

#endif