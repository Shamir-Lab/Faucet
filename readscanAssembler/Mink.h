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

float fpRate = .01;
int j = 0;

string read_load_file;
string read_scan_file;
string bloom_input_file;
string junctions_input_file;
string short_pair_filter_file;
string long_pair_filter_file;

int read_length;
uint64_t estimated_kmers;
bool two_hash = false;
bool from_bloom = false;
bool from_junctions = false;
bool just_load = false;
bool fastq = false;
bool node_graph = false;
bool paired_ends = false;
bool no_cleaning = false;
int maxSpacerDist = 100; //max is 128, smaller --> more frequent spacers, bigger --> less frequent.  Measured in base pairs
int64_t nb_reads;

set<kmer_type> all_kmers;
string file_prefix;

#endif