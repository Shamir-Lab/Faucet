#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <inttypes.h>
#include <cmath> // for log2f
#include <algorithm> // for max
#include <unistd.h> // for truncate
#include <string>
#include <unordered_map>
using std::string;
using std::unordered_map;

#ifndef DEBLOOM_H
#define DEBLOOM_H

#include "Bloom.h"
#include "Kmer.h"

#define DEBUGE(a)  //printf a

extern unordered_map<kmer_type, unsigned char *> junctionMap;
extern uint64_t  b1_size ;
extern uint64_t nbkmers_solid ;

void allocateJunctionMap(uint64_t size);
void junctionMapToFile(string filename);
bool jcheck(kmer_type kmer, int j, int strand, Bloom* bloo1);
bool find_next_junction(int* pos, kmer_type * kmer, string read, int j, Bloom* bloo1);
void smart_traverse_read(string read, Bloom* bloo1, int j);
int debloom_kpomerscan(char* solids_file, Bloom* bloo1, int j);
int debloom_readscan(char* solids_file, Bloom* bloo1, int j, int genome_size);
void end_debloom_partition(bool last_partition);

#endif
