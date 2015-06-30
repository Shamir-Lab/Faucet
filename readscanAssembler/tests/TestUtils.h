#include <string>

#ifndef TEST_UTILS
#define TEST_UTILS

#include "../../utils/Kmer.h"
#include "../../utils/Bloom.h"

using std::string;
extern kmer_type test_kmer;

kmer_type getKmerFromString(std::string kmerString);
bool kmer_matches_readseq(char* read, kmer_type kmer, int i);
bool kmer_matches_kmer(kmer_type kmer1, int i1, kmer_type kmer2, int i2);

Bloom* loadBloom(string list[], int numKmers, int k);

void fail(char* testName, char* errorMessage);
void fail(char* testName);

void succeed(char* testName);

#endif