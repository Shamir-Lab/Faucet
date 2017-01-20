#include "BloomFilter.hpp"
#include <random>
#include <iostream>
#include <vector>
#include <chrono>

int main(int argc, char *argv[])
{	 
/* Initialise. Do this once (not for every
     random number). */
  std::random_device rd;
  std::mt19937_64 gen(rd());

  /* This is where you define the number generator for unsigned long long: */
  std::uniform_int_distribution<uint64_t> dis;
	  
  int size = atoi(argv[1]);

  bloom_parameters param = bloom_parameters(size, .01);
  bloom_filter bloom(param);
  
  int num_kmers = 10000000;
  std::vector<uint64_t> randomKmers;
  randomKmers.reserve(num_kmers);

  /* A few random numbers: */ 
  std::cout << "Generating random kmers.\n"; 
  for (int n=0; n<num_kmers; ++n) 
  	randomKmers[n] = dis(gen);

  std::cout << "Inserting kmers into BF\n";
  auto begin = std::chrono::high_resolution_clock::now();
  for (int n=0; n<num_kmers; ++n) {
	bloom.insert(randomKmers[n]);
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::cout << "10M insertions: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/1000000 << "ms" << std::endl;

  std::cout << "Querying kmers against BF\n";
  begin = std::chrono::high_resolution_clock::now();
  for (int n=0; n<num_kmers; ++n) {
	bloom.contains(randomKmers[n]);
  }
  end = std::chrono::high_resolution_clock::now();
  std::cout << "10M positive queries: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/1000000 << "ms" << std::endl;

  std::cout << "Generating new random kmers.\n"; 
  for (int n=0; n<num_kmers; ++n) 
  	randomKmers[n] = dis(gen);

  std::cout << "Querying kmers against BF\n";
  begin = std::chrono::high_resolution_clock::now();
  for (int n=0; n<num_kmers; ++n) {
	bloom.contains(randomKmers[n]);
  }
  end = std::chrono::high_resolution_clock::now();
  std::cout << "10M negative queries: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/1000000 << "ms" << std::endl;

  std::cout << "Done.";

  return 0;
}