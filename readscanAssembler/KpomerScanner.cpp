#include "KpomerScanner.h"
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <time.h>
using namespace std;

KpomerScanner::KpomerScanner(string inputFile, Bloom* bloo1, JChecker* checker){
  kmers_file = inputFile;
  bloom = bloo1;
  jchecker = checker;
  junctionMap = new JunctionMap();
}


void KpomerScanner::printScanSummary(){
  printf("\n Distinct junctions: %lli \n", (uint64_t)junctionMap->getNumJunctions());
  printf(" Number of kmers that we j-checked: %lli \n", NbJCheckKmer);
  printf (" Number of reads with no junctions: %lli \n",NbNoJuncs);
  printf("Number of spacers: %lli \n", NbSpacers);
  printf("Number of processed kmers: %lli \n", NbProcessed);
}

JunctionMap*  KpomerScanner::getJunctionMap(){
  return junctionMap;
}

void KpomerScanner::scanKpomers(){
  printf("Kpomer scan not implemented yet!");
}