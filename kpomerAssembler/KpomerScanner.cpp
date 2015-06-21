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
  printf("Number of processed kmers: %lli \n", NbProcessed);
}

JunctionMap*  KpomerScanner::getJunctionMap(){
  return junctionMap;
}


void KpomerScanner::scan_kpomer(string kpomer){
  kmer_type thisone, next;
    printf("kpomer: %s\n", &kpomer[0]);
  for(int i = 0; i < 2; i++){
    getFirstKmerFromRead(&thisone,&kpomer[0]);
    int ext = NT2int(kpomer[sizeKmer]);
    //printf("%s \n", print_kmer(left));
    //printf("%s \n", print_kmer(right));
      
    //printf("Strand %d \n", i);
    printf("This: %s \n", print_kmer(thisone));
    if(junctionMap->contains(thisone)){
      junctionMap->getJunction(thisone)->update(ext, 1, 0);
    }
    else{
      for(int nt=0; nt<4; nt++) {
        if(nt != ext){
          next = next_kmer(thisone,nt, 0);
          printf("Next: %s \n", print_kmer(next));
          if(bloom->oldContains(get_canon(next))){ 
              NbCandKmer++;
              if(jchecker->jcheck(next)){
                junctionMap->createJunction(thisone, nt); 
              }                     
          }
          //printf("\n");
        }
      }
    } 
    revcomp_sequence(&kpomer[0], sizeKmer+1);
  }
}


int KpomerScanner::scan_all_kpomers(string kmers_file)
{

  ifstream solidKmers;
  solidKmers.open(kmers_file);

  string kpomer;
 
  // write all positive extensions in disk file
  while (getline(solidKmers, kpomer))
  {
    scan_kpomer(kpomer);
    NbProcessed++;
  }

  solidKmers.close();
}