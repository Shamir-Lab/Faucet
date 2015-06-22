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
    //printf("kpomer: %s\n", &kpomer[0]);
  for(int i = 0; i < 2; i++){
    getFirstKmerFromRead(&thisone,&kpomer[0]);
    int exte = NT2int(kpomer[sizeKmer]);
    //printf("%s \n", print_kmer(left));
    //printf("%s \n", print_kmer(right));
      
    //printf("Strand %d \n", i);
    //printf("This: %s \n", print_kmer(thisone));
    if(junctionMap->contains(thisone)){
      Junction* j = junctionMap->getJunction(thisone);
      // printf("Old junction: %s \n", print_kmer(thisone));
      // printf("Existing extensions: %d %d %d %d \n", j->ext[0], j->ext[1], j->ext[2], j->ext[3]);
      // printf("Added extension %d. \n", exte);
      // if(j->ext[exte] == 0){
      //   printf("%d made it complex\n", exte);
      // }
      j->update(exte, 1, 0);
      //printf("After extensions: %d %d %d %d \n", j->ext[0], j->ext[1], j->ext[2], j->ext[3]);
 
    }
    else{
      for(int nt=0; nt<4; nt++) {
        if(nt != exte){
          next = next_kmer(thisone,nt, 0);
         // printf("Next: %s \n", print_kmer(next));
          if(bloom->oldContains(get_canon(next))){ 
              NbCandKmer++;
              if(jchecker->jcheck(next)){
                //printf("New junction: %s \n", print_kmer(thisone));
                //printf("Real extension: %d \n", exte);
                //printf("Alternate extension: %d \n", nt);
                junctionMap->createJunction(thisone, exte); 
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