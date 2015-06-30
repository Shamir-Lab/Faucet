#include "ReadScanner.h"
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <time.h>
using namespace std;

ReadScanner::ReadScanner(string readFile, Bloom* bloo1, JChecker* checker){
  reads_file = readFile;
  bloom = bloo1;
  junctionMap = new JunctionMap();
  jchecker = checker;
}

void ReadScanner::resetHashes(kmer_type kmer){
  hash0 = bloom->get_rolling_hash(kmer,0);
  hash1 = bloom->get_rolling_hash(kmer,1);
}


void ReadScanner::printScanSummary(){
  printf("\n Distinct junctions: %lli \n", (uint64_t)junctionMap->getNumJunctions());
  printf(" Number of kmers that we j-checked: %lli \n", NbJCheckKmer);
  printf (" Number of reads with no junctions: %lli \n",NbNoJuncs);
  printf("Number of spacers: %lli \n", NbSpacers);
  printf("Total number of kmers: %lli \n", (uint64_t)(readLength-sizeKmer)*(uint64_t)readsProcessed*(uint64_t)2);
  printf("Number of processed kmers: %lli \n", NbProcessed);
}

JunctionMap* ReadScanner::getJunctionMap(){
  return junctionMap;
}

//Returns true if the kmer "real" is a junction. Uses as a reference "real_ext" since it knows that's one valid path
bool ReadScanner::testForJunction(DoubleKmer doubleKmer){
  kmer_type real_ext = doubleKmer.getRealExtension();
  kmer_type real = doubleKmer.getKmer();
  for(int nt=0; nt<4; nt++) {//for each extension
    kmer_type test_ext = doubleKmer.getExtension(nt); //get possible extension
    if(real_ext != test_ext && real_ext != real && test_ext != real){//if the branch has 3 distinct kmers
      if(bloom->oldContains(get_canon(test_ext)))//if the branch checks out initially
      { 
        NbJCheckKmer++;
        if(jchecker->jcheck(test_ext)){//if it j-checks- new junction! //To make rolling, fix this to use proper jcheck!
            return true;
        }
      }
    }
  }
}

//starts at position *pos, kmer *kmer on read, and scans till it either finds an existing junction
//or finds a new one.  Does not modify any junctions- simply updates doubleKmer to be at a junction or off the end of the read
//True if junction was found
bool ReadScanner::find_next_junction(DoubleKmer * doubleKmer){
  //Iterate forward through the read
  for (; doubleKmer->onRead(); doubleKmer->forward())
  {
      //check for an oldjunciton
      if(junctionMap->contains(doubleKmer->getKmer())){
        return true;
      }
      //check for a new junction
      if(testForJunction(*doubleKmer)){
        return true;
      }
  }
  //Didn't find a junction. Return null
  return false;
}

void ReadScanner::scan_forward(string read){
  int pos = 0;
  kmer_type kmer;
  DoubleKmer* doubleKmer = new DoubleKmer(&read);//stores current kmer throughout
  DoubleKmer* lastJunc;
  //printf("Scanning read: %s\n", &read[0]);
  while(find_next_junction(doubleKmer))
  {
    //printf("Pos: %d, Dir: %s, %s \n", doubleKmer->pos,doubleKmer->directionAsString(), print_kmer(doubleKmer->getKmer()));
    if(!junctionMap->contains(doubleKmer->getKmer())){//need to create a junction
      junctionMap->createJunction(doubleKmer->getKmer());
    }
    junctionMap->getJunction(doubleKmer->getKmer())->addCoverage(doubleKmer->getRealExtensionNuc());
    if(lastJunc){

      int ext1 = lastJunc->getRealExtensionNuc(), ext2 = doubleKmer->getRealExtensionNuc();
      if(lastJunc->direction == BACKWARD){
        ext1 = 4;
      }
      if(doubleKmer->direction == FORWARD){
        ext2 = 4;
      }
      int dist = 2*(doubleKmer->pos - lastJunc->pos) + doubleKmer->offset() - lastJunc->offset();
      junctionMap->linkJunctions(lastJunc->getKmer(),ext1, doubleKmer->getKmer(), ext2, dist);

    }
    else{ 
      //LINK JUNCTION TO SPACER IF NECESSARY
    }
    lastJunc = new DoubleKmer(doubleKmer);
    NbProcessed++;

    int skipDist;
    if(doubleKmer->direction == FORWARD){
      skipDist = max(1,junctionMap->getJunction(doubleKmer->getKmer())->dist[doubleKmer->getRealExtensionNuc()]);
    }
    else{
     skipDist = max(1,junctionMap->getJunction(doubleKmer->getKmer())->dist[4]);
    } 
    doubleKmer->advanceDist(skipDist);
  }   
  //LINK JUNCTION TO SPACER IF NECESSARY
  free(doubleKmer);
  free(lastJunc);
}


void ReadScanner::scanReads()
{
  NbCandKmer=0, NbRawCandKmer = 0, NbJCheckKmer = 0, NbNoJuncs = 0, 
  NbSkipped = 0, NbProcessed = 0, readsProcessed = 0, NbSolidKmer =0;

  time_t start;
  time_t stop;
  time(&start);

  ifstream solidReads;
  solidReads.open(reads_file);

  string read;

  // write all positive extensions in disk file
  printf("Weight before read scan: %f \n", bloom->weight());
  int lastSum = 0, thisSum = 0;
  while (getline(solidReads, read))
  {
    //lastSum = thisSum;
    readLength = read.length();
    scan_forward(read);
    if ((readsProcessed%10000)==0) fprintf (stderr,"Reads processed: %c %lld",13,(long long)readsProcessed);
    readsProcessed++;
  }

  solidReads.close();
  time(&stop);
  printf("Time in seconds for read scan: %f \n", difftime(stop,start));
}


