#include "ReadScanner.h"
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <time.h>
using namespace std;

ReadScanner::ReadScanner(JunctionMap* juncMap, string readFile, Bloom* bloo1, JChecker* checker){
  reads_file = readFile;
  bloom = bloo1;
  junctionMap = juncMap;
  jchecker = checker;
}

void ReadScanner::resetHashes(kmer_type kmer){
  hash0 = bloom->get_rolling_hash(kmer,0);
  hash1 = bloom->get_rolling_hash(kmer,1);
}

void ReadScanner::printScanSummary(){
  printf("\n Distinct junctions: %lli \n", (uint64_t)junctionMap->getNumJunctions());
  printf("Number of kmers that we j-checked: %lli \n", NbJCheckKmer);
  printf ("Number of reads with no junctions: %lli \n",NbNoJuncs);
  printf("Number of spacers: %lli \n", NbSpacers);
  printf("Number of processed kmers: %lli \n", NbProcessed);
  printf("Number of skipped kmers: %lli \n", NbSkipped);
}

JunctionMap* ReadScanner::getJunctionMap(){
  return junctionMap;
}

//Returns true if the kmer represented by doubleKmer is a junction. 
//Uses as a reference "real_ext" since it knows that's one valid path
bool ReadScanner::testForJunction(DoubleKmer doubleKmer){
  kmer_type real_ext = doubleKmer.getRealExtension();
  kmer_type real = doubleKmer.getKmer();
  for(int nt=0; nt<4; nt++) {//for each extension
    kmer_type test_ext = doubleKmer.getExtension(nt); //get possible extension
    if(real_ext != test_ext && real_ext != real && test_ext != real){//if the branch has 3 distinct kmers
      if(bloom->oldContains(get_canon(test_ext)))//if the branch checks out initially
      { 
        NbJCheckKmer++;
        if(jchecker->jcheck(test_ext)){
          return true;
        }
      }
    }
  }
  return false;
}

//starts at position *pos, kmer *kmer on read, and scans till it either finds an existing junction
//or finds a new one.  Does not modify any junctions- simply updates doubleKmer to be at a junction or off the end of the read
//True if junction was found
bool ReadScanner::find_next_junction(DoubleKmer * doubleKmer){
  //Iterate forward through the read
  for (; doubleKmer->getDistToEnd() > 0; doubleKmer->forward())
  {
      //check for an already found junciton
      if(junctionMap->contains(doubleKmer->getKmer())){
        return true;
      }
      //check for a new junction
      if(testForJunction(*doubleKmer)){
        return true;
      }
      NbProcessed++;
  }
  return false;
}

//should only be called on a read with no real junctions
//Adds a fake junction in the middle and points it to the two ends
void ReadScanner::add_fake_junction(string read){
  DoubleKmer* middleKmer = new DoubleKmer(&read, read.length()/2- sizeKmer/2, FORWARD);
  junctionMap->createJunction(middleKmer);
  Junction* junc = junctionMap->getJunction(middleKmer);
  junc->addCoverage(middleKmer->getExtensionIndex(FORWARD));
  junc->update(middleKmer->getExtensionIndex(BACKWARD), middleKmer->getTotalPos());
  junc->update(middleKmer->getExtensionIndex(FORWARD), middleKmer->getDistToEnd());
  free(middleKmer);
}

void ReadScanner::scan_forward(string read){
  int pos = 0;
  kmer_type kmer;
  DoubleKmer* doubleKmer = new DoubleKmer(&read);//stores current kmer throughout
  doubleKmer->forward();//don't scan the first position- points off the read
  DoubleKmer* lastJunc;
  
  //handle all junctions - scan forward, find and create junctions, and link adjacent junctions to each other
  while(find_next_junction(doubleKmer))
  {
    //create a junction at the current spot if none exists
    if(!junctionMap->contains(doubleKmer)){
      junctionMap->createJunction(doubleKmer);
    }
    junctionMap->getJunction(doubleKmer)->addCoverage(doubleKmer->getRealExtensionNuc());
    
    //if there was a last junction, link the two 
    if(lastJunc){
      junctionMap->directLinkJunctions(lastJunc, doubleKmer);
    }

    //If this is the first junction, link it to the beginning of the read.
    else{ 
      junctionMap->getJunction(doubleKmer)
        ->update(doubleKmer->getExtensionIndex(BACKWARD), doubleKmer->getTotalPos());
    }

    //bookkeeping work to advance to the next iteration of the loop
    free(lastJunc);
    lastJunc = new DoubleKmer(doubleKmer);
    int dist = max(1,junctionMap->getSkipDist(doubleKmer, FORWARD));
    doubleKmer->advanceDist(dist);

    NbProcessed++,  NbSkipped += dist-1;
  }   
  free(doubleKmer);

  //If there were no juncs- put a fake junction in the middle, point it to the ends.
  if(!lastJunc){
      NbNoJuncs++;
      add_fake_junction(read);
  }

  //If there was at least one junction, point the last junction found to the end of the read
  else {
    Junction* junc = junctionMap->getJunction(lastJunc->getKmer());
    junc->update(lastJunc->getExtensionIndex(FORWARD), lastJunc->getDistToEnd());
  }

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


