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

//Returns true if the kmer represented by readKmer is a junction. 
//Uses as a reference "real_ext" since it knows that's one valid path
bool ReadScanner::testForJunction(ReadKmer readKmer){
  kmer_type real_ext = readKmer.getRealExtension();
  
  int branchCount = 1;
  //if thejunction kmer doesn't even j-check backwards, not worth calling it a junction.
  if(readKmer.getMaxGuaranteedJ(!readKmer.direction) < jchecker->j){ //but we should only do the jcheck if the position on the read doesn't already guarantee it'll be fine
    if(!jchecker->jcheck(readKmer.getRevCompKmer())){
      return false;
    }
  }
  //if the real extension doesn't j-check forwards,start count at 0
  if(readKmer.getMaxGuaranteedJ(readKmer.direction) < jchecker->j){ //but we should only do the jcheck if the position on the read doesn't already guarantee it'll be fine
    if(!jchecker->jcheck(real_ext)){
      branchCount = 0;
    }
  }

  //Now given that the real extension is legit, see if there's also an alternate extension
  kmer_type real = readKmer.getKmer();
  for(int nt=0; nt<4; nt++) {//for each extension
    kmer_type test_ext = readKmer.getExtension(nt); //get possible extension
    if(real_ext != test_ext && real_ext != real && test_ext != real){//if the branch has 3 distinct kmers
      if(bloom->oldContains(get_canon(test_ext)))//if the branch checks out initially
      { 
        NbJCheckKmer++;
        if(jchecker->jcheck(test_ext)){//if the branch jchecks
          branchCount++;
          if(branchCount > 1){
            return true;
          }
        }
      }
    }
  }
  return false;
}

//starts at position *pos, kmer *kmer on read, and scans till it either finds an existing junction
//or finds a new one.  Does not modify any junctions- simply updates readKmer to be at a junction or off the end of the read
//True if junction was found
bool ReadScanner::find_next_junction(ReadKmer * readKmer){
  //Iterate forward through the read
  for (; readKmer->getDistToEnd() > 0; readKmer->forward())
  {
      //check for an already found junciton
      if(junctionMap->contains(readKmer->getKmer())){
        return true;
      }
      //check for a new junction
      if(testForJunction(*readKmer)){
        return true;
      }
    if(strcmp(print_kmer(readKmer->getKmer()), "TCACCTTGGCCTCCCAAAGTGCTGGGA")==0){
      printf("Kmer %s is not a junction =(\n", "TCACCTTGGCCTCCCAAAGTGCTGGGA");
    }
      NbProcessed++;
  }
  return false;
}

//should only be called on a read with no real junctions
//Adds a fake junction in the middle and points it to the two ends
void ReadScanner::add_fake_junction(string read){
  ReadKmer* middleKmer = new ReadKmer(&read, read.length()/2- sizeKmer/2, FORWARD);
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
  ReadKmer* readKmer = new ReadKmer(&read);//stores current kmer throughout
  readKmer->forward();//don't scan the first position- points off the read
  ReadKmer* lastJunc;
  
  //handle all junctions - scan forward, find and create junctions, and link adjacent junctions to each other
  while(find_next_junction(readKmer))
  {
    //create a junction at the current spot if none exists
    if(!junctionMap->contains(readKmer)){
      junctionMap->createJunction(readKmer);
    }
    junctionMap->getJunction(readKmer)->addCoverage(readKmer->getRealExtensionNuc());
    
    //if there was a last junction, link the two 
    if(lastJunc){
      junctionMap->directLinkJunctions(lastJunc, readKmer);
    }

    //If this is the first junction, link it to the beginning of the read.
    else{ 
      junctionMap->getJunction(readKmer)
        ->update(readKmer->getExtensionIndex(BACKWARD), readKmer->getTotalPos());
        //2*j because we don't know for sure there isn't a junction 
    }

    //bookkeeping work to advance to the next iteration of the loop
    free(lastJunc);
    lastJunc = new ReadKmer(readKmer);
    int dist = max(1,junctionMap->getSkipDist(readKmer, FORWARD));
    readKmer->advanceDist(dist); //CHANGE BACK TO DIST 

    NbProcessed++,  NbSkipped += dist-1;
  }   
  free(readKmer);

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


