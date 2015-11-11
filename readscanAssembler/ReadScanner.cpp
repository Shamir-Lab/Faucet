#include "ReadScanner.h"
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <time.h>
using namespace std;

ReadScanner::ReadScanner(JunctionMap* juncMap, string readFile, Bloom* bloo1, Bloom* junc_bloom, JChecker* checker){
  reads_file = readFile;
  bloom = bloo1;
  pair_filter = junc_bloom;
  junctionMap = juncMap;
  jchecker = checker;
}

void ReadScanner::printScanSummary(){
  printf("\nDistinct junctions: %lli \n", (uint64_t)junctionMap->getNumJunctions());
  //printf("Number of junction pairs that exist on reads: %d\n", juncPairSet.size());
  printf("Number of kmers that we j-checked: %lli \n", NbJCheckKmer);
  printf ("Number of reads with no junctions: %lli \n",NbNoJuncs);
  printf("Number of processed kmers: %lli \n", NbProcessed);
  printf("Number of skipped kmers: %lli \n", NbSkipped);
  printf("Reads without errors: %lli\n", readsNoErrors);
}

JunctionMap* ReadScanner::getJunctionMap(){
  return junctionMap;
}

//Returns true if the kmer represented by readKmer is a junction. 
//Uses as a reference "real_ext" since it knows that's one valid path
//Should only be called on a kmer at least j away from the end of the read
bool ReadScanner::testForJunction(ReadKmer readKmer){
  kmer_type real_ext = readKmer.getRealExtension();
  
  //Check alternate extensions, and if the total valid extension count is greater than 1, return true. 
  kmer_type real = readKmer.getKmer();
  for(int nt=0; nt<4; nt++) {//for each extension
    kmer_type test_ext = readKmer.getExtension(nt); //get possible extension
    if(real_ext != test_ext){//if the alternate and real extensions are different- note that I took out the other two distinction checks
      if(bloom->oldContains(get_canon(test_ext)))//if the branch checks out initially
      { 
        NbJCheckKmer++;
        if(jchecker->jcheck(test_ext)){//if the branch jchecks
            return true;
        }
      }
    }
  }
  return false;  
}

//starts at the given ReadKmer, and scans till it either finds an existing junction, finds a new junction, or hits the end of the read.
//Does not modify any junctions- simply updates ReadKmer to be at a junction or off the end of the read
//Returns true if junction was found
bool ReadScanner::find_next_junction(ReadKmer * readKmer){
  //Iterate forward to the end of the read
  for (; readKmer->getDistToEnd() > 2*jchecker->j; readKmer->forward()) //CHANGED TO 2*j from 0
  {
      //check for an already found junciton
      if(junctionMap->isJunction(readKmer->getKmer())){
        return true;
      }
      //check for a new junction
      if(testForJunction(*readKmer)){
        return true;
      }
      NbProcessed++;
  }
  return false;
}

//Should only be called on a read with no real junctions
//Adds a fake junction in the middle and points it to the two ends.  This ensures we have coverage of long linear regions, and that we capture
//sinks at the end of such regions.
void ReadScanner::add_fake_junction(string read){
  ReadKmer* middleKmer = new ReadKmer(&read, read.length()/2- sizeKmer/2, FORWARD);
  junctionMap->createJunction(middleKmer);
  Junction* junc = junctionMap->getJunction(middleKmer);
  junc->addCoverage(middleKmer->getRealExtensionNuc());
  junc->update(middleKmer->getExtensionIndex(BACKWARD), middleKmer->getTotalPos()-2*jchecker->j);
  junc->update(middleKmer->getExtensionIndex(FORWARD), middleKmer->getDistToEnd()-2*jchecker->j);
  delete(middleKmer);
}


//Scans a read. 
//Identifies all junctions on the read, and links adjacent junctions to each other.
//Also updates the relevant distance field on the first junction to point to the start of the read, and on the last
//Junction to point to the end of the read.
//If there are no junctions, add_fake_junction is called
std::list<kmer_type> ReadScanner::scan_forward(string read){
  std::list<kmer_type> result = {};
  ReadKmer* readKmer = new ReadKmer(&read);//stores current kmer throughout
  backwardSet = {};

  for(int i = 0; i < 2*jchecker->j + 1; i++){
    readKmer->forward();  
  }

  ReadKmer* backJunc = nullptr;
  ReadKmer* firstBackJunc = nullptr;
  ReadKmer* lastForwardJunc = nullptr;
  ReadKmer* lastKmer = nullptr;
  Junction* lastJunc = nullptr;
  Junction* junc = nullptr;
  int numBackward = 0;

  //handle all junctions - scan forward, find and create junctions, and link adjacent junctions to each other
  while(find_next_junction(readKmer))
  { 
    // mark first R junction and last F junction
    if(readKmer->direction == BACKWARD){
      if(!backJunc){
        backJunc = new ReadKmer(readKmer);
        firstBackJunc = new ReadKmer(readKmer);
      }
      else{
        *backJunc = readKmer;
      }
      result.push_back(backJunc->getRealExtension());
    }
    else{
      if(backJunc){
         pair_filter->addPair(backJunc->getRealExtension(), readKmer->getRealExtension());
      }
      if(!lastForwardJunc){
        lastForwardJunc = new ReadKmer(readKmer);
      }
      else{
        *lastForwardJunc = readKmer;
      }
    }
   
    junc = junctionMap->getJunction(readKmer);
    //create a junction at the current spot if none exists
    if(!junc){
      junctionMap->createJunction(readKmer);
      junc = junctionMap->getJunction(readKmer);
    }

    junc->addCoverage(readKmer->getRealExtensionNuc()); //add coverage of the junction
    
    //if there was a last junction, link the two 
    if(lastKmer){
      junctionMap->directLinkJunctions(lastKmer, readKmer, lastJunc, junc);
    }
    //If this is the first junction, link it to the beginning of the read.
    else{ 
      lastKmer = new ReadKmer(readKmer);
      junc
        ->update(readKmer->getExtensionIndex(BACKWARD), readKmer->getTotalPos()-2*jchecker->j);//-2*j ADDED
    }

    *lastKmer = *readKmer;
    lastJunc = junc;

    int index = readKmer->getExtensionIndex(FORWARD);
    int dist = max(1, (int)(junc->dist[index]));

    readKmer->advanceDist(dist); 

    NbProcessed++,  NbSkipped += dist-1;
  }   
  
  //If there were no junctions on the read, add a fake junction
  if(!lastKmer){
      NbNoJuncs++;
      add_fake_junction(read);
  }
  //If there was at least one junction, point the last junction found to the end of the read
  else {
    lastJunc->update(lastKmer->getExtensionIndex(FORWARD), lastKmer->getDistToEnd()-2*jchecker->j); //2*j ADDED
  }

  if(firstBackJunc && lastForwardJunc){
    //pair_filter->addPair(firstBackJunc->getRealExtension(), lastForwardJunc->getRealExtension());
  }

  delete(lastKmer);
  delete(readKmer);
  delete(firstBackJunc);
  delete(lastForwardJunc);

  return result;
}

bool ReadScanner::isValidRead(string read){
  //Move to the first valid kmer
  for(ReadKmer kmer = ReadKmer(&read); kmer.getDistToEnd()>= 0; kmer.forward(), kmer.forward()){
    if(!bloom->oldContains(kmer.getCanon())){
      return false;
    } 
  }

  return true;
}

//returns back junctions for use in paired end info
std::list<kmer_type> ReadScanner::scanInputRead(string read){
  std::list<kmer_type> result = {};
  std::list<string> readList = getUnambiguousReads(read);
    while(!readList.empty()){
        read = readList.front();
        readList.pop_front();
        if(read.length() >= sizeKmer + 2*jchecker->j + 1){
          unambiguousReads++;
            //printf("Checking for errors.\n");
          if(isValidRead(read)){
            //printf("None! Scanning\n");
            result.splice(result.end(),scan_forward(read));
            readsNoErrors++;
          }
        }
    }
    return result;
}

void ReadScanner::scanReads(bool fastq, bool paired_ends)
{
  NbCandKmer=0, NbRawCandKmer = 0, NbJCheckKmer = 0, NbNoJuncs = 0, 
  NbSkipped = 0, NbProcessed = 0, readsProcessed = 0, NbSolidKmer =0, 
  readsNoErrors = 0,  NbJuncPairs = 0, unambiguousReads = 0;

  time_t start;
  time_t stop;
  time(&start);

  ifstream solidReads;
  solidReads.open(reads_file);

  string read;

  // write all positive extensions in disk file
  printf("Weight before read scan: %f \n", bloom->weight());
  int lastSum = 0, thisSum = 0;
  bool firstEnd = true;
  std::list<kmer_type> backJuncs1;
  std::list<kmer_type> backJuncs2;
  while (getline(solidReads, read))
  {
    getline(solidReads, read);//since it's a fasta we skip the first of every pair of lines
    
    if(firstEnd){
      backJuncs1 = scanInputRead(read);
    }
    else{
      backJuncs2 = scanInputRead(read);
    }
    
    if(paired_ends && !firstEnd){
      for(auto it = backJuncs1.begin(); it != backJuncs1.end(); it++){
          kmer_type pair1 = *it;
          for(auto it2 = backJuncs2.begin(); it2 != backJuncs2.end(); it2++){
            kmer_type pair2 = *it2;
            pair_filter->addPair(pair1,pair2);
          }
      }
    }

    if ((readsProcessed%10000)==0) fprintf (stderr,"Reads processed: %c %lld",13,(long long)readsProcessed);
    readsProcessed++;
    if(fastq) getline(solidReads, read), getline(solidReads, read); //if fastq skip two more lines
    firstEnd = !firstEnd;
  }

  solidReads.close();
  time(&stop);
  printf("Reads processed: %lli\n", readsProcessed);
  printf("Unambiguous reads: %lli\n", unambiguousReads);
  printf("Time in seconds for read scan: %f \n", difftime(stop,start));
}


