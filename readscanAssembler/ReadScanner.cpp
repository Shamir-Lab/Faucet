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
            // if(doubleKmer.direction == BACKWARD){
            //   printf("Reverse complement ");
            // }
            // printf("Kmer %s, extension of ", print_kmer(test_ext));
            // printf("%s, jchecks.\n", print_kmer(real));
            // printf("The real extension was %s\n", print_kmer(real_ext));
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
  printf("Finding next junction, starting at %s\n", print_kmer(doubleKmer->getKmer()));
  for (; doubleKmer->getDistToEnd() > 0; doubleKmer->forward())
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
  printf("Read: %s\n", &read[0]);
  int pos = 0;
  kmer_type kmer;
  DoubleKmer* firstJunc;
  DoubleKmer* doubleKmer = new DoubleKmer(&read);//stores current kmer throughout
  doubleKmer->forward();//don't scan the first position- points off the read
  DoubleKmer* lastJunc;
  
  //handle all junctions - link to each other
  while(find_next_junction(doubleKmer))
  {
    if(!junctionMap->contains(doubleKmer->getKmer())){//need to create a junction
      junctionMap->createJunction(doubleKmer->getKmer());
    }
    junctionMap->getJunction(doubleKmer->getKmer())->addCoverage(doubleKmer->getRealExtensionNuc());
    if(lastJunc){
      junctionMap->directLinkJunctions(lastJunc, doubleKmer);
    }
    else{ 
      firstJunc = new DoubleKmer(doubleKmer);
      Junction* junc = junctionMap->getJunction(firstJunc->getKmer());
      if(junc->nextJunc[firstJunc->getExtensionIndex(BACKWARD)] == -1){
          junc->update(firstJunc->getExtensionIndex(BACKWARD), firstJunc->getTotalPos(), -1);
      }
    }
    free(lastJunc);
    lastJunc = new DoubleKmer(doubleKmer);
    NbProcessed++;
  
    int dist = max(1,junctionMap->getSkipDist(doubleKmer, FORWARD));
    printf("Skip distance %d\n", junctionMap->getSkipDist(doubleKmer, FORWARD));
    doubleKmer->advanceDist(dist);
  }   
  free(doubleKmer);

  //handle case where there are no juncs
  if(!lastJunc){
      
  }
  //handle case with at least one junction
  else {
    Junction* junc = junctionMap->getJunction(lastJunc->getKmer());
    if(junc->nextJunc[lastJunc->getExtensionIndex(FORWARD)] == -1){
      junc->update(lastJunc->getExtensionIndex(FORWARD), lastJunc->getDistToEnd(), -1);
  }
    //remove any caps between the first and last junction
    // for(DoubleKmer uncapper = *firstJunc; uncapper.getTotalPos() <= lastJunc->getTotalPos(); uncapper.forward()){
    //   junctionMap->removeCap(uncapper.getKmer());
    // }

    //Handle leftmost junction cap
    //if(junctionMap->getSkipDist(firstJunc,BACKWARD) < firstJunc->getTotalPos()){
      //Cap* = findFirstCap(read, firstJunc->getTotalPos() - junctionMap->getSkipDist(firstJunc, BACKWARD));
      //if(!Cap){
        //PUSH CAP TO BEGINNING
      //}
      //else{
        //LINK JUNC TO FARTHER CAP
      //}
    //}
  }

  free(lastJunc);
  free(firstJunc);
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


