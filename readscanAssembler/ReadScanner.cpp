#include "ReadScanner.h"
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
using namespace std;


ReadScanner::ReadScanner(JunctionMap* juncMap, string readFile, Bloom* bloo1, Bloom* short_filter, Bloom* long_filter, JChecker* checker, int maxDist){
  reads_file = readFile;
  bloom = bloo1;
  short_pair_filter = short_filter;
  long_pair_filter = long_filter;
  junctionMap = juncMap;
  jchecker = checker;
  maxSpacerDist = maxDist;
}

void ReadScanner::printScanSummary(){
  printf("\nDistinct junctions: %lli \n", (uint64_t)junctionMap->getNumJunctions());
  //printf("Number of junction pairs that exist on reads: %d\n", juncPairSet.size());
  printf("Number of kmers that we j-checked: %lli \n", NbJCheckKmer);
  printf("Number of reads with no junctions: %lli \n",NbNoJuncs);
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
  // if (real_ext == revcomp(real)){ return true; }

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
bool ReadScanner::find_next_junction(ReadKmer * readKmer, int lastJuncPos){

  //Iterate forward to the end of the read
  for (; readKmer->getDistToEnd() > 2*jchecker->j; readKmer->forward()) //CHANGED TO 2*j from 0
  {
      //check for an already found junciton
      if(junctionMap->isJunction(readKmer->getKmer())){
        // std::cout << "known junction found\n";
        return true;
      }
      //check for the max spacer dist
      if(readKmer->getTotalPos() - lastJuncPos >= 2*maxSpacerDist-1 ){ //|| testForJunction(*readKmer)){
        
        // if (junctionMap->isJunction(readKmer->getRevCompKmer())){ continue; } // avoid creating junction if RC already is one
        return true;
      }
      // check if new junction should be created
      if(testForJunction(*readKmer)){
        // std::cout << "added "<< readKmer->directionAsString() <<" junction with index " << readKmer->getTotalPos() << std::endl;
        return true;
      }
      NbProcessed++;
      // delete rcReadKmer;
  }
  return false;
}

//Should only be called on a read with no real junctions
//Adds a fake junction in the middle and points it to the two ends.  This ensures we have coverage of long linear regions, and that we capture
//sinks at the end of such regions.
//Returns the real extension of the fake junction, for junction pairing
kmer_type ReadScanner::add_fake_junction(string read){
  // std::cout << "added fake junction\n";
  ReadKmer* middleKmer = new ReadKmer(&read, read.length()/2- sizeKmer/2, FORWARD);
  kmer_type extension = middleKmer->getRealExtension();
  junctionMap->createJunction(middleKmer);
  Junction* junc = junctionMap->getJunction(middleKmer);
  junc->addCoverage(middleKmer->getRealExtensionNuc());
  junc->update(middleKmer->getExtensionIndex(BACKWARD), middleKmer->getTotalPos()-2*jchecker->j);
  junc->update(middleKmer->getExtensionIndex(FORWARD), middleKmer->getDistToEnd()-2*jchecker->j);
  delete middleKmer;
  middleKmer = nullptr;
  return extension;
}


//Scans a read. 
//Identifies all junctions on the read, and links adjacent junctions to each other.
//Also updates the relevant distance field on the first junction to point to the start of the read, and on the last
//Junction to point to the end of the read.
//If there are no junctions, add_fake_junction is called
std::list<kmer_type> ReadScanner::scan_forward(string read, bool no_cleaning){
  std::list<kmer_type> result = {};
  ReadKmer* readKmer = new ReadKmer(&read);//stores current kmer throughout

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
  int rev_pos = 0;
  int for_pos = 0;
  int lastJuncPos = 0;

  //handle all junctions - scan forward, find and create junctions, and link adjacent junctions to each other
  while(find_next_junction(readKmer, lastJuncPos))
  { 
     junc = junctionMap->getJunction(readKmer);
     lastJuncPos = readKmer->getTotalPos(); //need to record this before we advance any farther

    //create a junction at the current spot if none exists
    if(!junc){
      junctionMap->createJunction(readKmer);
      junc = junctionMap->getJunction(readKmer);
    }
    // else{
    //   std::cout << "found real junction\n";
      // std::cout << print_kmer(readKmer->doubleKmer.getCanon()) << std::endl;
    // }
    result.push_back(readKmer->getRealExtension());

    // mark first R junction and last F junction
    if(readKmer->direction == BACKWARD){
      if(!backJunc){
        backJunc = new ReadKmer(readKmer);
        firstBackJunc = new ReadKmer(readKmer);
        rev_pos = firstBackJunc->pos;
      }
      else{
        *backJunc = readKmer;
      }
    }
    else{
      // if(backJunc){
      //    short_pair_filter->addPair(JuncPair(backJunc->getRealExtension(), readKmer->getRealExtension()));
      // }
      if(!lastForwardJunc){
        lastForwardJunc = new ReadKmer(readKmer);
        for_pos = lastForwardJunc-> pos;
      }
      else{
        *lastForwardJunc = readKmer;
      }
    }

    junc->addCoverage(readKmer->getRealExtensionNuc()); //add coverage of the junction
    
    //if there was a last junction, link the two 
    if(lastKmer){
      junctionMap->directLinkJunctions(lastKmer, readKmer, lastJunc, junc);
    }
    //If this is the first junction, link it to the beginning of the read.
    else{ 
      // std::cout<< "first junc position is " << readKmer->getTotalPos()-2*jchecker->j <<std::endl;
      lastKmer = new ReadKmer(readKmer);
      junc
        ->update(readKmer->getExtensionIndex(BACKWARD), readKmer->getTotalPos()-2*jchecker->j);//-2*j ADDED
    }

    *lastKmer = *readKmer;
    lastJunc = junc;

    int index = readKmer->getExtensionIndex(FORWARD);
    int dist = max(1, (int)(junc->dist[index]));
    // std::cout << dist << std::endl;
    readKmer->advanceDist(dist); 

    NbProcessed++,  NbSkipped += dist-1;
  }   
  
  //If there were no junctions on the read, add a fake junction
  if(!lastKmer){
      NbNoJuncs++;
      result.push_back(add_fake_junction(read));
      // std::cout << "added fake junction" << std::endl;
  }
  else {
    //If there was at least one junction, point the last junction found to the end of the read
    lastJunc->update(lastKmer->getExtensionIndex(FORWARD), lastKmer->getDistToEnd()-2*jchecker->j); //2*j ADDED
    // std::cout << "some junction exists, connecting with "<< lastKmer->getDistToEnd()-2*jchecker->j << std::endl;
  }
  if (!no_cleaning){
    if(result.size()==2){
      if (firstBackJunc && lastForwardJunc && !(rev_pos > for_pos)){
        // printf("accepted rev_pos %d for_pos %d\n", rev_pos, for_pos);
        short_pair_filter->addPair(JuncPair(firstBackJunc->getRealExtension(), lastForwardJunc->getRealExtension()));
      }
    }
    else if(result.size()>2){ 
      // copy list to vector to be able to iterate over - not sure if this is optimal
      std::vector<kmer_type> v{ std::begin(result), std::end(result) };
      for (int i = 0; i< v.size()-2; i++){
        short_pair_filter->addPair(JuncPair(v[i], v[i+2]));            
      }
      v.clear();
    }
  }
  delete lastKmer;
  delete readKmer;
  delete firstBackJunc;
  delete lastForwardJunc;
  return result;
}

std::list<std::pair<string, int> > ReadScanner::getValidReads(string read){
  std::list<std::pair<string,int> > result = {};
  //Move to the first valid kmer
  int start = 0, end = 0;
  int minLength = sizeKmer;//only use valid reads with at least this many valid kmers

  for(ReadKmer kmer = ReadKmer(&read); kmer.getDistToEnd()>= 0; kmer.forward(), kmer.forward()){
    if(bloom->oldContains(kmer.getCanon())){
      end++;
    } 
    else{
      if(end >= start + minLength){ //buffer to ensure no reads exactly kmer size- might be weird edge cases there
        result.push_back(std::make_pair(read.substr(start, end-start + sizeKmer-1) , start));
        // std::cout << "Adding " << read.substr(start, end-start + sizeKmer) << " as valid string\n";
      }
      start = kmer.pos + 1;
      end = kmer.pos + 1;
    }
  }
   if(end >= start + minLength){ //buffer to ensure no reads exactly kmer size- might be weird edge cases there
        result.push_back(std::make_pair(read.substr(start, end-start+sizeKmer-1) , start));
        // std::cout << "Adding " << read.substr(start, end-start+sizeKmer) << " as valid string\n";
    }
  return result;
}

//returns back junctions for use in paired end info
std::list<kmer_type> ReadScanner::scanInputRead(string read, bool no_cleaning){
  std::list<kmer_type> result = {};
  std::list<string> readList = getUnambiguousReads(read);
  std::list<std::pair<string, int> > validReads;
  string validRead;
  bool mercy = true;
  while(!readList.empty()){
        read = readList.front();
        readList.pop_front();
        if(read.length() >= sizeKmer + 2*jchecker->j + 1){
          unambiguousReads++;
          validReads = getValidReads(read);

          // copy the list over to a vector, to make interation a bit simpler
          // std::vector<std::pair<string, int > > vr_copy{ std::begin(validReads), std::end(validReads) };

          // // retrieval of mercy kmers - when we have two adjacent solid regions
          // // seperated by a low coverage region, and no junctions flanking the 
          // // low region, we add the kmers of the low region to B2          
          // if (mercy && validReads.size()>=2){
          //   for (int i = 0; i < vr_copy.size() - 1; i++){
          //     int start = vr_copy[i].second + vr_copy[i].first.length();
          //     int end = vr_copy[i+1].second;
          //     kmer_type frontKmer;
          //     kmer_type backKmer;
              
          //   }
          // }

          while(!validReads.empty()){
            validRead = validReads.front().first;
            validReads.pop_front();
            result.splice(result.end(),scan_forward(validRead, no_cleaning));
            readsNoErrors++;
          }
        }
    }
            
    return result;
}

void ReadScanner::scanReads(bool fastq, bool paired_ends, bool no_cleaning)
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
  int emptyCount = 0, notEmptyCount = 0;
  while (getline(solidReads, read))
  {
    getline(solidReads, read);//since it's a fasta we skip the first of every pair of lines
    // std::cout << "read" << std::endl;
    if(firstEnd){
      backJuncs1 = scanInputRead(read, no_cleaning);
    }
    else{
      backJuncs2 = scanInputRead(read, no_cleaning);
    }
    
    if(paired_ends && !firstEnd){//Current logic: ensure each one in the first list has a pair in the second
      if(!backJuncs1.empty() && !backJuncs2.empty()){
        notEmptyCount++;
        //printf("Backjuncs1 and backjunc2 both not empty.\n");
        for(auto it = backJuncs1.begin(); it != backJuncs1.end(); it++){
            //printf("Processing one back junc 1\n");
            kmer_type pair1 = *it;
            bool paired = false;
            if (!no_cleaning){
              for(auto it2 = backJuncs2.begin(); it2 != backJuncs2.end(); it2++){
                kmer_type pair2 = *it2;
                 //printf("Adding pair\n");
                 //pair_filter->addPair(JuncPair(pair1, pair2));
              
                if(long_pair_filter->containsPair(JuncPair(pair1,pair2))){
                    paired = true;
                }
              }
              if (!paired){
                long_pair_filter->addPair(JuncPair(pair1, *backJuncs2.begin()));
                // long_pair_filter->addPair(JuncPair(pair1, backJuncs2.back()));
              }
          }
        }
      }
      else{emptyCount++;}
    }
    if ((readsProcessed%100000)==0) fprintf (stdout," reads scanned: %c %lld\n",13,(long long)readsProcessed);
    readsProcessed++;
    if(fastq) getline(solidReads, read), getline(solidReads, read); //if fastq skip two more lines
    firstEnd = !firstEnd;
  }
  printf("Empty count: %d, not empty count: %d\n", emptyCount, notEmptyCount);

  solidReads.close();
  time(&stop);
  printf("Reads processed: %lli\n", readsProcessed);
  printf("Unambiguous reads: %lli\n", unambiguousReads);
  printf("Time in seconds for read scan: %f \n", difftime(stop,start));
}

