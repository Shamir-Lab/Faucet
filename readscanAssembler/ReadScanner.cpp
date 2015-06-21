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
  spacerDist = 30;//default.
  jchecker = checker;
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

//starts at position *pos, kmer *kmer on read, and if it returns true, pos and *kmer should be the
//pos and kmer of the next branch point. If it returns false, no guarantee- it ran off the end and we can handle that.
bool ReadScanner::find_next_junction(int* pos, kmer_type * kmer, string read){
  
  /*
  int lastSpacer = *pos;
  kmer_type potentialSpacer;
  int spacerExt;
  */

  for (; *pos < read.length()-sizeKmer; (*pos)++)//for each pos in read
  {
    kmer_type next_real = next_kmer_in_read(*kmer,*pos,&read[0], 0);//get next real kmer
  
    for(int nt=0; nt<4; nt++) {//for each extension
      kmer_type next_test = next_kmer(*kmer,nt, 0);//get possible extension
      if(next_real != next_test && next_real != *kmer && next_test != *kmer){//if the branch has 3 distinct kmers
        nextHash0 = bloom->roll_hash(hash0, NT2int(read[*pos]), nt, 0);
        nextHash1 = bloom->roll_hash(hash1, NT2int(read[*pos]), nt, 1);
        if(bloom->contains(nextHash0, nextHash1))//if the branch checks out initially
        { 
            if(!junctionMap->contains(*kmer)){//if it's nnot an already found junction
              NbJCheckKmer++;
              if(jchecker->jcheck(&read[*pos+1],nextHash0, nextHash1)){//if it j-checks- new junction!
                  junctionMap->createJunction(*kmer, NT2int(read[*pos+sizeKmer]));
                  return true;
              }
            }
            else{//if it is an already found junction
                return true;
            }
          }
      }
    }
    //handle potential spacers if there's no junction.
    
    /*if(*pos - lastSpacer == spacerDist){ //if we've scanned one spacer worth, this could be a spacer!
        potentialSpacer = *kmer;
        spacerExt = NT2int(read[*pos + sizeKmer]);
    }
    if(*pos - lastSpacer == 2*spacerDist){//make the current potential spacer, this is the new potential spacer
        junctionMap->createJunction(potentialSpacer, spacerExt);
        NbSpacers++;
        potentialSpacer = *kmer;
        spacerExt =NT2int(read[*pos + sizeKmer]);
        lastSpacer = *pos;
    }
    */

    shift_kmer(kmer, NT2int(read[*pos+sizeKmer]), 0); 
    bloom->advance_hash(&read[0], &hash0, &hash1, *pos, *pos+1);//advance hash values

    NbProcessed++;
    //if ((NbProcessed % 10000)==0) fprintf (stderr,"%c Deblooming kmers: %d",13,NbProcessed);
  }
  return false;
}

void ReadScanner::smart_traverse_read(string read){

  int pos = 0, lastJuncPos = 0, lastJuncExt = 0; 
  Junction * lastJunc = NULL;
  kmer_type kmer;
  bool noJuncs = true;
  getFirstKmerFromRead(&kmer,&read[0]);//stores current kmer throughout
  hash0 = bloom->get_rolling_hash(kmer, 0);
  hash1 = bloom->get_rolling_hash(kmer, 1);

  while(pos < read.length()-sizeKmer)
  {
      //assumption: kmer is set, pos is set, lastJunc info is all set
      juncInfo = junctionMap->getJunction(kmer);
      if(juncInfo) // is a seen junction!
      { 

        juncInfo->update(NT2int(read[pos+sizeKmer]), 1, 0);

        if(lastJunc)//update info on last junction if it exists
        {
          lastJunc->update(lastJuncExt, pos - lastJuncPos, 0);
          juncInfo->update(0, 0, pos - lastJuncPos);//doesn't actually make sense
        } 

        lastJunc = juncInfo; //this is now the last junc
        lastJuncExt = NT2int(read[pos+sizeKmer]); //the extension is at pos+k in the read
        lastJuncPos = pos; //this is now the last junc position
       
        advance_kmer(&read[0], &kmer, pos, pos + (int)lastJunc->ext[lastJuncExt]);//advance kmer to jump forward
        bloom->advance_hash(&read[0], &hash0, &hash1, pos, pos+(int)lastJunc->ext[lastJuncExt]);//advance hash values

        pos += (int)lastJunc->ext[lastJuncExt]; //set new pos appropriately 
        NbProcessed++;
      }   
      else// not at a seen junction, need to scan! 
      {         
        if(!find_next_junction(&pos, &kmer, read)) {  //scanned off the end of the read
          pos = read.length()-sizeKmer; //pos is now = length so we'll be done after this
          if(lastJunc){
            lastJunc->update(lastJuncExt, pos-lastJuncPos, 0);
          }
          continue;
        }
        noJuncs = false;//found a junction, so info was updated.  Proceed.
      }
  }
  NbNoJuncs += noJuncs; //add one t
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
    smart_traverse_read(read);
    revcomp_sequence(&read[0], read.length());
    smart_traverse_read(read);
    if ((readsProcessed%10000)==0) fprintf (stderr,"Reads processed: %c %lld",13,(long long)readsProcessed);
    readsProcessed++;
  }

  solidReads.close();
  time(&stop);
  printf("Time in seconds for read scan: %f \n", difftime(stop,start));
}


