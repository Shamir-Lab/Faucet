#include "ReadScanner.h"
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <time.h>
using namespace std;

ReadScanner::ReadScanner(string readFile, Bloom* bloo1){
  reads_file = readFile;
  bloom = bloo1;
  last = new kmer_type[20000];
  nextList = new kmer_type[20000];
  j = 0; //default no j-checking
}

void ReadScanner::printScanSummary(){
  printf("\n Distinct junctions: %d \n", junctionMap.size());
  printf(" Number of kmers that we j-checked: %lli \n", NbJCheckKmer);
  printf (" Number of reads with no junctions: %lli \n",NbNoJuncs);
  printf("Number of skipped kmers: %lli \n", ((readLength-sizeKmer)*readsProcessed*2) - NbProcessed);
  printf("Number of processed kmers: %lli \n", NbProcessed);
}

void ReadScanner::setJ(int jVal){
  j = jVal;
}

//creates a new junction with first found extension nucExt
void ReadScanner::createJunction(kmer_type kmer, int nucExt){
    junction* junc = new junction;
    updateJunction(junc, nucExt, 1, 0);                  
    junctionMap[kmer] = *junc;
}

//Updates the junc info based on finding a path of length length from the extension nucExt
void ReadScanner::updateJunction(junction * junc, int nucExt, int lengthFor, int lengthBack){
      junc->ext[nucExt] = max(junc->ext[nucExt],(unsigned char) lengthFor);
      junc->back = max(junc->back, (unsigned char) lengthBack);
}

//j = 0 always returns true
//j > 0 checks extensions up to j deep from kmer, and returns true if there is a sequence of j extensions
//which returns all positive in the bloom filter
inline bool ReadScanner::jcheck(kmer_type kmer, int strand){
  if(j == 0){
    return true;
  }
  kmer_type this_kmer, nextKmer;
  int lastCount, nextCount;
  lastCount = 1;
  last[0] = kmer;

  for(int i = 0; i < j; i++){//for each level up to j
    nextCount = 0; //have found no extensions yet
    for(int k = 0; k < lastCount; k++){ //for each kmer in the last level
      this_kmer = last[k]; //set the working kmer
      for(int nt = 0; nt < 4; nt++){ //for each possible extension
        nextKmer = next_kmer(this_kmer, nt, strand);//set the extension kmer
        
        if(bloom->contains(get_canon(nextKmer))){ //add to next level if it's in the bloom filter
          if(i == (j-1)){
            return true;//if this is the last level return true after the first check
          }
          nextList[nextCount] = nextKmer;
          nextCount++;
        }
      }
    }

    if(nextCount == 0){ //if there are no kmers in the list now, return false
      return false;
    }
    //reset counts and lists for next level of th search
    lastCount = nextCount;
    temp = last; 
    last = nextList;
    nextList = temp;
  }
}


//starts at position *pos, kmer *kmer on read, and if it returns true, pos and *kmer should be the
//pos and kmer of the next branch point. If it returns false, no guarantee- it ran off the end and we can handle that.
bool ReadScanner::find_next_junction(int* pos, kmer_type * kmer, string read){
  
  for (; *pos < read.length()-sizeKmer; *pos++)//for each pos in read
  {
    kmer_type next_real = next_kmer_in_read(*kmer,*pos,&read[0], 0);//get next real kmer
  
    for(int nt=0; nt<4; nt++) {//for each extension
      kmer_type next_test = next_kmer(*kmer,nt, 0);//get possible extension
      if(next_real != next_test && next_real != *kmer && next_test != *kmer){//if the branch has 3 distinct kmers
        if(bloom->contains(get_canon(next_test)))//if the branch checks out initially
        { 
            juncIt = junctionMap.find(*kmer);
            if(juncIt == junctionMap.end()){//if it's nnot an already found junction
              NbJCheckKmer++;
              if(jcheck(next_test, 0)){//if it j-checks- new junction!
                  createJunction(*kmer, NT2int(read[*pos+sizeKmer]));
                  return true;
              }
            }
            else{//if it is an already found junction
                return true;
            }
          }
      }
    }
    shift_kmer(kmer, NT2int(read[*pos+sizeKmer]), 0);
    NbProcessed++;
    //if ((NbProcessed % 10000)==0) fprintf (stderr,"%c Deblooming kmers: %d",13,NbProcessed);
  }
  return false;
}

void ReadScanner::smart_traverse_read(string read){

  int pos = 0, lastJuncPos = 0, lastJuncExt = 0; 
  junction * lastJunc = NULL;
  kmer_type kmer;
  bool noJuncs = true;
  getFirstKmerFromRead(&kmer,&read[0]);//stores current kmer throughout

  while(pos < read.length()-sizeKmer)
  {
      //assumption: kmer is set, pos is set, lastJunc info is all set
      juncIt = junctionMap.find(kmer);
      if(juncIt != junctionMap.end()) // is a seen junction!
      { 
        juncInfo = &juncIt->second;

        updateJunction(juncInfo, NT2int(read[pos+sizeKmer]), 1, 0);

        if(lastJunc)//update info on last junction if it exists
        {
          updateJunction(lastJunc, lastJuncExt, pos - lastJuncPos, 0);
          updateJunction(juncInfo, 0, 0, pos - lastJuncPos);//doesn't actually make sense
        } 

        lastJunc = juncInfo; //this is now the last junc
        lastJuncExt = NT2int(read[pos+sizeKmer]); //the extension is at pos+k in the read
        lastJuncPos = pos; //this is now the last junc position
       
        advance_kmer(&read[0], &kmer, pos, pos + (int)lastJunc->ext[lastJuncExt]);//advance kmer to jump forward
        pos += (int)lastJunc->ext[lastJuncExt]; //set new pos appropriately 
        NbProcessed++;
      }   
      else// not at a seen junction, need to scan! 
      {         
        if(!find_next_junction(&pos, &kmer, read)) {  //scanned off the end of the read
          pos = read.length()-sizeKmer; //pos is now = length so we'll be done after this
          if(lastJunc){
            updateJunction(lastJunc, lastJuncExt, pos-lastJuncPos, 0);
          }
          continue;
        }
        noJuncs = false;//found a junction, so info was updated.  Proceed.
      }
  }
  NbNoJuncs += noJuncs; //add one t
}

void ReadScanner::scanReads(int genome_size)
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
  printf("Weight before read scan: %li \n", bloom->weight());
  int lastSum = 0, thisSum = 0;
  while (getline(solidReads, read))
  {
    //lastSum = thisSum;
    readLength = read.length();
    smart_traverse_read(read);
    revcomp_sequence(&read[0], readLength);
    smart_traverse_read(read);
    if ((readsProcessed%10000)==0) fprintf (stderr,"Reads processed: %c %lld",13,(long long)readsProcessed);
    readsProcessed++;
  }

  solidReads.close();
  time(&stop);
  printf("Time in seconds for read scan: %f \n", difftime(stop,start));
}

void ReadScanner::writeJunction(ofstream*jFile, junction toPrint){
  for(int i = 0; i < 4; i++){
    *jFile << (int)toPrint.ext[i] << " " ;
  }
  *jFile << (int)toPrint.back << " ";
  *jFile << "\n";
}

void ReadScanner::junctionMapToFile(string filename){
    ofstream jFile;
    jFile.open(filename);
    junction toPrint;

    printf("Writing to junction file\n");
    kmer_type kmer;
    for(juncIt = junctionMap.begin(); juncIt != junctionMap.end(); juncIt++){
        kmer = juncIt->first;
        jFile << print_kmer(kmer) << " " ;
        writeJunction(&jFile, juncIt->second);    
    }
    printf("Done writing to junction file\n");
    jFile.close();
}