#include "JunctionMap.h"
#include <string>
#include <fstream>
using std::ofstream;
using std::string;

int JunctionMap::getNumComplexJunctions(){
  int count = 0;
  for(auto juncIt = junctionMap.begin(); juncIt != junctionMap.end(); juncIt++){
     if(juncIt->second.numPathsOut() != 1){
        count++;
     }  
  }
  return count;
}

int JunctionMap::getSkipDist(DoubleKmer* doubleKmer, bool direction){
    int index = doubleKmer->getExtensionIndex(direction);
    return getJunction(doubleKmer->getKmer())->dist[index];
}

//returns the junction if it exists or a null pointer otherwise
Junction* JunctionMap::getJunction(kmer_type kmer){
  auto juncIt = junctionMap.find(kmer);
  if(juncIt == junctionMap.end()){
    return NULL;
  }
  else{
    return &(juncIt->second);
  }
}

void JunctionMap::directLinkJunctions(DoubleKmer* kmer1, DoubleKmer* kmer2){
    int ext1 = kmer1->getExtensionIndex(FORWARD);
    int ext2 = kmer2->getExtensionIndex(BACKWARD);
    
    int dist = kmer2->getTotalPos() - kmer1->getTotalPos();

    getJunction(kmer1->getKmer())->update(ext1, dist, kmer2->getKmer());
    getJunction(kmer2->getKmer())->update(ext2, dist, kmer1->getKmer());
}

void JunctionMap::extendJunctionCap(DoubleKmer* kmer, bool dir){
    // int juncExt = kmer->getExtensionIndex(dir);
    // Junction* junc = getJunction(kmer->getKmer());
    // removeCap(junc->nextJunc[juncExt]);
    
    // kmer_type capKmer;
    // int newDist;
    // string read = *(kmer->read);
    // if(dir == BACKWARD){
    //   getFirstKmerFromRead(&capKmer, &read[0]);
    //   capKmer = revcomp(capKmer);
    //   newDist = kmer->getTotalPos();
    // }
    // else{
    //   getFirstKmerFromRead(&capKmer,&read[read.length() - sizeKmer]);
    //   newDist = 2*read.length()-kmer->getTotalPos()-1;
    // }
    // junc->nextJunc[juncExt] = capKmer;
    // junc->dist[juncExt] = newDist;

    // addCap(capKmer, new Cap(newDist, kmer->getKmer()));
}

void superLinkJunctionToCap(DoubleKmer* kmer, Cap* cap, bool dir){

}

void superLinkCapToCap(Cap* leftCap, Cap* rightCap){

}
  
//creates a new junction with first found extension nucExt
//returns it using the pointer supplied.
void JunctionMap::createJunction(kmer_type kmer){  
  junctionMap.insert(std::pair<kmer_type, Junction>(kmer, *(new Junction())));
}

void JunctionMap::writeToFile(string filename){
    ofstream jFile;
    jFile.open(filename);

    printf("There are %d complex junctions.\n", getNumComplexJunctions());
    printf("Writing to junction file\n");
    kmer_type kmer;
    for(auto juncIt = junctionMap.begin(); juncIt != junctionMap.end(); juncIt++){
        kmer = juncIt->first;
        jFile << print_kmer(kmer) << " " ;
        juncIt->second.writeToFile(&jFile);
        jFile << "\n";    
    }
    printf("Done writing to junction file\n");
    jFile.close();
}

int JunctionMap::getNumJunctions(){
    return junctionMap.size();
}

bool JunctionMap::contains(kmer_type kmer){
    return junctionMap.find(kmer) != junctionMap.end();
}

JunctionMap::JunctionMap(){
  junctionMap = {};
  capMap = {};
}