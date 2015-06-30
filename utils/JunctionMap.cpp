#include "JunctionMap.h"
#include <string>
#include <fstream>
using std::ofstream;
using std::string;

int JunctionMap::getNumComplexJunctions(){
  int count = 0;
  for(auto juncIt = junctionMap.begin(); juncIt != junctionMap.end(); juncIt++){
     if(juncIt->second.numPathsOut() > 1){
        count++;
     }  
  }
  return count;
}

//returns the junction if it exists or a null pointer otherwise
Junction* JunctionMap::getJunction(kmer_type kmer){
  auto juncIt = junctionMap.find(kmer);
  if(juncIt == junctionMap.end()){
    return NULL;
  }
  else{
    return &juncIt->second;
  }
}

void JunctionMap::linkJunctions(kmer_type kmer1, int ext1, kmer_type kmer2, int ext2, int dist){
   getJunction(kmer1)->update(ext1, dist, kmer2);
   getJunction(kmer2)->update(ext2, dist, kmer1);
}

//returns the junction if it exists or a null pointer otherwise
Cap* JunctionMap::getCap(kmer_type kmer){
  auto capIt = capMap.find(kmer);
  if(capIt == capMap.end()){
    return NULL;
  }
  else{
    return &capIt->second;
  }
}

//creates a new junction with first found extension nucExt
void JunctionMap::createJunction(kmer_type kmer){  
    junctionMap.insert(
      std::pair<kmer_type, Junction>(kmer, *(new Junction()))); 
}

void JunctionMap::addCap(kmer_type kmer, Cap cap){
  capMap.insert(std::pair<kmer_type, Cap>(kmer, cap));
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