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

//creates a new junction with first found extension nucExt
void JunctionMap::createJunction(kmer_type kmer, int nucExt){  
    junctionMap.insert(
      std::pair<kmer_type, Junction>(kmer, *(new Junction(nucExt)))); 
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
}
