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

//returns the junction if it exists or a null pointer otherwise
Junction* JunctionMap::getJunction(DoubleKmer kmer){
  return getJunction(kmer.getKmer());
}

void JunctionMap::directLinkJunctions(DoubleKmer* kmer1, DoubleKmer* kmer2){
    int ext1 = kmer1->getExtensionIndex(FORWARD);
    int ext2 = kmer2->getExtensionIndex(BACKWARD);
    Junction* junc1 = getJunction(kmer1);
    Junction* junc2 = getJunction(kmer2);
    
    int dist = kmer2->getTotalPos() - kmer1->getTotalPos();

    junc1->update(ext1, dist);
    junc2->update(ext2, dist);
}

//void JunctionMap::extendJunctionCap(DoubleKmer* kmer, bool dir){
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
//}

void JunctionMap::finishGraphWriteBasicContigs(){

}

void JunctionMap::createJunction(DoubleKmer* doubleKmer){  
  createJunction(doubleKmer->getKmer());
}

//creates a new junction with first found extension nucExt
//returns it using the pointer supplied.
void JunctionMap::createJunction(kmer_type kmer){  
  junctionMap.insert(std::pair<kmer_type, Junction>(kmer, *(new Junction())));
}


void JunctionMap::killJunction(kmer_type kmer){
  junctionMap.erase(kmer);
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

bool JunctionMap::contains(DoubleKmer* doubleKmer){
  return contains(doubleKmer->getKmer());
}

bool JunctionMap::contains(kmer_type kmer){
    return junctionMap.find(kmer) != junctionMap.end();
}

JunctionMap::JunctionMap(Bloom* bloo1){
  junctionMap = {};
  bloom = bloo1;
}