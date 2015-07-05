#include "JunctionMap.h"
#include <string>
#include <fstream>
using std::ofstream;
using std::string;

//NOT IMPLEMENTED
//NOT TESTED
//Scans forward from junction junc at index i with bloom filter
//If it hits another junction at or before the distance specified by the given junction, it returns null
//If it does not, it returns the kmer at that distance, since that must be a sink.
kmer_type * JunctionMap::findSink(Junction junc, kmer_type startKmer, int index){
    DoubleKmer doubleKmer(startKmer);
    kmer_type kmer;
    int scanDist = 0;
    int maxDist = junc.dist[index];

    //get to the first kmer that has a backward and forward kmer not equal to the start point
    if(index == 4){

    }
    else{

    }

    //then scan!
    while(scanDist < maxDist){
        if(contains(doubleKmer.revcompKmer)){
            return NULL;
        }
        scanDist++;
        if(scanDist > maxDist){
            return new kmer_type(doubleKmer.kmer);
        }
        kmer = doubleKmer.kmer;
        if(contains(kmer)){
            //return 
        }
    }

}

//IMPLEMENTED
//NOT TESTED
//Iterates through all of the junctions
//For each, bloom scans in every direction till another junction is hit or the stored distance is reached
//If there is no junction within that distance, a sink has been reached, and is added to the set of sinks.
set<kmer_type> JunctionMap::getSinks(){
    set<kmer_type> sinks = {};
    kmer_type kmer;
    Junction junction;
    for(auto it = junctionMap.begin(); it != junctionMap.end(); it++){
        kmer = it->first;
        junction = it->second;
        for(int i = 0; i < 4; i++){
            if(!junction.linked[i]){
                kmer_type* sink = findSink(junction, kmer, i);
                if(sink){
                    sinks.insert(*sink);
                }
            }
        }
    }
    return sinks;
}

//IMPLEMENTED
//NOT TESTED
//Iterates through all of the junctions, cleaning the non-complex ones.
//Every time a non-complex junction is found, the relevant cFPs are added to the cFP set and the junction is destroyed
set<kmer_type> JunctionMap::getCFPs(){
    set<kmer_type> cFPs = {};
    kmer_type kmer;
    Junction junction;
    for(auto it = junctionMap.begin(); it != junctionMap.end(); it++){
        kmer = it->first;
        junction = it->second;
        if (junction.numPathsOut() == 1){
            for(int i = 0; i < 4; i++){
                if(junction.cov[i] == 0){
                    kmer_type nextKmer = next_kmer(kmer, i, FORWARD);
                    if(bloom->oldContains(nextKmer)){
                        cFPs.insert(nextKmer);
                    }
                }
            }
            killJunction(kmer);
        }
    }
    return cFPs;
}

int JunctionMap::getNumComplexJunctions(){
  int count = 0;
  for(auto juncIt = junctionMap.begin(); juncIt != junctionMap.end(); juncIt++){
     if(juncIt->second.numPathsOut() != 1){
        count++;
     }  
  }
  return count;
}

int JunctionMap::getSkipDist(ReadKmer* readKmer, bool direction){
    int index = readKmer->getExtensionIndex(direction);
    return getJunction(readKmer->getKmer())->dist[index];
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
Junction* JunctionMap::getJunction(ReadKmer kmer){
  return getJunction(kmer.getKmer());
}

void JunctionMap::directLinkJunctions(ReadKmer* kmer1, ReadKmer* kmer2){
    int ext1 = kmer1->getExtensionIndex(FORWARD);
    int ext2 = kmer2->getExtensionIndex(BACKWARD);
    Junction* junc1 = getJunction(kmer1);
    Junction* junc2 = getJunction(kmer2);
    
    int dist = kmer2->getTotalPos() - kmer1->getTotalPos();

    junc1->update(ext1, dist);
    junc2->update(ext2, dist);
    junc1->link(ext1);
    junc2->link(ext2);
}

void JunctionMap::finishGraphWriteBasicContigs(){

}

void JunctionMap::createJunction(ReadKmer* readKmer){  
  createJunction(readKmer->getKmer());
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

bool JunctionMap::contains(ReadKmer* readKmer){
  return contains(readKmer->getKmer());
}

bool JunctionMap::contains(kmer_type kmer){
    return junctionMap.find(kmer) != junctionMap.end();
}

JunctionMap::JunctionMap(Bloom* bloo1){
  junctionMap = {};
  bloom = bloo1;
}