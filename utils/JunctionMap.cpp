#include "JunctionMap.h"
#include <string>
#include <fstream>
using std::ofstream;
using std::string;

//IMPLEMENTED
//NOT TESTED
//Gets the valid extension of the given kmer based on the bloom filter and cFPs.  Must JCheck! so this cuts off tips
//Assume its not a junction
//Returns -1 if there is no valid extension
//ASSUMES NO CFP SET- since this is only done in findSinks, BEFORE the cFPs are found
int JunctionMap::getValidJExtension(DoubleKmer kmer){
    //printf("Getting valid extension of %s\n", print_kmer(kmer.kmer));
    kmer_type nextKmer;
    for(int i = 0; i < 4; i++){
        nextKmer = kmer.getExtension(i, FORWARD);
        //printf("Testing extension %s\n", print_kmer(nextKmer));
        if(bloom->oldContains(get_canon(nextKmer))){
            if(jchecker->jcheck(nextKmer)){
                return i;
            }
        }
    }
    //No extension found!
    return -1;
}

//IMPLEMENTED
//NOT TESTED
//Scans forward from junction junc at index i with bloom filter
//If it hits another junction at or before the distance specified by the given junction, it returns null
//If it does not, it keeps scanning until it hits another junction or an actual sink
//If it hits a sink, it returns it.  If it hits a junction, it tests how far that junction points along the path.
//Based on the indicated overlap, it either decides the entire intermediate sequence is real or the connection is a 
//false positive connection
kmer_type * JunctionMap::findSink(Junction junc, kmer_type startKmer, int index){
    DoubleKmer doubleKmer(startKmer);
    kmer_type kmer;
    int scanDist;
    int maxDist = junc.dist[index];

    //get to the first kmer from which we need to bloom scan.  This is different for forwards and backward extensions
    if(index == 4){
        //in this case that's the reverse kmer 
        if(contains(doubleKmer.revcompKmer)){
            return NULL;
        }
        scanDist = 1;
        doubleKmer.reverse(); 
    }
    else{
        //in this case thats the forward version of the kmer extension specified by the index
        doubleKmer.forward(index);
        if(contains(doubleKmer.kmer)){
            return NULL;
        }
        scanDist = 2;
    }
    
    kmer_type nextJunc = -1; //stores the junction found
    int nextJuncDist, nextJuncExtIndex;
    int maxReadLength = 100; //the maximum distance we expect any junction to point- helps us know when to stop scanning
    int lastNuc =first_nucleotide(doubleKmer.revcompKmer); //stores the last nuc so we know which extension we came from

    kmer_type sinkKmer = -1;
    //if we're at the position where the sink would be, record the value for later use
    if(scanDist == maxDist){
        sinkKmer = doubleKmer.kmer;
    }

    //Scan forward until there's no chance of finding a junction that indicates an overlapping kmer 
    while(scanDist < maxDist + maxReadLength){ 
        //move forward if possible
        int validExtension = getValidJExtension(doubleKmer);
        if(validExtension == -1){
            if(sinkKmer == -1){
                return new kmer_type(doubleKmer.kmer);
            }
            return new kmer_type(sinkKmer);
        }
        doubleKmer.forward(validExtension); 
        scanDist++;

        //handle backward junction case
        if(contains(doubleKmer.revcompKmer)){
            nextJunc = doubleKmer.revcompKmer;
            nextJuncDist = scanDist;
            nextJuncExtIndex = lastNuc;
            break;
        }
        scanDist++;

        //handle forward junction case
        if(contains(doubleKmer.kmer)){
            nextJunc = doubleKmer.kmer;
            nextJuncDist = scanDist;
            nextJuncExtIndex = 4;
            break; 
        }

        //if we're at the position where the sink would be, record the value for later use
        if(scanDist == maxDist){
            sinkKmer = doubleKmer.kmer;
        }
        lastNuc = first_nucleotide(doubleKmer.revcompKmer);
    }

    //if no junction was found, must be a sink!
    if(nextJunc == -1){
        return new kmer_type(sinkKmer);
    }
    int otherMaxDist = getJunction(nextJunc)->dist[nextJuncExtIndex];
    int overlap = otherMaxDist + maxDist - nextJuncDist;
    if(overlap >= 0){
        return NULL; //if there is indicated overlap, this is not a sink.
    }
    return new kmer_type(sinkKmer); //if there is no indicated overlap between this and the next junction, we found a sink!
}

//IMPLEMENTED
//NOT TESTED
//Iterates through all of the junctions
//For each, bloom scans in every direction till another junction is hit or the stored distance is reached
//If there is no junction within that distance, a sink has been reached, and is added to the set of sinks.
//To be used FIRST after read scan
set<kmer_type>* JunctionMap::getSinks(){
    printf("Getting sinks.\n");
    set<kmer_type>* sinks = new set<kmer_type>({});
    kmer_type kmer;
    Junction junction;
    int juncsTested = 0;
    for(auto it = junctionMap.begin(); it != junctionMap.end(); it++, juncsTested++){
        kmer = it->first;
        junction = it->second;
        for(int i = 0; i < 5; i++){
            if( !junction.linked[i] && (i == 4 || junction.cov[i] > 0) ){
                kmer_type* sink = findSink(junction, kmer, i);
                if(sink){
                    sinks->insert(*sink);
                }
                else{
                }
            }
        }
        if(juncsTested % 10000 == 0) printf("Tested %d/%d junctions for sinks.\n", juncsTested, junctionMap.size());
    }
    return sinks;
}

//IMPLEMENTED
//NOT TESTED
//Iterates through all of the junctions, cleaning the non-complex ones.
//Every time a non-complex junction is found, the relevant cFPs are added to the cFP set and the junction is destroyed
//To be used AFTER findSinks
set<kmer_type>* JunctionMap::getCFPs(){
    set<kmer_type>* cFPs = new set<kmer_type>({});
    kmer_type kmer;
    Junction junction;
    int juncsTested = 0;
    int juncSize = junctionMap.size();
    for(auto it = junctionMap.begin(); it != junctionMap.end(); it++,  juncsTested++){
        kmer = it->first;
        junction = it->second;
        if (junction.numPathsOut() == 1){
            for(int i = 0; i < 4; i++){
                if(junction.cov[i] == 0){
                    kmer_type nextKmer = next_kmer(kmer, i, FORWARD);
                    if(bloom->oldContains(get_canon(nextKmer))){
                        if(jchecker->jcheck(nextKmer)){
                            cFPs->insert(nextKmer);
                        }
                    }
                }
            }
            killJunction(kmer);
        }
        if(juncsTested % 10000 == 0) printf("Tested %d/%d junctions for cFPs.\n", juncsTested, juncSize);
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

JunctionMap::JunctionMap(Bloom* bloo1, JChecker* jcheck){
  junctionMap = {};
  bloom = bloo1;
  jchecker = jcheck;
}