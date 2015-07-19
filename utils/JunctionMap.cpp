#include "JunctionMap.h"
#include <string>
#include <fstream>
using std::ofstream;
using std::string;

//Gets the valid extension of the given kmer based on the bloom filter and cFPs.  Uses JChecking! so this cuts off tips
//Assume the given kmer is not a junction
//Returns -1 if there is no valid extension
//Returns -2 if there are multiple
//ASSUMES NO CFP SET- since this is only done in findSinks, BEFORE the cFPs are found
int JunctionMap::getValidJExtension(DoubleKmer kmer){
    kmer_type nextKmer;
    int answer = -1;
    for(int i = 0; i < 4; i++){
        nextKmer = kmer.getExtension(i, FORWARD);
        if(bloom->oldContains(get_canon(nextKmer))){
            if(jchecker->jcheck(nextKmer)){
                if(answer != -1){
                    //Found multiple valid extensions!
                    return -2;
                }
                answer = i; 
            }
        }
    }
    return answer;
}

//Scans forward from junction junc at index i with bloom filter
//If it hits another junction at or before the distance specified by the given junction, it returns null
//If it does not, it keeps scanning until it hits another junction or an actual sink
//If it hits a sink, it returns it.  If it hits a junction, it tests how far that junction points along the path.
//Based on the indicated overlap, it either decides the entire intermediate sequence is real or the connection is a 
//false positive connection
//Assumes that the junction at the given index is not linked.  If its linked there's no need to do a search.
kmer_type * JunctionMap::findSink(Junction junc, kmer_type startKmer, int index){
    DoubleKmer doubleKmer(startKmer);
    kmer_type kmer;
    int scanDist;
    int maxDist = junc.dist[index];

    //First, process the first 1-2 kmers in order to reach the first kmer from which can properly bloom scan.  
    //This is different for forwards and backward extensions.
    if(index == 4){
        //If we're searching backwards, we only need to specially process the reverse kmer, and then scan from there
        if(isJunction(doubleKmer.revcompKmer)){
            return NULL;
        }
        scanDist = 1;
        doubleKmer.reverse(); 
    }
    else{
        //if we're searching forwards, we need to specially handle the first extension since the junction tells us where to go.
        //From there, we can scan normally.
        doubleKmer.forward(index);
        if(isJunction(doubleKmer.kmer)){
            return NULL;
        }
        scanDist = 2;
    }
    
    //After the scan loop, these will store the needed info about the junction found
    kmer_type nextJunc = -1; 
    int nextJuncDist, nextJuncExtIndex;
    int lastNuc; //stores the last nuc so we know which extension we came from. 

    kmer_type sinkKmer;
    //if we're at or past the position where the sink would be, record the value for later use
    if(scanDist >= maxDist ){ //REMOVED THE - 2 * jchecker->j
        if(scanDist > maxDist){
            printf("Error: scanDist %d is greater than maxDist %d.\n", scanDist, maxDist);
        }
        sinkKmer = doubleKmer.kmer;
    }

    //Scan forward until there's no chance of finding a junction that indicates an overlapping kmer 
    while(scanDist < maxDist + maxReadLength){ 

        //move forward if possible
        int validExtension = getValidJExtension(doubleKmer);
        if(validExtension == -1){//Then we're at the end of the line! Found a sink.
            return new kmer_type(sinkKmer);
        }
        if(validExtension == -2){ //this ambiguity only happens when a sink has two false extensions- but it's still a sink!
            return new kmer_type(doubleKmer.kmer);
        }
        lastNuc = first_nucleotide(doubleKmer.revcompKmer); //must update this before advancing
        doubleKmer.forward(validExtension); 
        scanDist++;

        //handle backward junction case
        if(isJunction(doubleKmer.revcompKmer)){
            nextJunc = doubleKmer.revcompKmer;
            nextJuncDist = scanDist;
            nextJuncExtIndex = lastNuc;
            break;
        }
        scanDist++;

        //handle forward junction case
        if(isJunction(doubleKmer.kmer)){
            nextJunc = doubleKmer.kmer;
            nextJuncDist = scanDist;
            nextJuncExtIndex = 4;
            break; 
        }

        //if we're at the position where the sink would be, record the value for later use
        //REMOVED THE -2*jchecker->j
        if(scanDist == maxDist){ //subtract 2*j (which corresponds to j positions in the forward direction) to ensure jchecking works
            sinkKmer = doubleKmer.kmer;
        }
    }

    //if no junction was found, must be a sink!
    if(nextJunc == -1){
        return new kmer_type(sinkKmer);
    }

    //Otherwise, we found a junction, but it may or may not indicate an actual link.  Calculate overlap distance to verify.
    int otherMaxDist = getJunction(nextJunc)->dist[nextJuncExtIndex];
    int overlap = otherMaxDist + maxDist - nextJuncDist;
    if(overlap >= 0){
        return NULL; //if there is indicated overlap, this is not a sink.
    }
    return new kmer_type(sinkKmer); //if there is no indicated overlap between this and the next junction, we found a sink!
}

//Iterates through all of the junctions
//For each, bloom scans in every direction till another junction is hit or the stored distance is reached
//If there is no junction within that distance, a sink has been reached, and is added to the set of sinks.
//To be used FIRST after read scan
unordered_set<kmer_type>* JunctionMap::getSinks(){
    printf("Getting sinks.\n");
    unordered_set<kmer_type>* sinks = new unordered_set<kmer_type>({});
    kmer_type kmer;
    Junction junction;
    int juncsTested = 0;
    for(auto it = junctionMap.begin(); it != junctionMap.end(); it++, juncsTested++){
        kmer = it->first;
        junction = it->second;
        for(int i = 0; i < 5; i++){
            if( !junction.linked[i] && (i == 4 || junction.cov[i] > 0) ){ //need the || since we want to scan backwards, but we don't store backwards coverage
                kmer_type* sink = findSink(junction, kmer, i);
                if(sink){
                    kmer_type copy = *sink;
                    sinks->insert(copy);
                    delete(sink);
                }
            }
        }
        if(juncsTested % 10000 == 0) printf("Tested %d/%d junctions for sinks.\n", juncsTested, junctionMap.size());
    }
    return sinks;
}

//Iterates through all of the junctions, cleaning the non-complex ones.
//Every time a non-complex junction is found, the one real extension is added to the RealExtension set and the junction is destroyed
//To be used AFTER findSinks
unordered_map<kmer_type,int>* JunctionMap::getRealExtensions(){
    unordered_map<kmer_type,int>* realExtensions = new unordered_map<kmer_type,int>({});
    kmer_type kmer;
    Junction junction;
    int juncsTested = 0;
    int juncSize = junctionMap.size();
    for(auto it = junctionMap.begin(); it != junctionMap.end(); juncsTested++){
        kmer = it->first;
        junction = it->second;
        it++; //must do this before potentially erasing the element
        if (junction.numPathsOut() == 1){
            for(int i = 0; i < 4; i++){
                if(junction.cov[i] > 0){
                    //Must record both the base kmer and the valid extension to uniquely identify what's happening, since
                    //the same kmer can be a valid extension of one junction but an invalid extension of another junction.
                    (*realExtensions)[kmer] = i; 
                }
            }
            junctionMap.erase(kmer); 
        }
        if(juncsTested % 10000 == 0) printf("Tested %d/%d junctions for cFPs.\n", juncsTested, juncSize);
    }
    return realExtensions;
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

int JunctionMap::getNumSolidJunctions(int i){
  int count = 0;
  for(auto juncIt = junctionMap.begin(); juncIt != junctionMap.end(); juncIt++){
     if(juncIt->second.isSolid(i)){
        count++;
     }  
  }
  return count;
}

//Reads the distance from the junction corresponding to readKmer
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

//Assumes the two kmers are adjacent junctions on the same read.
//Links them!
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

void JunctionMap::createJunction(ReadKmer* readKmer){  
  createJunction(readKmer->getKmer());
}

void JunctionMap::createJunction(kmer_type kmer){  
  Junction newJunc;
  junctionMap.insert(std::pair<kmer_type, Junction>(kmer, newJunc));
}

void JunctionMap::killJunction(kmer_type kmer){
  junctionMap.erase(kmer);
}

//File format:
//One line for each junction.  On each line, the kmer is printed, then the junction is printed.  
//See Junction.h for junction print documentation.
void JunctionMap::writeToFile(string filename){
    ofstream jFile;
    jFile.open(filename);

    for(int i = 0; i < 5; i++){
        printf("There are %d junctions with solidity at least %d.\n", getNumSolidJunctions(i), i);
    }
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

bool JunctionMap::isJunction(ReadKmer* readKmer){
  return isJunction(readKmer->getKmer());
}

bool JunctionMap::isJunction(kmer_type kmer){
    return junctionMap.find(kmer) != junctionMap.end();
}

JunctionMap::JunctionMap(Bloom* bloo1, JChecker* jcheck, int read_length){
  junctionMap = {};
  bloom = bloo1;
  jchecker = jcheck;
  maxReadLength = read_length;
}