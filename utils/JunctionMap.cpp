#include "JunctionMap.h"
#include <string>
#include <fstream>
#include <sstream>
#include <assert.h>
using std::ofstream;
using std::ifstream;
using std::istringstream;
using std::string;

ContigGraph* JunctionMap::buildContigGraph(){
    ContigGraph* contigGraph = new ContigGraph();
    printf("Building contig graph.\n");
    
    printf("Building branching regions.\n");    
    buildBranchingPaths(contigGraph);

    printf("Destroying complex junctions.\n");
    destroyComplexJunctions();

    printf("Building linear regions.\n");
    buildLinearRegions(contigGraph);

    printf("Checking graph.\n");
    contigGraph->checkGraph();
    
    // printf("Switching to node vector.\n");
    // contigGraph->switchToNodeVector();

    printf("Done building contig graph.\n");
    return contigGraph;
}

void JunctionMap::buildBranchingPaths(ContigGraph* contigGraph){
    for(auto it = junctionMap.begin(); it != junctionMap.end(); it++){ //for each junction
        kmer_type kmer = it->first;
        Junction junction = it->second;
        if(junction.numPathsOut() > 1){ //if the junction is complex
            //printf("Kmer %s\n", print_kmer(kmer));
            ContigNode* startNode = contigGraph->createContigNode(kmer, junction);
            for(int i = 0; i < 5; i++){ 
                if(junction.getCoverage(i) > 0 && !startNode->contigs[i]){ //for every valid path out which doesn't already have a contig
                    //printf("Building contig from index %d\n", i);
                    Contig* contig = getContig(junction, kmer, i);
                    ContigNode* otherNode = nullptr;
                    kmer_type far_kmer = contig->getSideKmer(2);
                    if (isJunction(far_kmer)){

                        if(far_kmer!=kmer){ //complex junction- should create a contig node
                            //printf("Path builder found a node: %s\n", print_kmer(far_kmer));
                            Junction* far_junc = getJunction(far_kmer);
                            otherNode = contigGraph->createContigNode(far_kmer, *far_junc);//create a contig on the other side if it doesn't exist yet
                            contig->setEnds(startNode, contig->ind1, otherNode, contig->ind2);
                        }else{ // inverted repeat - looped back to source junction
                            contig->setEnds(startNode, contig->ind1, startNode, contig->ind2);                            
                        }
                    }
                }
            }
        }
    }  
}

void JunctionMap::destroyComplexJunctions(){
    for(auto it = junctionMap.begin(); it != junctionMap.end(); ){ //for each junction
        kmer_type kmer = it->first;
        Junction junction = it->second;
        if ( junction.numPathsOut() > 1 ) { //if the junction is complex
            it++;
            killJunction(kmer);
        }
        else{
            it++;
        }
    }  
}

void JunctionMap::buildLinearRegions(ContigGraph* contigGraph){
    for(auto it = junctionMap.begin(); it != junctionMap.end(); it++){ //for each junction
        kmer_type kmer = it->first;
        Junction junction = it->second;
        if ( junction.numPathsOut() == 1 ) { //if the junction is not complex
            Contig* forwardContig;
            Contig* backwardContig;
            for ( int i = 0 ; i < 4 ; i++ ){ 
                if ( junction.getCoverage(i) > 0 ) { //for every valid path out- should only be 1
                    forwardContig = getContig(junction, kmer, i);
                }
            }
            backwardContig = getContig(junction, kmer, 4);

            contigGraph->addIsolatedContig(*backwardContig->concatenate(forwardContig, 1, 1));

            delete(forwardContig);
            delete(backwardContig);
        }
        else{
            printf("ERROR: shouldn't be any complex junctions left during linear region build. %d paths out\n", junction.numPathsOut());
        }
    }  
}

//Gets the contig from this junction to the next complex junction or sink
Contig* JunctionMap::getContig(Junction startJunc, kmer_type startKmer, int startIndex){
    //just for debugging
    // std::list<BfSearchResult> results = {};
    // std::list<Junction> junctions = {};
    // junctions.push_back(startJunc);

    //needed for loop
    std::list<kmer_type> kmers_to_destroy = {};
    Junction junc = startJunc;
    kmer_type kmer = startKmer;
    std::deque<unsigned char> coverages;
    std::deque<unsigned char> distances;
    int index = startIndex;

    coverages.push_back(junc.getCoverage(startIndex));
    string contigString(print_kmer(kmer));

    if(index == 4) contigString = print_kmer(revcomp(kmer));
    BfSearchResult result;

    bool done = false;
    Contig* contig = new Contig();

    while(!done){
        result = findNeighbor(junc, kmer, index);

        if(result.contig.length() > sizeKmer){
            contigString += result.contig.substr(sizeKmer); //trim off the first k chars to avoid repeats 
        }
        distances.push_back((unsigned char) result.distance);

        if(result.isNode){

            Junction nextJunc = *getJunction(result.kmer);
           
            coverages.push_back(nextJunc.getCoverage(result.index));

            if (nextJunc.numPathsOut() == 1){
                junc = nextJunc;
                kmer = result.kmer;
                index = junc.getOppositeIndex(result.index);
                // kmers_to_destroy.push_back(kmer); 

                if (result.kmer == startKmer){
                    done = true;
                }
                else{ // isolated self-loop --> nothing to destroy
                    kmers_to_destroy.push_back(kmer);
                }
            }
            else{
                done = true;
            }
        }
        else{
            done = true;
        }
    }

    contig->setContigJuncs(ContigJuncList(contigString, distances, coverages));
    if(result.isNode){
        contig->setIndices(startIndex, result.index);
    }
    else{
        contig->setIndices(startIndex, 4);
    }

    //destroy all kmers found
    for(auto it = kmers_to_destroy.begin(); it != kmers_to_destroy.end(); it++){
        kmer_type toDestroy = *it;
        if (toDestroy != startKmer){
            killJunction(toDestroy);
        }
    }

    kmers_to_destroy.clear();

     return contig;
}


//TESTED that it always returns the correct type of answer, and that it's correct if it's a sink.\
//Scans forward from junction junc at index i with bloom filter
//If it hits another junction at or before the distance specified by the given junction, returns a "node" result with that junction
//If it does not, it keeps scanning until it hits another junction or an actual sink
//If it hits a sink, it returns it.  If it hits a junction, it tests how far that junction points along the path.
//Based on the indicated overlap, it either decides the entire intermediate sequence is real or the connection is a 
//false positive connection.  Then returns either a sink or a node result.
BfSearchResult JunctionMap::findNeighbor(Junction junc, kmer_type startKmer, int index){
    DoubleKmer doubleKmer(startKmer);
    kmer_type kmer;
    // represents distance from the start junction- 1 is the reverse kmer after the junction,
    // 2 is the first forward extension, 3 is the second reverse kmer, etc.
    int dist = 1; 
    // Junctions record the farthest distance any read has extended past each extension- this limits the search space of the 
    // findNeighbor call
    int maxDist = junc.dist[index];

    string contig("");

    int lastNuc; //stores the last nuc so we know which extension we came from. 

    //First, process the first 1-2 kmers in order to reach the first kmer from which can properly bloom scan.  
    //This is different for forwards and backward extensions.
    if(index == 4){ // the backwards extension of the junction
        doubleKmer.reverse(); // turn the input kmer around to start the backward search
        for(int i = 0; i < sizeKmer; i++){ // add the nucleotides of the start kmer to the found contig
            contig.append(1,getNucChar(code2nucleotide(doubleKmer.kmer, i)));
        }
        //If we're searching backwards, we only need to specially process the reverse kmer, and then scan from there
        if(isJunction(doubleKmer.kmer)){ // if the initial reverse kmer is already a junction, we found a contig! 
            return BfSearchResult(doubleKmer.kmer, 
                true,  // is a node
                4, // the junctions point away from each other; thus, the start junction is on the backward extension of the found junction
                1, // Distance 1
                contig); 
        }
    }
    else{ // forward extension case
        // here the scan startpoint is the next forward kmer- but since we're at a junction we must get there manually using the 
        // given index, no bloom scan possible till after that first extension.
        lastNuc = first_nucleotide(doubleKmer.revcompKmer);  // save the original last nucleotide, since it will be lost by advancing forward
        contig += getNucChar(code2nucleotide(doubleKmer.kmer, 0)); // add nucleotide to contig
        doubleKmer.forward(index); // advance the kmer! Now it's at the extension of the junction instead of the junction itself
        for(int i = 0; i < sizeKmer; i++){
            contig += getNucChar(code2nucleotide(doubleKmer.kmer, i)); // Add the current k characters to the contig
        }
        if(isJunction(doubleKmer.revcompKmer)){
            return BfSearchResult(doubleKmer.revcompKmer, 
                true, // is a node
                lastNuc, // this nucleotide is the extension which we would follow from the found junction to get back to the start
                1, // distance 1
                contig);
        }
        if (maxDist == 1) { // if maxDist is 1 but no junction found, make a sink instead

        }
        dist = 2; // if we arrived here, we didn't return above, and can search at distance 2.
        if(isJunction(doubleKmer.kmer)){
            return BfSearchResult(doubleKmer.kmer, 
                true,  // is a node
                4,  // backward extension
                2, // distance 2
                contig);
        }
    }
    
    //After the scan loop, these will store the needed info about the junction found
    BfSearchResult sinkResult, juncResult;
    int nextJuncDist, nextJuncExtIndex;

    //if we're at or past the position where the sink would be, record the value for later use
    if(dist >= maxDist){ //REMOVED THE - 2 * jchecker->j

        if(dist > maxDist){ // dist >= maxDist should really only occur iff dist == maxDist. Otherwise maxDist pointed to a reverse kmer as a sink, which is wrong.
            printf("Error: dist %d is greater than maxDist %d.\n", dist, maxDist);
            std::cout << "Searching from kmer " << print_kmer(startKmer) << "\n";
            printDistAndExtension(dist, maxDist, index);
        }
        sinkResult = BfSearchResult(doubleKmer.kmer, false, 5, dist, contig); // set up sink result if we are at dist == maxDist -> the sink
    }

    // Scan forward until maxDist
    // Here we know we should not ever find a sink, and we can also not find a junction, until maxDist, without violating our assumptions!
    while(dist < maxDist){ 

        //move forward if possible
        int validExtension = getValidJExtension(doubleKmer);
        assert(validExtension != -1 && validExtension != -2); // We cannot find a sink closer than maxDist!!
    
        lastNuc = first_nucleotide(doubleKmer.revcompKmer); //must update this before advancing
        doubleKmer.forward(validExtension); 
        contig += getNucChar(validExtension); //include this in the contig regardless of which way the end junction faces
        
        dist++;
        //handle backward junction case
        if (dist == maxDist) { 
            break; //if at maxDist, break the loop.
        }
        if (isJunction(doubleKmer.revcompKmer)) {
            std::cout << "Found a sequence overlap with no indication of a linking read! Assuming no real connection.\n";
            printDistAndExtension(dist, maxDist, index);
        }


        dist++;
        //handle forward junction case
        if (dist == maxDist) {
            break; // break if at maxDist
        }
        if (isJunction(doubleKmer.kmer)) {
            std::cout << "Found a sequence overlap with no indication of a linking read! Assuming no real connection.\n";
            printDistAndExtension(dist, maxDist, index);
        }
    }

    // If we get here, we must not have found a junction. Make sure it looks like a sink, and return it.
    assert(dist == maxDist); // This just checks the loop logic
    
    if(isJunction(doubleKmer.kmer)){ // if found a junction, return it!
        return BfSearchResult(doubleKmer.kmer, true, 4, dist, contig);
    } 

    // If there was no junction, this must be a sink. Check that the parity works out for the sink to point away from the junction
    if(index == 4) { // These tests make sure sink points away from junction, since at this point we know we found a sink.
        if (dist % 2 != 1) {
            printDistAndExtension(dist, maxDist, index);
        }
        assert(dist % 2 == 1);  
    } else{
        if (dist % 2 != 0) {
            printDistAndExtension(dist, maxDist, index);
        }
        assert(dist % 2 == 0);
    }

    return BfSearchResult(doubleKmer.kmer, false, 5, dist, contig); // return the sink!
}

void JunctionMap::printDistAndExtension(int dist, int maxDist, int index) {
    std::cout << "Dist: " << dist << ", maxDist: " << maxDist << "\n";
    std::cout << "Ext: " << index << "\n";
}

//Gets the valid extension of the given kmer based on the bloom filter and cFPs.  Uses JChecking! so this cuts off tips
//Assume the given kmer is not a junction
//Returns -1 if there is no valid extension
//Returns -2 if there are multiple
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

//Returns true if multiple extensions of the given kmer jcheck
//Assumes the given kmer is in the BF
bool JunctionMap::isBloomJunction(kmer_type kmer){
    kmer_type ext;
    int pathCount = 0;
    for(int i = 0; i < 4; i++){
        ext = next_kmer(kmer, i, FORWARD);
        if(jchecker->jcheck(ext)){
            pathCount++;
        }
    }
    return pathCount > 1;
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
    return nullptr;
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
//Assumes kmer1 corresponds to junc1, and kmer2 to junc2
//Links them!
void JunctionMap::directLinkJunctions(ReadKmer* kmer1, ReadKmer* kmer2, Junction* junc1, Junction* junc2){
    int ext1 = kmer1->getExtensionIndex(FORWARD);
    int ext2 = kmer2->getExtensionIndex(BACKWARD);
    
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
        jFile << juncIt->second.toString();
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

//builds junction map from junction map file.
//Assumes JunctionMap was just initialized
void JunctionMap::buildFromFile(string junction_file){
    junctionMap = {};
    ifstream jFile(junction_file);

    printf("Reading from Junction file to build junction map.\n");
    kmer_type kmer;
    Junction junc;
    string word, line;
    while(getline(jFile, line)){
        istringstream iss(line);

        iss >> word;
        getFirstKmerFromRead(&kmer, &word[0]);

        getline(iss, word);
        junc = Junction(word);

        junctionMap[kmer] = junc;
    }
    jFile.close();
}