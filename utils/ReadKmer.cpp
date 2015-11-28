#include "ReadKmer.h"
#include <iostream>
#include <algorithm>
#include <string>

using std::string;


char* ReadKmer::directionAsString(){
    if(direction == FORWARD){
        return (char*)"forward";
    }
    else{
        return (char*)"backward";
    }
}

int ReadKmer::getMaxGuaranteedJ(bool dir){
    if(dir == FORWARD){
        return getDistToEnd()/2-1;
    }
    else{
        return getTotalPos()/2-1;
    }
}

//Returns number of forward operations needed to move to the last kmer on the read
int ReadKmer::getDistToEnd(){
 return read->length()*2- getTotalPos() - 2*sizeKmer+1;
}

//Returns the number of forward operations needed to go from the beginning of the read to this ReadKmer
int ReadKmer::getTotalPos(){
    return 2*pos + offset();
}

bool ReadKmer::onRead(){
    return (getTotalPos() >= 1) && (getDistToEnd()>=1);
}

kmer_type ReadKmer::getRevCompKmer(){
    if(direction == FORWARD){
        return doubleKmer.revcompKmer;
    }
    else{
        return doubleKmer.kmer;
    }
}

kmer_type ReadKmer::getKmer(){
    if(direction == FORWARD){
        return doubleKmer.kmer;
    }
    else{
        return doubleKmer.revcompKmer;
    }
}

//returns the offset of this kmer from the one at the same position facing backward.  
//Backward: 0
//Forward: 1
//used for calculating distances.
int ReadKmer::offset(){
    if(direction == FORWARD){
        return 1;
    }
    else{
        return 0;
    }
}

void ReadKmer::forward(){
    direction = !direction;
    if(direction == FORWARD){
        return; //switching from facing backward to forward doesn't entail a shift
    }
    int newNuc = 0;
    if(pos + sizeKmer <  read->length()){
        newNuc = NT2int((*read)[pos + sizeKmer]);   
    }
    doubleKmer.forward(newNuc);
    pos++;
}

void ReadKmer::advanceDist(int dist){
    for(int i = 0; i < dist; i++){
        forward();
    }
}

kmer_type ReadKmer::getCanon(){
    return doubleKmer.getCanon();
}

int ReadKmer::getExtensionIndex(bool dir){
    if(dir != direction){
        return 4; //backward index
    }
    return getRealExtensionNuc();
}

kmer_type ReadKmer::getExtension(int newNuc){
    return doubleKmer.getExtension(newNuc, direction);
}

//can be used as an index for the junction
int ReadKmer::getRealExtensionNuc(){
    if(direction == FORWARD){
        return NT2int((*read)[sizeKmer + pos]);  
    } 
    else{
        return revcomp_int(NT2int((*read)[pos-1]));
    }
}

kmer_type ReadKmer::getRealExtension(){
    return getExtension(getRealExtensionNuc());  
}

//Starts all the way at the front- facing off the read
ReadKmer::ReadKmer(string* theRead): doubleKmer(0){
    read = theRead;
    kmer_type kmer = 0;
    getFirstKmerFromRead(&kmer,&((*read)[0]));
    doubleKmer = DoubleKmer(kmer);
    pos = 0;
    direction = BACKWARD;
}

//Creates a double kmer corresponding to the given read, the index into the read, and the direction
ReadKmer::ReadKmer(string* theRead, int index, bool dir): doubleKmer(0){
    read = theRead;
    kmer_type kmer;
    getFirstKmerFromRead(&kmer,&((*read)[index]));
    doubleKmer = DoubleKmer(kmer);
    pos = index;
    direction = dir;
}

ReadKmer::ReadKmer(ReadKmer* toCopy): doubleKmer(toCopy->doubleKmer){
    read = toCopy->read;
    doubleKmer = toCopy->doubleKmer;
    pos = toCopy->pos;
    direction = toCopy->direction;
}

