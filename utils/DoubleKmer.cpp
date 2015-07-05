#include "DoubleKmer.h"

#include <algorithm>
#include <string>

using std::string;


char* DoubleKmer::directionAsString(){
    if(direction == FORWARD){
        return (char*)"forward";
    }
    else{
        return (char*)"backward";
    }
}

int DoubleKmer::getDistToEnd(){
 return read->length()*2- getTotalPos() - 2*sizeKmer+1;
}

int DoubleKmer::getTotalPos(){
    return 2*pos + offset();
}

bool DoubleKmer::onRead(){
    return (getTotalPos() >= 1) && (getDistToEnd()>=1);
}

kmer_type DoubleKmer::getKmer(){
    if(direction == FORWARD){
        return kmer;
    }
    else{
        return revcompKmer;
    }
}

//returns the offset of this kmer from the one at the same position facing backward.  
//Backward: 0
//Forward: 1
//used for calculating distances.
int DoubleKmer::offset(){
    if(direction == FORWARD){
        return 1;
    }
    else{
        return 0;
    }
}

void DoubleKmer::forward(){
    direction = !direction;
    if(direction == FORWARD){
        return; //switching from facing backward to forward doesn't entail a shift
    }
    int newNuc = NT2int((*read)[pos + sizeKmer]);
    kmer = next_kmer(kmer, newNuc, FORWARD);
    revcompKmer = next_kmer(revcompKmer, revcomp_int(newNuc), BACKWARD);
    pos++;
}

void DoubleKmer::backward(){
    int newNuc = NT2int((*read)[pos-1]);
    kmer = next_kmer(kmer, newNuc, BACKWARD);
    revcompKmer = next_kmer(revcompKmer, revcomp_int(newNuc), FORWARD);
    pos--;
}

void DoubleKmer::advanceDist(int dist){
    for(int i = 0; i < dist; i++){
        forward();
    }
}

kmer_type DoubleKmer::getCanon(){
    return std::min(kmer, revcompKmer);
}

int DoubleKmer::getExtensionIndex(bool dir){
    if(dir != direction){
        return 4; //backward index
    }
    return getRealExtensionNuc();
}

kmer_type DoubleKmer::getExtension(int newNuc){
    if(direction == FORWARD){
        return  next_kmer(kmer, newNuc, FORWARD);
    }
    else{
        return next_kmer(revcompKmer, newNuc, FORWARD); 
    }
}

//can be used as an index for the junction
int DoubleKmer::getRealExtensionNuc(){
    if(direction == FORWARD){
        return NT2int((*read)[sizeKmer + pos]);  
    } 
    else{
        return revcomp_int(NT2int((*read)[pos-1]));
    }
}

kmer_type DoubleKmer::getRealExtension(){
    return getExtension(getRealExtensionNuc());  
}

//Starts all the way at the front- facing off the read
DoubleKmer::DoubleKmer(string* theRead){
    read = theRead;

    getFirstKmerFromRead(&kmer,&((*read)[0]));
    revcompKmer = revcomp(kmer);
    pos = 0;
    direction = BACKWARD;
}

//Creates a double kmer corresponding to the given read, the index into the read, and the direction
DoubleKmer::DoubleKmer(string* theRead, int index, bool dir){
    read = theRead;
    getFirstKmerFromRead(&kmer,&((*read)[index]));
    revcompKmer = revcomp(kmer);
    pos = index;
    direction = dir;
}

DoubleKmer::DoubleKmer(DoubleKmer* toCopy){
    read = toCopy->read;
    kmer = toCopy->kmer;
    revcompKmer = toCopy->revcompKmer;
    pos = toCopy->pos;
    direction = toCopy->direction;
}
