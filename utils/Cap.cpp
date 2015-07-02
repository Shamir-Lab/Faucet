#include <fstream>
#include "Cap.h"
using std::ofstream;
using std::max;

void writeToFile(ofstream* jFile);

//Kmer, then "ext,ext,ext,ext" then "cov,cov,cov,cov" for each of A,C,T,G in order.
void Cap::writeToFile(ofstream*jFile){
  *jFile <<"Distance: " << dist << " ";
  *jFile <<"Last ID: " << (long long) lastJunc << " ";
}

Cap Cap::extend(int extraDistance, kmer_type juncID){
    return *(new Cap(extraDistance+dist, juncID));
}

//explicitly set if it's a spacer or not
Cap::Cap(int distance, kmer_type juncID){
    dist=distance;
    lastJunc = juncID;
}

