#include <fstream>
#include "Junction.h"
using std::ofstream;
using std::max;

void writeToFile(ofstream* jFile);

int Junction::numPathsOut(){
  int numPaths = 0;
  for(int i = 0; i < 4; i++){
    if(cov[i] > 0){
      numPaths++;
    }
  }
  return numPaths;
}

void Junction::link(int nucExt){
  linked[nucExt] = true;
}

void Junction::addCoverage(int nucExt){
  cov[nucExt] = cov[nucExt] + 1;
}

//Updates the junc info based on finding a path of length length from the extension nucExt
void Junction::update(int nucExt, unsigned char lengthFor){
      dist[nucExt] = max(dist[nucExt], lengthFor);
}


//Kmer, then "ext,ext,ext,ext" then "cov,cov,cov,cov" for each of A,C,T,G in order.
void Junction::writeToFile(ofstream*jFile){
  *jFile <<"Distances: ";
  for(int i = 0; i < 5; i++){
    *jFile << (int)dist[i] << "," ;
  }
  *jFile << " Coverages: " ;
  for(int i = 0; i < 5; i++){
    *jFile << (int)cov[i] << "," ;
  }
  *jFile << " Linked: ";
  for(int i = 0; i < 5; i++){
    *jFile << linked[i] << "," ;
  }
  *jFile << "\n";
}

//explicitly set if it's a spacer or not
Junction::Junction(){
  for(int i  = 0; i < 5; i++){
    dist[i] = 0;
    cov[i] = 0;
    linked[i] = false;
  }
}

