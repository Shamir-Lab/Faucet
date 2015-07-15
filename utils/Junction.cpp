#include "Junction.h"
#include <fstream>
#include <limits.h>

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

bool Junction::isSolid(int threshold){
  int pathsOut = 0;
  for(int i = 0; i < 4; i++){
    if(cov[i] >= threshold){
      pathsOut++;
    }
  }
  return pathsOut > 1;
}

void Junction::addCoverage(int nucExt){
  cov[nucExt] = cov[nucExt] + 1;

  //handle overflow
  if(cov[nucExt] == 0){
    cov[nucExt] = UCHAR_MAX; 
  }
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
}

//explicitly set if it's a spacer or not
Junction::Junction(){
  for(int i  = 0; i < 5; i++){
    dist[i] = 0;
    cov[i] = 0;
    linked[i] = false;
  }
}

