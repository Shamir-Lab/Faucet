#include <fstream>
#include "Junction.h"
using std::ofstream;
using std::max;

void writeToFile(ofstream* jFile);

int Junction::numPathsOut(){
  int numPaths = 0;
  for(int i = 0; i < 4; i++){
    if(ext[i] > 0){
      numPaths++;
    }
  }
  return numPaths;
}

//Updates the junc info based on finding a path of length length from the extension nucExt
void Junction::update(int nucExt, int lengthFor, int lengthBack){
      ext[nucExt] = max(ext[nucExt],(unsigned char) lengthFor);
      cov[nucExt] = cov[nucExt]+1;
      back = max(back, (unsigned char) lengthBack);
}


//Kmer, then "ext,ext,ext,ext" then "cov,cov,cov,cov" for each of A,C,T,G in order.
void Junction::writeToFile(ofstream*jFile){
  for(int i = 0; i < 4; i++){
    *jFile << (int)ext[i] << "," ;
  }
  *jFile << " " ;
  for(int i = 0; i < 4; i++){
    *jFile << (int)cov[i] << "," ;
  }
  *jFile << " ";
  *jFile << (int)back << " ";
  *jFile << "\n";
}

Junction::Junction(int nucExt){
    ext[0] = 0, ext[1] = 0, ext[2] = 0, ext[3] = 0;
    cov[0] = 0, cov[1] = 0, cov[2] = 0, cov[3] = 0;
    back = 0;
    ext[nucExt] = 1;       
    cov[nucExt] = 1;    
}
