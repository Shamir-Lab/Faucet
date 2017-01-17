#include "Junction.h"
#include <fstream>
#include <limits.h>
#include <sstream>
#include <iostream>
using std::ofstream;
using std::max;
using std::istringstream;
using std::stringstream;

int Junction::getOppositeIndex(int index){
  if(index < 4){
    return 4;
  }
  else{
    for(int i = 0; i < 4; i++){
      if(cov[i] > 0){
        return i;
      }
    }
  }
}

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

int Junction::getCoverage(int nucExt){
  if(nucExt < 4){
    return (int)cov[nucExt];
  }
  return (int)cov[0] + (int)cov[1] + (int)cov[2] + (int)cov[3];
}

void Junction::setCoverage(int nucExt, int coverage){
  cov[nucExt] = coverage;
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
  if (nucExt < 0 || nucExt > 4) {
    std::cout << "ERROR: Invalid nucExt in Junction::update: " << nucExt << "\n";
  }
  dist[nucExt] = max(dist[nucExt], lengthFor);

}

//"dist dist dist dist dist  cov cov cov cov cov  linked linked linked linked linked " for each of A,C,T,G,Back, in order.
string Junction::toString(){
  stringstream stream;
  for(int i = 0; i < 5; i++){
    stream << (int)dist[i] << " " ;
  }
  stream << " ";
  for(int i = 0; i < 5; i++){
    stream << getCoverage(i) << " " ;
  }
  stream << " ";
  for(int i = 0; i < 5; i++){
    stream << linked[i] << " " ;
  }
  return stream.str();
}

//explicitly set if it's a spacer or not
Junction::Junction(){
  for(int i  = 0; i < 4; i++){
    dist[i] = 0;
    cov[i] = 0;
    linked[i] = false;
  }
  dist[4] = 0;
  linked[4] = false;
}

//Get junction from string printout
Junction::Junction(string juncString){
  istringstream iss(juncString);
  string val;
  for(int i = 0; i < 5; i++){
    iss >> val;
    dist[i] = stoi(val);
  }
  for(int i = 0; i < 4; i++){
    iss >> val;
    cov[i] = stoi(val);
  }
  iss >> val;
  for(int i = 0; i < 5; i++){
    iss >> val;
    linked[i] = stoi(val);
  }
}

