#include <fstream>
#include "Node.h"
using std::ofstream;
using std::max;

void writeToFile(ofstream* jFile);

void Node::deletePath(int index){
  dist[index] = 0;
  cov[index] = 0;
  nextJunc[index] = -1;
  backIndex[index] = -1; //the index along which the next node points back to this one at this extension
}

int Node::numPathsOut(){
  int numPaths = 0;
  for(int i = 0; i < 4; i++){
    if(cov[i] > 0){
      numPaths++;
    }
  }
  return numPaths;
}

//Updates the junc info based on finding a path of length length from the extension nucExt
//With given backIndex
void Node::update(int nucExt, int lengthFor,  kmer_type jID, int bI){
      dist[nucExt] = max(dist[nucExt], lengthFor);
      nextJunc[nucExt] = jID;
      backIndex[nucExt] = (char)bI;
}


//Kmer, then "ext,ext,ext,ext" then "cov,cov,cov,cov" for each of A,C,T,G in order.
void Node::writeToFile(ofstream*jFile){
  *jFile <<"Distances: ";
  for(int i = 0; i < 5; i++){
    *jFile << (int)dist[i] << "," ;
  }
  *jFile << " Coverages: " ;
  for(int i = 0; i < 5; i++){
    *jFile << (int)cov[i] << "," ;
  }
  *jFile << " Back indices: " ;
  for(int i = 0; i < 5; i++){
    *jFile << (int)backIndex[i] << "," ;
  }
  *jFile << " IDs: ";
  for(int i = 0; i < 5; i++){
    if(nextJunc[i] == -1){
      *jFile << "x,";
    }
    else{
      *jFile << print_kmer(nextJunc[i]) << "," ;
    }
  }
}

Node::Node(Junction junc){
  for(int i  = 0; i < 5; i++){
    dist[i] = junc.dist[i];
    cov[i] = junc.getCoverage(i);
    nextJunc[i] = -1;
    backIndex[i] = -1;
  }
}

Node::Node(){
  for(int i  = 0; i < 5; i++){
    dist[i] = 0;
    cov[i] = 0;
    nextJunc[i] = -1;
    backIndex[i] = -1;
  }
}

