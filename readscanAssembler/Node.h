#ifndef NODE
#define NODE

#include <iostream>
#include "../utils/Kmer.h"
#include "../utils/Junction.h"
using std::ofstream;

class Node{
public:
    int dist[5];
    unsigned char cov[5];
    kmer_type nextJunc[5];

    int numPathsOut();
    void writeToFile(ofstream* jFile);
    void update(int nucExt, int length, kmer_type jID); //if lengthFor is greater than the distance to the next junction bad things will happen
    Node(Junction junc);
    Node();

};

#endif