#ifndef JUNCTION
#define JUNCTION

#include <iostream>
#include "Kmer.h"
using std::ofstream;

class Junction{
public:
    unsigned char dist[5];
    unsigned char cov[5];

    int numPathsOut();
    void writeToFile(ofstream* jFile);
    void update(int nucExt, unsigned char length); //if lengthFor is greater than the distance to the next junction bad things will happen
    void addCoverage(int nucExt);
    Junction();

};

#endif