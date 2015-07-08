#ifndef JUNCTION
#define JUNCTION

#include <iostream>
#include "Kmer.h"
using std::ofstream;

class Junction{
public:
    //The following three fields are indexed by a value from 0-4.  If the value is from 0-3, it indicates the forward extension from adding
    //A, C, T, or G, respectively. If the index is 4, it referes to the backwards direction.
    unsigned char dist[5]; //the distance to the next adjacent junction, or the farthest scanned as of yet without hitting another junction
    unsigned char cov[5]; //the number of reads along this extension
    bool linked[5]; //whether or not we found another junction along this extension

    int numPathsOut(); //Returns the number of forward paths out of the junction with positive coverage

    //Format:
    //One line for the whole junction.
    //"Distances: 0,1,2,3,4, Coverages: 0,1,2,3,4, Linked: 0,1,2,3,4," 
    void writeToFile(ofstream* jFile); 
 
    //Updates the junction to point to given distance, if it's greater than the current distance stored.
    void update(int nucExt, unsigned char length);

    //Increments the coverage along the given extension by 1.
    void addCoverage(int nucExt);

    //Sets linked to true along the given extension.
    void link(int nucExt);

    //Initializes with 0 coverage, 0 distance, and linked false.
    Junction();
};

#endif