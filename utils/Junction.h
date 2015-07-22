#ifndef JUNCTION
#define JUNCTION

#include <iostream>
#include <string>
#include "Kmer.h"
using std::string;
using std::ofstream;

class Junction{
public:
    //The following three fields are indexed by a value from 0-4.  If the value is from 0-3, it indicates the forward extension from adding
    //A, C, T, or G, respectively. If the index is 4, it refers to the backwards direction.
    unsigned char dist[5]; //the distance to the next adjacent junction, or the farthest scanned as of yet without hitting another junction
    unsigned char cov[5]; //the number of reads along this extension
    bool linked[5]; //whether or not we found another junction along this extension

    int numPathsOut(); //Returns the number of forward paths out of the junction with positive coverage
    bool isSolid(int threshold); //Returns true if at least 2 paths out of the junction have at least a threshold coverage.  
    
    
    //"dist dist dist dist dist  cov cov cov cov cov  linked linked linked linked linked " for each of A,C,T,G,Back, in order.
    string toString();
 
    //Updates the junction to point to given distance, if it's greater than the current distance stored.
    void update(int nucExt, unsigned char length);

    //Increments the coverage along the given extension by 1.
    void addCoverage(int nucExt);

    //Sets linked to true along the given extension.
    void link(int nucExt);

    //Initializes with 0 coverage, 0 distance, and linked false.
    Junction();

    //Get junction from string printout
    Junction(string juncString);
};

#endif