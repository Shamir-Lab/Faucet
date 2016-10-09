#ifndef JUNCTION
#define JUNCTION

#include <iostream>
#include <string>
#include "Kmer.h"
using std::string;
using std::ofstream;

class Junction{
private:   
    unsigned char cov[4]; //the number of reads along this extension

public:
    //The following three fields are indexed by a value from 0-4.  If the value is from 0-3, it indicates the forward extension from adding
    //A, C, T, or G, respectively. If the index is 4, it refers to the backwards direction.
    unsigned char dist[5]; //the distance to the next adjacent junction, or the farthest scanned as of yet without hitting another junction
    bool linked[5]; //whether or not we found another junction along this extension

    //Returns an index that points to a valid path in the direction opposite the direction of the given index
    //If input is 4, it returns the valid path of 0-4
    //If input is not 0,1,2,3, returns 3
    int getOppositeIndex(int index);
    int numPathsOut(); //Returns the number of forward paths out of the junction with positive coverage
    bool isSolid(int threshold); //Returns true if at least 2 paths out of the junction have at least a threshold coverage.  
    
    
    //"dist dist dist dist dist  cov cov cov cov cov  linked linked linked linked linked " for each of A,C,T,G,Back, in order.
    string toString();
 
    //Updates the junction to point to given distance, if it's greater than the current distance stored.
    void update(int nucExt, unsigned char length);


    //Returns coverage along given extnsion
    //If nucExt == 4, returns the sum of the four coverage fields
    int getCoverage(int nucExt);

    void setCoverage(int nucExt, int coverage);


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