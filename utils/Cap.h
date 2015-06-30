#ifndef CAP
#define CAP

#include <iostream>
#include "Kmer.h"
using std::ofstream;

class Cap{
public:
    int dist;
    kmer_type lastJunc;

    void writeToFile(ofstream* jFile);

    Cap extend(int dist, kmer_type juncID);

    Cap(int distance,kmer_type juncID);
};

#endif