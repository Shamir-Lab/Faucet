#ifndef JUNCTION_MAP
#define JUNCTION_MAP

#include <unordered_map>
#include <unordered_set>
#include <string>
#include "Kmer.h"
#include "Junction.h"
#include "Cap.h"
#include "ReadKmer.h"
#include "Bloom.h"
#include "JChecker.h"
#include "Kmer.h"
#include <fstream>
using std::ofstream;
using std::unordered_map;
using std::string;
using std::unordered_set;

class JunctionMap{

private: 
    Bloom* bloom;
    JChecker* jchecker; 
    int maxReadLength; //needed for finding sinks properly- tells you when to stop scanning   

public:
    unordered_map<kmer_type,Junction> junctionMap;  //stores the junctions themselves

    //Does bloom scans from the junctions that are not linked to other junctions.  
    //This identifies all the sinks in the graph and returns them as a set.
    unordered_set<kmer_type>* getSinks(); 

    //Assumes sinks were already found.
    //For each junction that has only one real extension, replace it with a marker that stores only the base kmer of that junction
    //and the one valid extension.  This serves the purpose of the cFP set for Minia.
    unordered_map<kmer_type, int>* getRealExtensions();

    //Scans forward from a junction at the specified index i.
    //If a sink is found, it is returned.  Otherwise, returns NULL
    kmer_type * findSink(Junction junc, kmer_type kmer ,int i);

    //Gets the valid extension of the given kmer based on the bloom filter and cFPs.  Uses JChecking! so this cuts off tips
    //Assume the given kmer is not a junction
    //Returns -1 if there is no valid extension
    //Returns -2 if there are multiple
    //ASSUMES NO CFP SET- since this is only done in findSinks, BEFORE the cFPs are found
    int getValidJExtension(DoubleKmer kmer);


    //File format:
    //One line for each junction.  On each line, the kmer is printed as a string, then the junction is printed.  
    //See Junction.h for junction print documentation.
    void writeToFile(string filename); 

    //Finds the junction associated with the given kmer and returns how far we can skip in the given direction from that junction
    int getSkipDist(ReadKmer* readKmer, bool direction);

    //Directly links two adjacent junctions from the same read
    void directLinkJunctions(ReadKmer* kmer1, ReadKmer* kmer2, Junction* junc1, Junction* junc2);

    int getNumComplexJunctions(); //Gets the number of junctions with more than one valid extension
    int getNumSolidJunctions(int i); //Gets the number of solid complex junctions, multiple valid extensions of coverage at least i
    int getNumJunctions();

    void createJunction(kmer_type kmer);
    void createJunction(ReadKmer* readKmer);
    bool isJunction(kmer_type kmer); //returns true if there is a junction at the given kmer
    bool isJunction(ReadKmer* readKmer); //same as above
    Junction* getJunction(ReadKmer kmer); //returns the junction located at the given kmer, or NULL if there is none
    Junction* getJunction(kmer_type kmer); //same as above
    void killJunction(kmer_type kmer); //removes the junction at the specified kmer, if there is one

    JunctionMap(Bloom* bloo, JChecker* jchecker, int maxReadLength);
};
#endif