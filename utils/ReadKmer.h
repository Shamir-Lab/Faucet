#ifndef READ_KMER
#define READ_KMER

#include <string>

#include "Kmer.h"
#include "DoubleKmer.h"

using std::string;

//Used to represent a kmer with relation to a read.
//Stores the read, the kmer and revcomp as a DoubleKmer and the position on the read (left end of the kmer). 
class ReadKmer{
public:
    //basic fields
    string* read;
    DoubleKmer doubleKmer;
    int pos;
    bool direction;

    //The real extension of the ReadKmer definitely j-checks to at least the return value, based on its position on the read.
    int getMaxGuaranteedJ(bool dir);

    int getDistToEnd(); //returns dist to end
    int getTotalPos(); //returns dist to start 

    //returns true if the kmer is at a valid position on the read.  Does not include the first backward kmer or the last forward kmer.
    bool onRead();    

    kmer_type getKmer();
    kmer_type getRevCompKmer();

    //moves one position forward. If facing BACKWARD, simply changes the kmer to face  FORWARD.
    //This may entail going from facing BACKWARD to FORWARD
    void forward(); 
    void advanceDist(int dist); //calls forward repeatedly

    //Gets the index on a corresponding junction which corresponds to going along the read in the given direction
    //e.g. if the ReadKmer faces forward and direction is BACKWARD, returns 4.
    //If the ReadKmer faces forward and the direction is FORWARD and the next real nucleotide on the read is 'T', returns 2
    int getExtensionIndex(bool direction); 
    kmer_type getExtension(int newNuc);// Gets the next extension in the direction its facing, for the given nucleotide extension
    kmer_type getRealExtension();// gets the next real kmer in the direction the ReadKmer is pointing
    int getRealExtensionNuc(); //gets the next real nucleotide in the direction the ReadKmer is pointing
    kmer_type getCanon(); 
    int offset(); //gets the contribution of the direction to the distance- 0 if BACKWARD, 1 if FORWARD

    char* directionAsString(); //for printing

    ReadKmer(string* theRead); //initializes the DoubleKmer to refer to the first kmer in the read
    ReadKmer(string* theRead, int index, bool dir);//Creates a double kmer corresponding to the given read, the index into the read, and the direction
    ReadKmer(ReadKmer* toCopy); //copy construct

};



#endif