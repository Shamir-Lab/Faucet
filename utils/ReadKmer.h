#ifndef READ_KMER
#define READ_KMER

#include <string>

#include "Kmer.h"
#include "DoubleKmer.h"

using std::string;

//Used to represent a kmer with relation to a read.
//Stores the read, the kmer, revcomp kmer, and the position on the read (left end of the kmer). 
class ReadKmer{
public:
    string* read;
    DoubleKmer doubleKmer;
    int pos;
    bool direction;

    int getMaxGuaranteedJ(bool dir);//if this is the junction, the real extension definitely j-checks to at least the return value.
    int getDistToEnd();
    int getTotalPos();
    char* directionAsString();
    bool onRead();
    kmer_type getKmer();
    kmer_type getRevCompKmer();
    void forward();
    void advanceDist(int dist);
    int getExtensionIndex(bool direction);
    kmer_type getExtension(int newNuc);// facing forward, for jchecking
    kmer_type getRealExtension();// facing forward, for jchecking
    int getRealExtensionNuc();
    kmer_type getCanon();
    int offset();

    //initializes the DoubleKmer to refer to the first kmer in the read
    ReadKmer(string* theRead);
    //Creates a double kmer corresponding to the given read, the index into the read, and the direction
    ReadKmer(string* theRead, int index, bool dir);
    ReadKmer(ReadKmer* toCopy);
};



#endif