#ifndef READ_KMER
#define READ_KMER

#include <string>

#include "Kmer.h"

using std::string;

//Used to represent a kmer with relation to a read.
//Stores the read, the kmer, revcomp kmer, and the position on the read (left end of the kmer). 
class DoubleKmer{
public:
    string* read;
    kmer_type kmer;
    kmer_type revcompKmer;
    int pos;
    bool direction;

    int getDistToEnd();
    int getTotalPos();
    char* directionAsString();
    bool onRead();
    kmer_type getKmer();
    void forward();
    void backward();
    void advanceDist(int dist);
    int getExtensionIndex(bool direction);
    kmer_type getExtension(int newNuc);// facing forward, for jchecking
    kmer_type getRealExtension();// facing forward, for jchecking
    int getRealExtensionNuc();
    kmer_type getCanon();
    int offset();

    //initializes the DoubleKmer to refer to the first kmer in the read
    DoubleKmer(string* theRead);
    //Creates a double kmer corresponding to the given read, the index into the read, and the direction
    DoubleKmer(string* theRead, int index, bool dir);
    DoubleKmer(DoubleKmer* toCopy);
};



#endif