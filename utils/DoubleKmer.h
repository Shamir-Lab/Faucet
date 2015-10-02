#ifndef DOUBLE_KMER
#define DOUBLE_KMER

#include "Kmer.h"

class DoubleKmer{

public:
    kmer_type kmer;
    kmer_type revcompKmer;
    
    void forward(int nuc);
    
    //Takes as input the nucleotide extension and the direction
    //If the direction is BACKWARD, the nucleotide should be given as seen in the reverse direction 
    //It will not be complemented within the function.
    kmer_type getExtension(int nuc, bool dir);

    kmer_type getCanon();

    void reverse();
    
    DoubleKmer(kmer_type forwardKmer);
};
#endif 