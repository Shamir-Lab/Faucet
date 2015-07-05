#include "DoubleKmer.h"

#include <algorithm>

void DoubleKmer::forward(int nuc){
    kmer = next_kmer(kmer, nuc, FORWARD);
    revcompKmer = next_kmer(revcompKmer, revcomp_int(nuc), BACKWARD);
}

kmer_type DoubleKmer::getExtension(int nuc, bool direction){
    if(direction == FORWARD){
        return  next_kmer(kmer, nuc, FORWARD);
    }
    else{
        return next_kmer(revcompKmer, nuc, FORWARD); 
    }
}

kmer_type DoubleKmer::getCanon(){
    return std::min(kmer, revcompKmer);
}

DoubleKmer::DoubleKmer(kmer_type forwardKmer){
    kmer = forwardKmer;
    revcompKmer = revcomp(kmer);
}