#ifndef JUNC_PAIR_SEARCH
#define JUNC_PAIR_SEARCH

#include "Kmer.h"
#include "ContigJuncList.h"
#include <sstream>
#include <deque>
#include <algorithm> 
#include <iostream>  
using std::stringstream;

struct JuncPair{
    JuncPair(kmer_type km1, kmer_type km2): kmer1(km1), kmer2(km2){}
    kmer_type kmer1;
    kmer_type kmer2;
    friend bool operator==(JuncPair a, JuncPair b) { 
        return a.kmer1 == b.kmer1 && a.kmer2 == b.kmer2; 
    };

};

namespace std {
  template <> struct hash<JuncPair>
  {
    size_t operator()(const JuncPair & x) const
    {
        return (std::hash<uint64_t>()(x.kmer1) ^ (std::hash<uint64_t>()(x.kmer2) << 1) >> 1);
    }
  };
}


//Stores all the info involved in finding a junction candidate by searching from a node
class JuncResult{
public:
    JuncResult(kmer_type km, int dist, int cov): kmer(km), distance(dist), coverage(cov){} 
    kmer_type kmer;
    int distance;
    int coverage;
    friend bool operator<(JuncResult a, JuncResult b);
    friend bool operator>(JuncResult a, JuncResult b);
};

//Stores all the info involved in finding a junction candidate by searching from a node
struct JuncPairResult{
    JuncPairResult(JuncPair p, int dist, int cov): pair(p), distance(dist), coverage(cov){} 
    JuncPair pair;
    int distance;
    int coverage;
};

#endif