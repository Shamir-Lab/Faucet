//
//  Bloom.h
//
//  Created by Guillaume Rizk on 9/02/12.
//

#ifndef Bloom_h
#define Bloom_h
#include <stdlib.h>
#include <inttypes.h>
#include <stdint.h>
#include <set>
#include <string>
#include "Kmer.h"
#include "ReadKmer.h"
// not using kmer_type from Kmer.h because I don't want this class to depend on Kmer.h
#ifdef _largeint
#include "LargeInt.h"
typedef LargeInt<KMER_PRECISION> bloom_elem;
#else
#ifdef _ttmath
#include "ttmath/ttmath.h"
typedef ttmath::UInt<KMER_PRECISION> bloom_elem;
#else
#if (! defined kmer_type) || (! defined _LP64)
typedef uint64_t bloom_elem;
#else
typedef kmer_type bloom_elem;
#endif
#endif
#endif

#define NSEEDSBLOOM 10
#define CUSTOMSIZE 1

static const int bits_per_char = 0x08;    // 8 bits in 1 char(unsigned)
static const unsigned char bit_mask[bits_per_char] = {
    0x01,  //00000001
    0x02,  //00000010
    0x04,  //00000100
    0x08,  //00001000
    0x10,  //00010000
    0x20,  //00100000
    0x40,  //01000000
    0x80   //10000000
};


static const uint64_t rbase[NSEEDSBLOOM] =
{
    0xAAAAAAAA55555555ULL, 
    0x33333333CCCCCCCCULL,
    0x6666666699999999ULL,
    0xB5B5B5B54B4B4B4BULL,
    0xAA55AA5555335533ULL,
    0x33CC33CCCC66CC66ULL,
    0x6699669999B599B5ULL,
    0xB54BB54B4BAA4BAAULL,
    0xAA33AA3355CC55CCULL,
    0x33663366CC99CC99ULL
};


class Bloom{
    
protected:
    
#ifdef _largeint
    inline uint64_t hash_func(LargeInt<KMER_PRECISION> elem, int num_hash);
#endif
#ifdef _ttmath
    inline uint64_t hash_func(ttmath::UInt<KMER_PRECISION> elem, int num_hash);
#endif
#ifdef _LP64
    inline uint64_t hash_func(__uint128_t key, int num_hash);
#endif
    inline uint64_t hash_func(uint64_t key, int num_hash);
    inline void generate_hash_seed(); 
    uint64_t user_seed;
    uint64_t seed_tab[NSEEDSBLOOM];
    uint64_t char_hash[2][4];
  
    uint64_t getCharHash(int key, int num_hash);
    uint64_t getLastCharHash(uint64_t key, int num_hash);

    int n_hash_func;
    uint64_t nchar;
    int k;

    //only relevant for a fake bloom
    bool fake;

    int hashSize;
    uint64_t bloomMask;
    std::set<bloom_elem> valid_set;
    std::set<uint64_t> valid_hash0;
    std::set<uint64_t> valid_hash1;

public:
    int getHashSize();
    uint64_t getBloomMask();

    unsigned char * blooma;


    /**********************************************************************************
    These are the important things that are currently being used.
    ***********************************************************************************/
    
    float weight(); //returns the proportion of 1's in the filter.  So should be between 0.0 and 1.0
    

    Bloom* create_bloom_filter_2_hash(uint64_t estimated_items, float fpRate); //creates for two hash functions and given fpRate

    Bloom* create_bloom_filter_optimal(uint64_t estimated_items, float fpRate); //creates for smallest size given the fpRate

    //loads all the kmers in the reads file into the bloom filter.
    //Input is assumed to be a raw string for each read, one per line.
    void load_from_reads(const char* reads_filename); 

    //loads all the kmers in the kmers file into the bloom filter.
    //Input is assumed to be one kmer per line as a string.
    void load_from_kmers(const char* kmers_filename); 
    
    //The hash function minia used and we are now using.
    inline uint64_t oldHash(uint64_t key, int num_hash){
      uint64_t hash = seed_tab[num_hash];
      hash ^= (hash <<  7) ^  key * (hash >> 3) ^ (~((hash << 11) + (key ^ (hash >> 5))));
      hash = (~hash) + (hash << 21); // hash = (hash << 21) - hash - 1;
      hash = hash ^ (hash >> 24);
      hash = (hash + (hash << 3)) + (hash << 8); // hash * 265
      hash = hash ^ (hash >> 14);
      hash = (hash + (hash << 2)) + (hash << 4); // hash * 21
      hash = hash ^ (hash >> 28);
      hash = hash + (hash << 31);
      return hash &= bloomMask;
    }

    //Add an element using the old hash function
    inline int oldAdd(bloom_elem elem)
    {
        uint64_t hA,hB;

        hA = oldHash(elem, 0);
        hB = oldHash(elem, 1);

        add(hA, hB);
    }

    //Check whether an element is contained using the old hash function
    inline int oldContains(bloom_elem elem)
    {
        if(fake){
            return (valid_set.find(elem) != valid_set.end());
        }
        uint64_t hA,hB;

        hA = oldHash(elem, 0);
        hB = oldHash(elem, 1);

        return contains(hA, hB);
    }


    /**********************************************************************************
    Most of the below is not currently used.  Much of it is for incremental hashing.
    ***********************************************************************************/

    //rotates hash to the right by dist. Assume 0 < dist < hashSize
    inline uint64_t  rotate_right(uint64_t  hash, int dist){
        dist %= hashSize;
        return ((hash >> dist) | (hash << (hashSize - dist))) & bloomMask;
    }

    //rotates hash to the right by dist. Assume 0 < dist < hashSize
    inline uint64_t  rotate_left(uint64_t  hash, int dist){
        dist %= hashSize;
        return ((hash << dist) | (hash >> (hashSize - dist))) & bloomMask;
    }

    //only for num_hash = 0 or 1
    uint64_t get_rolling_hash(uint64_t key, int num_hash);

    inline uint64_t roll_hash(uint64_t oldHash, int oldC, int newC, int num_hash){
      return rotate_left(oldHash ^ getCharHash(oldC, num_hash), 1) ^ rotate_right(getLastCharHash(newC, num_hash), k-1);
    }
    
    inline void advance_hash(char* read, uint64_t * hash0, uint64_t * hash1, int startPos, int endPos){
        for(int i = startPos; i < endPos; i++){
            *hash0 = roll_hash(*hash0, NT2int(read[i]), NT2int(read[i+sizeKmer]), 0);
            *hash1 = roll_hash(*hash1, NT2int(read[i]), NT2int(read[i+sizeKmer]), 1);
        }
    }

    inline void add(bloom_elem elem)
    {
        uint64_t hA,hB;

        hA = get_rolling_hash(elem, 0);
        hB = get_rolling_hash(elem, 1);

        add(hA, hB);    
    }


    inline void add(uint64_t h0, uint64_t h1)
    {
        uint64_t h = h0;
        for(int i=0; i<n_hash_func; i++, h += h1)
        {
            h %= tai;
            blooma [h >> 3] |= bit_mask[h & 7];
        }
    }


    inline int contains(bloom_elem elem)
    {
        if(fake){
            return (valid_set.find(elem) != valid_set.end());
        }
        uint64_t hA,hB;

        hA = get_rolling_hash(elem, 0);
        hB = get_rolling_hash(elem, 1);

        return contains(hA, hB);
    }

    inline int contains(uint64_t h0, uint64_t h1)
    { 
      if(fake){
        return (valid_hash0.find(h0) != valid_hash0.end()) 
          && (valid_hash1.find(h1) != valid_hash1.end());
      }
        uint64_t h = h0;
        for(int i=0; i<n_hash_func; i++, h = (h+h1)%tai)
        {
            if ((blooma[h >> 3 ] & bit_mask[h & 7]) != bit_mask[h & 7]){
                return 0;
            }
        }
        return 1;
        
    }


    /**********************************************************************************
    This is either very basic functions or for testing.
    ***********************************************************************************/

    //makes this a fake bloom filter that returns true only on specified kmers
    void fakify(std::set<bloom_elem> valid_kmers);    

    void setSeed(uint64_t seed) ;

    void set_number_of_hash_func(int i) ;
    
    /*void add(bloom_elem elem);
    int  contains(bloom_elem elem);
    void add(uint64_t hash0, uint64_t hash1);
    int contains(uint64_t hash0, uint64_t hash1);*/
    
    uint64_t tai;
    uint64_t nb_elem;
    
    void dump(char * filename);
    void load(char * filename);

    Bloom(uint64_t tai_bloom, int k);
    Bloom(int tai_bloom);
    Bloom(uint64_t tai_bloom);

    Bloom();
    
    ~Bloom();
};

void load_two_filters(Bloom* bloo1, Bloom* bloo2, std::string reads_filename);

#endif

