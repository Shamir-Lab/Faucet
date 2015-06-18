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
#include "Kmer.h"
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


static const int cpt_per_char = 2;    
static const unsigned char cpt_mask[cpt_per_char] = {
    0x0F,  //00001111
    0xF0,  //11110000
};

static const uint64_t cpt_mask21[21] = {
  0x0000000000000007ULL,//00000....00000111
  0x0000000000000038ULL,
  0x00000000000001C0ULL,
  0x0000000000000E00ULL,
  0x0000000000007000ULL,
  0x0000000000038000ULL,
  0x00000000001C0000ULL,
  0x0000000000E00000ULL,
  0x0000000007000000ULL,
  0x0000000038000000ULL,
  0x00000001C0000000ULL,
  0x0000000E00000000ULL,
  0x0000007000000000ULL,
  0x0000038000000000ULL,
  0x00001C0000000000ULL,
  0x0000E00000000000ULL,
  0x0007000000000000ULL,
  0x0038000000000000ULL,
  0x01C0000000000000ULL,
  0x0E00000000000000ULL,
  0x7000000000000000ULL

};


static const uint64_t cpt_mask32[32] = {
  0x0000000000000003ULL,//00000....00000011
  0x000000000000000CULL,
  0x0000000000000030ULL,//00000....000110000
  0x00000000000000C0ULL,
  0x0000000000000300ULL,
  0x0000000000000C00ULL,
  0x0000000000003000ULL,
  0x000000000000C000ULL,
  0x0000000000030000ULL,
  0x00000000000C0000ULL,
  0x0000000000300000ULL,
  0x0000000000C00000ULL,
  0x0000000003000000ULL,
  0x000000000C000000ULL,
  0x0000000030000000ULL,
  0x00000000C0000000ULL,
  0x0000000300000000ULL,
  0x0000000C00000000ULL,
  0x0000003000000000ULL,
  0x000000C000000000ULL,
  0x0000030000000000ULL,
  0x00000C0000000000ULL,
  0x0000300000000000ULL,
  0x0000C00000000000ULL,
  0x0003000000000000ULL,
  0x000C000000000000ULL,
  0x0030000000000000ULL,
  0x00C0000000000000ULL,
  0x0300000000000000ULL,
  0x0C00000000000000ULL,
  0x3000000000000000ULL,
  0xC000000000000000ULL
};



/* static const unsigned char incr_cpt_table[2][255] = 
 { 
 {1, 2,3}, 
 {3, 4,3}, 
 }; 
 */

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


/*
 
 0x2E7E5A8996F99AA5, 
 0x74B2E1FB222EFD24,
 0x8BBE030F6704DC29,
 0x6D8FD7E91C11A014, 
 0xFC77642FF9C4CE8C, 
 0x318FA6E7C040D23D, 
 0xF874B1720CF914D5, 
 0xC569F575CDB2A091, 
 */

//static uint64_t pri1=0x5AF3107A401FULL;
//static uint64_t pri2 =0x78C27CE77ULL;

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

    long weight();

    Bloom(uint64_t tai_bloom, int k);
    Bloom(int tai_bloom);
    Bloom(uint64_t tai_bloom);

    Bloom();
    
    
    ~Bloom();
    
    Bloom* create_bloom_filter_2_hash(uint64_t estimated_items, float fpRate);

    Bloom* create_bloom_filter_optimal(uint64_t estimated_items, float fpRate);

    void load_from_reads(const char* reads_filename);

    void load_from_kmers(const char* kmers_filename);
};



class BloomCpt: public Bloom { 
    public : 
    BloomCpt(int tai_bloom); 
    BloomCpt();

    void add(bloom_elem elem);
    
    int  contains_n_occ(bloom_elem elem, int nks);

    
}; 


class BloomCpt3: public BloomCpt { 
    public : 
    BloomCpt3(int tai_bloom); 
        ~BloomCpt3();

    uint64_t * blooma3;
    void add(bloom_elem elem);
    
    int  contains_n_occ(bloom_elem elem, int nks);
    
}; 



class BloomCpt2: public BloomCpt { 
    public : 
    BloomCpt2(int tai_bloom); 
    ~BloomCpt2();

    uint64_t * blooma2;
    void add(bloom_elem elem);
    
    int  contains_n_occ(bloom_elem elem, int nks);
    
}; 




#endif

