//
//  Bloom.cpp
//
//  Created by Guillaume Rizk on 9/02/12.
//

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <math.h>
#include "Bloom.h"
#include <set>
#include <list>
using std::ifstream;
using std::string;


/**********************************************************************************
These are the important things that are currently being used.
***********************************************************************************/


#include <iostream>
#include <cmath>
#include <algorithm>
#include <functional>

// root finder, taken from https://rosettacode.org/wiki/Roots_of_a_function#Brent.27s_Method
double brents_fun(std::function<double (double)> f, double lower, double upper, double tol, unsigned int max_iter)
{
    double a = lower;
    double b = upper;
    double fa = f(a);   // calculated now to save function calls
    double fb = f(b);   // calculated now to save function calls
    double fs = 0;      // initialize 
 
    if (!(fa * fb < 0))
    {
        std::cout << "Signs of f(lower_bound) and f(upper_bound) must be opposites" << std::endl; // throws exception if root isn't bracketed
        return -11;
    }
 
    if (std::abs(fa) < std::abs(b)) // if magnitude of f(lower_bound) is less than magnitude of f(upper_bound)
    {
        std::swap(a,b);
        std::swap(fa,fb);
    }
 
    double c = a;           // c now equals the largest magnitude of the lower and upper bounds
    double fc = fa;         // precompute function evalutation for point c by assigning it the same value as fa
    bool mflag = true;      // boolean flag used to evaluate if statement later on
    double s = 0;           // Our Root that will be returned
    double d = 0;           // Only used if mflag is unset (mflag == false)
 
    for (unsigned int iter = 1; iter < max_iter; ++iter)
    {
        // stop if converged on root or error is less than tolerance
        if (std::abs(b-a) < tol)
        {
            std::cout << "After " << iter << " iterations the root is: " << s << std::endl;
            return s;
        } // end if
 
        if (fa != fc && fb != fc)
        {
            // use inverse quadratic interopolation
            s =   ( a * fb * fc / ((fa - fb) * (fa - fc)) )
                + ( b * fa * fc / ((fb - fa) * (fb - fc)) )
                + ( c * fa * fb / ((fc - fa) * (fc - fb)) );
        }
        else
        {
            // secant method
            s = b - fb * (b - a) / (fb - fa);
        }
 
            // checks to see whether we can use the faster converging quadratic && secant methods or if we need to use bisection
        if (    ( (s < (3 * a + b) * 0.25) || (s > b) ) ||
                ( mflag && (std::abs(s-b) >= (std::abs(b-c) * 0.5)) ) ||
                ( !mflag && (std::abs(s-b) >= (std::abs(c-d) * 0.5)) ) ||
                ( mflag && (std::abs(b-c) < tol) ) ||
                ( !mflag && (std::abs(c-d) < tol))  )
        {
            // bisection method
            s = (a+b)*0.5;
 
            mflag = true;
        }
        else
        {
            mflag = false;
        }
 
        fs = f(s);  // calculate fs
        d = c;      // first time d is being used (wasnt used on first iteration because mflag was set)
        c = b;      // set c equal to upper bound
        fc = fb;    // set f(c) = f(b)
 
        if ( fa * fs < 0)   // fa and fs have opposite signs
        {
            b = s;
            fb = fs;    // set f(b) = f(s)
        }
        else
        {
            a = s;
            fa = fs;    // set f(a) = f(s)
        }
 
        if (std::abs(fa) < std::abs(fb)) // if magnitude of fa is less than magnitude of fb
        {
            std::swap(a,b);     // swap a and b
            std::swap(fa,fb);   // make sure f(a) and f(b) are correct after swap
        }
 
    } // end for
 
    std::cout<< "The solution does not converge or iterations are not sufficient" << std::endl;
 
} // end brents_fun


void Bloom::addPair(JuncPair pair){
    bloom_elem elem1 = pair.kmer1;
    bloom_elem elem2 = pair.kmer2;
    uint64_t hA,hB;

    hA = oldHash(elem1, 0) + oldHash(elem2, 0);
    hB = oldHash(elem1, 1) + oldHash(elem2, 1);

    //printf("Adding %lli, %lli\n", hA, hB);
    add(hA, hB);
}

int Bloom::containsPair(JuncPair pair){
    bloom_elem elem1 = pair.kmer1;
    bloom_elem elem2 = pair.kmer2;
    uint64_t hA,hB;

    hA = oldHash(elem1, 0) + oldHash(elem2, 0);
    hB = oldHash(elem1, 1) + oldHash(elem2, 1);

    //printf("Checking %lli, %lli\n", hA, hB);

    return contains(hA, hB);
}

void Bloom::fakify(){
    fake = true;
}

void Bloom::addFakeKmers(std::set<bloom_elem> valid_kmers) {
    valid_set.insert(valid_kmers.begin(), valid_kmers.end());
}


 Bloom::Bloom(uint64_t tai_bloom, int kVal)
 {
    fake = false;
     //printf("custom construc \n");
    k = kVal;
    n_hash_func = 4 ;//def
    user_seed =0;
    nb_elem = 0;
    hashSize = (int) log2(tai_bloom)+1;
    //printf("Hash size: %d \n", hashSize);
    tai = pow(2, hashSize);
    // tai = tai_bloom;
    //printf("Tai: %lli \n", tai);
    if(tai == 0){
       tai = 1;
    }
    bloomMask = tai-1;
    //printf("Mask: %lli \n", bloomMask);
    nchar = (tai/8LL);
    blooma =(unsigned char *)  malloc( nchar *sizeof(unsigned char)); // 1 bit per elem
    //printf("Allocation for filter: %lli bits. \n",nchar *sizeof(unsigned char)*8);
    memset(blooma,0,nchar *sizeof(unsigned char));
    //fprintf(stderr,"malloc bloom %lli MB \n",(tai/8LL)/1024LL/1024LL);
    this->generate_hash_seed();
 }

float Bloom::weight()
{
    // return the number of 1's in the Bloom, nibble by nibble
    const unsigned char oneBits[] = {0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4};
    long weight = 0;
    for(uint64_t index = 0; index < nchar; index++)
    {
        unsigned char current_char = blooma[index];
        weight += oneBits[current_char&0x0f];
        weight += oneBits[current_char>>4];
    }
    return (float)weight/(float)tai;
}


Bloom* Bloom::create_bloom_filter_2_hash(uint64_t estimated_items, float fpRate){
     
    Bloom * bloo1;
    int bits_per_item = 2*(int)(1/pow(fpRate,.5));// needed to process argv[5]

    printf("Bits per kmer: %d \n", bits_per_item);

    // int estimated_bloom_size = max( (int)ceilf(log2f(nb_reads * NBITS_PER_KMER )), 1);
    uint64_t estimated_bloom_size = (uint64_t) (estimated_items*bits_per_item);
    printf("Estimated items: %lli \n", estimated_items);
    printf("Estimated bloom size: %lli .\n", estimated_bloom_size);
    
    printf("BF memory: %f MB\n", (float)(estimated_bloom_size/8LL /1024LL)/1024);
    bloo1 = new Bloom(estimated_bloom_size, sizeKmer);


    printf("Number of hash functions: %d \n", 2);
    bloo1->set_number_of_hash_func(2);

    return bloo1;
}


Bloom* Bloom::create_bloom_filter_optimal(uint64_t estimated_items, float fpRate){
     
    Bloom * bloo1;
    int bits_per_item = -log(fpRate)/log(2)/log(2); // needed to process argv[5]

    printf("Bits per kmer: %d \n", bits_per_item);
    // int estimated_bloom_size = max( (int)ceilf(log2f(nb_reads * NBITS_PER_KMER )), 1);
    uint64_t estimated_bloom_size = (uint64_t) (estimated_items*bits_per_item);
    //printf("Estimated items: %lli \n", estimated_items);
    //printf("Estimated bloom size: %lli.\n", estimated_bloom_size);
    
    printf("BF memory: %f MB\n", (float)(estimated_bloom_size/8LL /1024LL)/1024);
    bloo1 = new Bloom(estimated_bloom_size, sizeKmer);

    printf("Number of hash functions: %d \n", (int)floorf(0.7*bits_per_item));
    bloo1->set_number_of_hash_func((int)floorf(0.7*bits_per_item));

    return bloo1;
}

void load_two_filters(Bloom* bloo1, Bloom* bloo2, string reads_filename, bool fastq){
    ifstream solidReads;
    solidReads.open(reads_filename);

    int readsProcessed = 0;

    string read;
    time_t start, stop;
    time(&start);
    printf("Weights before load: %f, %f \n", bloo1->weight(), bloo2->weight());
    uint64_t hashA, hashB;
    kmer_type canonKmer;
    uint64_t unambiguousReads = 0;
    while (getline(solidReads, read))
    {
        getline(solidReads, read);//since it's a fasta or fastq we skip the first of every pair of lines
        std::list<string> readList = getUnambiguousReads(read);
        while(!readList.empty()){
            read = readList.front();
            readList.pop_front();
            unambiguousReads++;
            for(ReadKmer kmer = ReadKmer(&read); kmer.getDistToEnd() >= 0 ; kmer.forward(), kmer.forward()){
                canonKmer = kmer.getCanon();
                hashA = bloo1->oldHash(canonKmer, 0);
                hashB = bloo1->oldHash(canonKmer, 1);
                if(bloo1->contains(hashA, hashB)){
                    bloo2->add(hashA, hashB);
                }
                else{
                    bloo1->add(hashA, hashB);
                }
            }    
        }
        readsProcessed++;
        if ((readsProcessed%100000)==0) fprintf (stdout,"reads consumed: %c %lld\n",13,(long long)readsProcessed);
        
        if(fastq) getline(solidReads, read),getline(solidReads, read); //ignore two more lines if its fastq
    }

    solidReads.close();
    printf("\n");
    printf("Weights after load: %f, %f \n", bloo1->weight(), bloo2->weight());
    printf("Reads processed: %d\n", readsProcessed);
    printf("Unambiguous reads: %lli\n", unambiguousReads);
    time(&stop);
    printf("Time to load: %f \n", difftime(stop,start));
}

void load_single_filter(Bloom* bloo1, string reads_filename, bool fastq){
    ifstream solidReads;
    solidReads.open(reads_filename);
    int readsProcessed = 0;
    string read;
    time_t start, stop;
    time(&start);
    printf("Weight before load: %f\n", bloo1->weight());
    uint64_t hashA, hashB;
    kmer_type canonKmer;
    uint64_t unambiguousReads = 0;

    while (getline(solidReads, read))
    {
        getline(solidReads, read);//since it's a fasta or fastq we skip the first of every pair of lines
        std::list<string> readList = getUnambiguousReads(read);
        while(!readList.empty()){
            read = readList.front();
            readList.pop_front();
            unambiguousReads++;
            for(ReadKmer kmer = ReadKmer(&read); kmer.getDistToEnd() >= 0 ; kmer.forward(), kmer.forward()){
                canonKmer = kmer.getCanon();
                hashA = bloo1->oldHash(canonKmer, 0);
                hashB = bloo1->oldHash(canonKmer, 1);                                
                bloo1->add(hashA, hashB);                
            }    
        }
        readsProcessed++;
        if ((readsProcessed%10000)==0) fprintf (stderr,"%c %lld",13,(long long)readsProcessed);        
        if(fastq) getline(solidReads, read),getline(solidReads, read); //ignore two more lines if its fastq
    }
    solidReads.close();
    printf("\n");
    printf("Weight after load: %f\n", bloo1->weight());
    printf("Reads processed: %d\n", readsProcessed);
    printf("Unambiguous reads: %lli\n", unambiguousReads);
    time(&stop);
    printf("Time to load: %f \n", difftime(stop,start));
}

void Bloom::load_from_reads(const char* reads_filename){
    ifstream solidReads;
    solidReads.open(reads_filename);

    int readsProcessed = 0;

    //should fix this.. just because right now there's no empty constructor
    string fake = "FAKE";
    ReadKmer kmer(&fake);

    string read;
    time_t start, stop;
    time(&start);
    printf("Weight before load: %f \n", weight());
    while (getline(solidReads, read))
    {
        getline(solidReads, read);//since it's a fasta we skip the first of every pair of lines
        for(kmer = ReadKmer(&read); kmer.getDistToEnd() >= 0 ; kmer.forward(), kmer.forward()){
            oldAdd(kmer.getCanon());
        }    
        readsProcessed++;
        if ((readsProcessed%10000)==0) fprintf (stderr,"%c %lld",13,(long long)readsProcessed);
    }

    solidReads.close();
    printf("\n");
    printf("Weight after load: %f \n", weight()); 
    time(&stop);
    printf("Time to load: %f \n", difftime(stop,start));
}

void Bloom::load_from_kmers(const char* kmers_filename){
    ifstream solidKmers;
    solidKmers.open(kmers_filename);

    int kmersProcessed = 0;
    kmer_type left,right;
    string kpomer;

    printf("Weight before load: %f \n", weight());
    int badKmerCount = 0;
    while (getline(solidKmers, kpomer))
    {
        if(kpomer.length() != sizeKmer+1){
            badKmerCount++;
            continue;
        }
      // printf("kpomer %s \n", &kpomer[0]);
        getFirstKmerFromRead(&left,&kpomer[0]);
        right = next_kmer(left, NT2int(kpomer[sizeKmer]),FORWARD);
        //printf("left %s \n", print_kmer(left));
        //printf("right %s \n", print_kmer(right));
        oldAdd(get_canon(left));
        oldAdd(get_canon(right));
        kmersProcessed++;
        if ((kmersProcessed%10000)==0) fprintf (stderr,"%c %lld",13,(long long)kmersProcessed);
    }
    printf("\n");
    solidKmers.close();
    printf("Number of wrong length kmers: %d \n", badKmerCount);
    printf("Weight after load: %f \n", weight()); 
}

int Bloom::getNumHash(){
    return n_hash_func;
}

int Bloom::getHashSize(){
    return hashSize;
}

uint64_t Bloom::getBloomMask(){
    return bloomMask;
}

/**********************************************************************************
Most of the below is not currently used, or has to do with generating hash seeds.  
Also much of it is for incremental hashing.
***********************************************************************************/

Bloom::Bloom()
{
    //empty default constructor
    nb_elem = 0;
    blooma = NULL;
}

void Bloom::setSeed(uint64_t seed)
{
    if(user_seed==0)
    {
        user_seed = seed;
        this->generate_hash_seed(); //regenerate the hash with the new seed
    }
    else{
        fprintf(stderr,"Warning! you should not change the seed a second time!, resuming with previous seed %llu \n",(unsigned long long)user_seed);
    }
}

void Bloom::set_number_of_hash_func(int i)
{
    if(i>NSEEDSBLOOM || i<1){
        fprintf(stderr,"%i is not a valid value for number of hash funcs, should be in [1-%i], resuming wild old value %i\n",i,NSEEDSBLOOM,n_hash_func );
        return;
    }  
    n_hash_func = i;
}

void Bloom::generate_hash_seed()
{
    unsigned int i;
    for ( i = 0; i < NSEEDSBLOOM; ++i)
    {
        seed_tab[i]= rbase[i];
    }
    for ( i = 0; i < NSEEDSBLOOM; ++i)
    {
        seed_tab[i]= seed_tab[i] * seed_tab[(i+3) % NSEEDSBLOOM] + user_seed ;
        //printf("%lli \n", seed_tab[i]);
    }

    for ( i = 0; i < 4; ++i)
    {
        char_hash[0][i]= rbase[i];
    }
    for ( i = 0; i < 4; ++i)
    {
         char_hash[0][i]=  char_hash[0][i] *  char_hash[0][(i+3) % 4] + user_seed ;
    }
    for ( i = 0; i < 4; ++i)
    {
        char_hash[0][i] &= bloomMask;
       // printf("%lli \n", char_hash[0][i]); //seems random!
    }

    for ( i = 0; i < 4; ++i)
    {
        char_hash[1][i]= rbase[i+4];
    }
    for ( i = 0; i < 4; ++i)
    {
         char_hash[1][i]=  char_hash[1][i] *  char_hash[1][(i+3) % 4] + user_seed ;
    }
    for ( i = 0; i < 4; ++i)
    {
        char_hash[1][i] &= bloomMask;
       // printf("%lli \n", char_hash[1][i]); //seems random
    }
}


uint64_t Bloom::getCharHash(int key, int num_hash){
    return char_hash[num_hash][key];
}

uint64_t Bloom::getLastCharHash(uint64_t key, int num_hash){
    return char_hash[num_hash][(int)(key & 3)];
}

//only for num_hash = 0 or 1
uint64_t Bloom::get_rolling_hash(uint64_t key, int num_hash)
{
    uint64_t hash = getLastCharHash(key, num_hash);
    for(int i = 1; i < k; i++){
        hash = rotate_right(hash, 1);
        key >>= 2;
        hash ^= getLastCharHash(key, num_hash);
    }
    return hash;
}


Bloom::~Bloom()
{
    valid_set.clear();
  if(blooma!=NULL) 
    free(blooma);
}

void Bloom::dump(char * filename)
{
 FILE *file_data;
 file_data = fopen(filename,"wb");
 fwrite(blooma, sizeof(unsigned char), nchar, file_data); //1+
 printf("bloom dumped \n");

}

void Bloom::load(char * filename)
{
 FILE *file_data;
 file_data = fopen(filename,"rb");
 printf("loading bloom filter from file, nelem %lli \n",nchar);
 int a = fread(blooma, sizeof(unsigned char), nchar, file_data);// go away warning..
 printf("bloom loaded\n");
}
