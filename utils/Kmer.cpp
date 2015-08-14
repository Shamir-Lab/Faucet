#ifndef ASSERTS
#define NDEBUG // disable asserts, they're computationnally intensive
#endif 

#include <stdio.h>
#include <assert.h>
#include <algorithm> // for min
#include "Kmer.h"
#include "lut.h"
#include <string>
#include <list>

using namespace std;

int sizeKmer;
uint64_t nsolids = 0;
kmer_type kmerMask;    
kmer_type kmerMaskm1;
const bool FORWARD = true;
const bool BACKWARD = false;
char bin2NT[4] = {'A','C','T','G'};

char getNucChar(int nucIndex){
    return bin2NT[nucIndex];
}

void setSizeKmer(int k){
    sizeKmer = k;
    if (sizeKmer>((int)sizeof(kmer_type)*4))
    {
        printf("Max kmer size on this compiled version is %u\n",sizeof(kmer_type)*4);
        exit(1);
    }
    if (sizeKmer == (int)(sizeof(kmer_type)*4))
        kmerMask = -1;
    else
        kmerMask=(((kmer_type)1)<<(sizeKmer*2))-1;
}

bool isValidNuc(char nt){
     switch(nt){
        case 'A':
        case 'C': 
        case 'T':
        case 'G':
            return true;
        default:
            return false;
    }
}

//returns a list of valid nucleotide sequences of length at least k
//Invalid chars and valid sequences shorter than k are discarded
std::list<std::string> getUnambiguousReads(string read){
    std::list<std::string> resultList;
    std::string result;
    int index = 0;
    while(index < read.length()){
        result = "";
        while(!isValidNuc(read[index]) && index < read.length()){
            index++;
        }
        while(isValidNuc(read[index]) && index < read.length()){
            result += read[index];
            index++;
        }
        if(result.length() >= sizeKmer) resultList.push_front(result);
    }
    return resultList;
}

int NT2int(char nt)
{
    int i;
    i = nt;
    i = (i>>1)&3; // that's quite clever, guillaume.
    return i;
}

int revcomp_int(int nt_int)
{   
    return (nt_int<2)?nt_int+2:nt_int-2;
}


unsigned char  code4NT(char *seq)
{
    int i;
    unsigned char x;
    
    x=0;
    for (i=0; i<4; ++i)
    {
        x = x*4 + NT2int(seq[i]);
    }
    return x;
}

unsigned char  code_n_NT(char *seq, int nb)
{
    int i;
    unsigned char x;
    
    x=0;
    for (i=0; i<nb; ++i)
    {
        x = x*4 + NT2int(seq[i]);
    }
    x = x << ((4-nb)*2) ;
    return x;
}


kmer_type  codeSeed(char *seq)
{
    return codeSeed(seq, sizeKmer, kmerMask);
}

kmer_type  codeSeed(char *seq, int sizeKmer, kmer_type kmerMask)
{
    int i;
    kmer_type x;
    
    x=0;
    for (i=0; i<sizeKmer; ++i)
    {
        x = x*4 + NT2int(seq[i]);
    }
    return x;
}
kmer_type  codeSeedRight(char *seq, kmer_type  val_seed, bool new_read)
{
    return codeSeedRight(seq, val_seed, new_read, sizeKmer, kmerMask);
}

kmer_type  codeSeedRight(char *seq, kmer_type  val_seed, bool new_read, int sizeKmer, kmer_type kmerMask)
{
    if (new_read) return codeSeed(seq, sizeKmer, kmerMask);
    // shortcut    
    return (val_seed * 4 + NT2int(seq[sizeKmer-1])) & kmerMask ;
}

kmer_type  codeSeedRight_revcomp(char *seq, kmer_type  val_seed, bool new_read)
{
    return codeSeedRight_revcomp(seq, val_seed, new_read, sizeKmer, kmerMask);    
}

kmer_type  codeSeedRight_revcomp(char *seq, kmer_type  val_seed, bool new_read, int sizeKmer, kmer_type kmerMask)
{
    
  if (new_read) return revcomp( codeSeed(seq, sizeKmer, kmerMask), sizeKmer );
  // shortcut    
  return ((val_seed >> 2) +  ( ((kmer_type) comp_NT[NT2int(seq[sizeKmer-1])]) <<  (2*(sizeKmer-1))  )  ) & kmerMask;
    
}

// warning: only call this function for sequential enumeration of kmers (no arbitrary position)
kmer_type extractKmerFromRead(char *readSeq, int position, kmer_type *graine, kmer_type *graine_revcomp)
{
    return extractKmerFromRead(readSeq, position, graine, graine_revcomp, true);
}

kmer_type extractKmerFromRead(char *readSeq, int position, kmer_type *graine, kmer_type *graine_revcomp, bool sequential)
{
    return extractKmerFromRead(readSeq, position, graine, graine_revcomp, sequential, sizeKmer, kmerMask);
}

kmer_type extractKmerFromRead(char *readSeq, int position, kmer_type *graine, kmer_type *graine_revcomp, bool sequential, int sizeKmer, kmer_type kmerMask)
{
    assert(graine != graine_revcomp); // make sure two different pointers
    bool new_read = (position == 0) || (!sequential); // faster computation for immediately overlapping kmers

        *graine = codeSeedRight(&readSeq[position], *graine, new_read, sizeKmer, kmerMask);
        *graine_revcomp = codeSeedRight_revcomp(&readSeq[position], *graine_revcomp, new_read, sizeKmer, kmerMask);

    return  min(*graine,*graine_revcomp);
}


//return the rightmost nucleotide of the kmer
int first_nucleotide(kmer_type kmer)
{
    int result;
#ifdef _largeint
    LargeInt<KMER_PRECISION> t = kmer;
    result = t.toInt()&3;
#else
#ifdef _ttmath
    ttmath::UInt<KMER_PRECISION> t = kmer&3;
    t.ToInt(result);
#else
    result = kmer&3;
#endif
#endif
    return result;
}

int code2seq (kmer_type code, char *seq)
{
    return code2seq (code, seq, sizeKmer, kmerMask);
}

int code2seq (kmer_type code, char *seq, int sizeKmer, kmer_type kmerMask)
{
    int i;
    kmer_type temp = code;
    char bin2NT[4] = {'A','C','T','G'};

    for (i=sizeKmer-1; i>=0; i--)
    {
        seq[i]=bin2NT[first_nucleotide(temp&3)];
        temp = temp>>2;
    }
    //printf("sizeKmer = %d\n", sizeKmer);
    seq[sizeKmer]='\0';
    return sizeKmer;
}

// return the i-th nucleotide of the kmer_type kmer
// i = 0 corresponds to the leftmost nucleotide
int code2nucleotide( kmer_type code, int which_nucleotide)
{
    kmer_type temp = code;
    temp = temp >> (2*(sizeKmer-1-which_nucleotide));
    return first_nucleotide(temp&3);
}

uint64_t revcomp(uint64_t x, int size) {
  int i;

  uint64_t revcomp = x;
  //  printf("x %x    revcomp  %x \n",x,revcomp);
  unsigned char * kmerrev  = (unsigned char *) (&revcomp);
  unsigned char * kmer  = (unsigned char *) (&x);

  for (i=0; i<8; ++i)
    {
       kmerrev[7-i] =     revcomp_4NT[kmer[i]];
    }

    return (revcomp >> (2*( 4*sizeof(uint64_t) - size))  ) ;
}

uint64_t revcomp(uint64_t x) {
  return revcomp(x,sizeKmer);
}

#ifdef _largeint
LargeInt<KMER_PRECISION> revcomp(LargeInt<KMER_PRECISION> x, int size) {
    int i;

    kmer_type revcomp = x;
    //  printf("x %x    revcomp  %x \n",x,revcomp);
    unsigned char * kmerrev  = (unsigned char *) (&(revcomp.array[0]));
    unsigned char * kmer  = (unsigned char *) (&(x.array[0]));

    for (i=0; i<8*KMER_PRECISION; ++i)
    {
        kmerrev[8*KMER_PRECISION-1-i] =     revcomp_4NT[kmer[i]];
    }

    return (revcomp >> (2*( 32*KMER_PRECISION - size))  ) ;
}

LargeInt<KMER_PRECISION> revcomp(LargeInt<KMER_PRECISION> x) {
  return revcomp(x,sizeKmer);
}
#endif

#ifdef _ttmath
ttmath::UInt<KMER_PRECISION> revcomp(ttmath::UInt<KMER_PRECISION> x, int size) {
    int i;

    kmer_type revcomp = x;
    //  printf("x %x    revcomp  %x \n",x,revcomp);
    unsigned char * kmerrev  = (unsigned char *) (&revcomp);
    unsigned char * kmer  = (unsigned char *) (&x);

    for (i=0; i<4*KMER_PRECISION; ++i)
    {
        kmerrev[4*KMER_PRECISION-1-i] =     revcomp_4NT[kmer[i]];
    }

    return (revcomp >> (2*( 16*KMER_PRECISION - size))  ) ;
}

ttmath::UInt<KMER_PRECISION> revcomp(ttmath::UInt<KMER_PRECISION> x) {
  return revcomp(x,sizeKmer);
}
#endif

#ifdef _LP64
__uint128_t revcomp(__uint128_t x, int size)
{
    //                  ---64bits--   ---64bits--
    // original kmer: [__high_nucl__|__low_nucl___]
    //
    // ex:            [         AC  | .......TG   ]
    //
    //revcomp:        [         CA  | .......GT   ]
    //                 \_low_nucl__/\high_nucl/
    uint64_t high_nucl = (uint64_t)(x>>64);
    int nb_high_nucl = size>32?size - 32:0;
    __uint128_t revcomp_high_nucl = revcomp(high_nucl, nb_high_nucl);
    if (size<=32) revcomp_high_nucl = 0; // srsly dunno why this is needed. gcc bug? uint64_t x ---> (x>>64) != 0
    
    uint64_t low_nucl = (uint64_t)(x&((((__uint128_t)1)<<64)-1));
    int nb_low_nucl = size>32?32:size;
    __uint128_t revcomp_low_nucl = revcomp(low_nucl, nb_low_nucl);

    return (revcomp_low_nucl<<(2*nb_high_nucl)) + revcomp_high_nucl;
}	

__uint128_t revcomp(__uint128_t x) {
  return revcomp(x,sizeKmer);
}
#endif

char revcomp_char(char c){
    switch(c){
        case 'A':
            return 'T';
        case 'C': 
            return 'G';
        case 'T':
            return 'A';
        case 'G':
            return 'C';
        default:
            printf("ERROR: tried to revcomp an invalid char\n");
            return 'X';
    }
}

std::string canon_contig(std::string contig){
    std::string revcomp_contig = revcomp_string(contig);
    if(contig <= revcomp_contig){
        return contig;
    }
    return revcomp_contig;
}

std::string revcomp_string(std::string s){
    string revcomp = "";
    for(int i = 0; i < s.length(); i++){
        revcomp += revcomp_char(s[s.length()-1-i]);
    }
    return revcomp;
}

// will be used by assemble()
void revcomp_sequence(char * s, int len)
{
#define CHAR_REVCOMP(a,b) {switch(a){\
	case 'A': b='T';break;case 'C': b='G';break;case 'G': b='C';break;case 'T': b='A';break;default: b=a;break;}}
		  int i;
		  unsigned char t;
		  for (i=0;i<len/2;i++)
		  {
			  t=s[i];
			  CHAR_REVCOMP(s[len-i-1],s[i]);
			  CHAR_REVCOMP(t,s[len-i-1]);
		  }
		  if (len%2==1)
			  CHAR_REVCOMP(s[len/2],s[len/2]);

}

//DONT USE THIS
//Legacy from initial minia code: use the other versioN!
//this one returns the canonical version which just confuses repeated extensions
kmer_type next_kmer(kmer_type graine, int added_nt, int *strand)
{
    printf("Using legacy next_kmer!");
    assert(added_nt<4);
    assert(graine<=revcomp(graine));
    assert((strand == NULL) || (*strand<2));

    kmer_type new_graine;
    kmer_type temp_graine;

    if (strand != NULL && *strand == 1)// the kmer we're extending is actually a revcomp sequence in the bidirected debruijn graph node
        temp_graine = revcomp(graine);
    else
        temp_graine = graine;

    new_graine = (((temp_graine) * 4 )  + added_nt) & kmerMask;
    //new_graine = (((graine) >> 2 )  + ( ((kmer_type)added_nt) << ((sizeKmer-1)*2)) ) & kmerMask; // previous kmer
    kmer_type revcomp_new_graine = revcomp(new_graine);

    if (strand != NULL)
        *strand = (new_graine < revcomp_new_graine)?0:1;

    return min(new_graine,revcomp_new_graine);
}


//tested!
//Strand tell is to shift right (strand = 0) or left (strand = 1)
void shift_kmer(kmer_type* graine, int added_nt, bool direction){

    assert(added_nt<4);
    assert(strand<2);

    if (direction == BACKWARD){// the kmer we're extending is actually a revcomp sequence in the bidirected debruijn graph node
        *graine >>= 2 ;
        *graine += (((uint64_t)added_nt) << (2*sizeKmer-2));
        *graine &= kmerMask;
    }
    else {
        *graine <<= 2;
        *graine += added_nt;
        *graine &= kmerMask;
    }
}


//Tested!
void getFirstKmerFromRead(kmer_type *kmer, char* read){
      for(int i = 0; i < sizeKmer; i++){
        shift_kmer(kmer, NT2int(read[i]), FORWARD);
      }
}

kmer_type getKmerFromRead(std::string read, int index){
    kmer_type kmer(0);
    for(int i = 0; i < sizeKmer; i++){
        shift_kmer(&kmer, NT2int(read[i+index]), FORWARD);
    }
    return kmer;
}

//Tested!
//USE THIS ONE
//preserves the orientation of the kmer, considers it as a revcomp if strand == 1, or normal if strand == 0
kmer_type next_kmer(kmer_type graine, int added_nt,  bool direction)
{
    kmer_type new_graine = graine;
    shift_kmer(&new_graine, added_nt, direction);
    return new_graine;
}

//Tested!
//Assumes index_in_read is in fact the index of kmer in the read!
kmer_type next_kmer_in_read(kmer_type kmer, int index_in_read, char* read, int strand){
    int nucPos = index_in_read + sizeKmer;
    if(strand == 1){
        nucPos = index_in_read-1;
    }
    return next_kmer(kmer,NT2int(read[nucPos]), strand);
}

kmer_type maskKmer(kmer_type kmer){
    return kmer & kmerMask;
}
// void revcompread(char** read, int length){
//     char* firstChar = *read;
//     char temp;
//     for (int i = 0; i <= length/2; i++){
//         temp = firstChar[i];
//         firstChar[i] = firstChar[length-1-i];
//         firstChar[length-1-i] = temp;
//     }
// }


//////////////////////////funcs for binary readSeq


kmer_type  codeSeed_bin(char *seq)
{
    int i;
    kmer_type x;
    
    x=0;
    for (i=0; i<sizeKmer; ++i)
    {
        x = x*4 + (seq[i]);
    }
    return x;
}

inline kmer_type  codeSeedRight_bin(char *seq, kmer_type  val_seed, bool new_read)
{
    
    if (new_read) return codeSeed_bin(seq);
    // shortcut
    return (val_seed * 4 + (seq[sizeKmer-1])) & kmerMask ;
    
}

inline kmer_type  codeSeedRight_revcomp_bin(char *seq, kmer_type  val_seed, bool new_read)
{
    
    if (new_read) return revcomp(codeSeed_bin(seq));
    // shortcut
  return ((val_seed >> 2) +  ( ((kmer_type) comp_NT[(int)(seq[sizeKmer-1])]) <<  (2*(sizeKmer-1))  )  ) & kmerMask;
    
}



kmer_type extractKmerFromRead_bin(char *readSeq, int position, kmer_type *graine, kmer_type *graine_revcomp, bool use_compressed)
{
    assert(graine != graine_revcomp); // make sure two different pointers

    bool new_read = (position == 0);
    if(!use_compressed)
    {
        *graine = codeSeedRight(&readSeq[position], *graine, new_read);
        *graine_revcomp = codeSeedRight_revcomp(&readSeq[position], *graine_revcomp, new_read);
    }
    else
    {
        *graine = codeSeedRight_bin(&readSeq[position], *graine, new_read);
        *graine_revcomp = codeSeedRight_revcomp_bin(&readSeq[position], *graine_revcomp, new_read);
    }
    return  min(*graine,*graine_revcomp);
}

kmer_type get_canon(kmer_type kmer){
    return min(kmer, revcomp(kmer));
}

kmer_type get_kmer_from_read(char* read, int pos){
    kmer_type kmer = 0;
    for (int i =pos; i < pos + sizeKmer; i++){
        shift_kmer(&kmer, NT2int(read[i]), FORWARD);
    }
    return kmer;
}

//Tested!
kmer_type advance_kmer(char* read, kmer_type* kmer,  int startPos, int endPos){
    if(endPos < startPos){
        printf("advance_kmer does not yet support backwards advances!\n");
    }
    for (int i = startPos + sizeKmer; i < endPos + sizeKmer; i++){
        shift_kmer(kmer, NT2int(read[i]), FORWARD);
    }
}

// debug only: convert a kmer_type to char*
char debug_kmer_buffer[1024];
char* print_kmer(kmer_type kmer)
{
    return print_kmer(kmer,sizeKmer,kmerMask);
}

char* print_kmer(kmer_type kmer, int sizeKmer, kmer_type kmerMask)
{
    code2seq(kmer,debug_kmer_buffer, sizeKmer, kmerMask);
    return debug_kmer_buffer;
}