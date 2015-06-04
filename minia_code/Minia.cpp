#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <inttypes.h>
#include <stdint.h>
#include <algorithm> // for max/min
#include <vector> // for sorting_kmers
#include <sys/time.h>

#define NNKS 4 // default minimal abundance for solidity
#define MIN_CONTIG_SIZE (2*sizeKmer+1)


float fpRate = .01;
float max_memory; // the most memory one should alloc at any time, in MB

int order = 0; // deblooming order; 0 = debloom everything; 1 = don't debloom 1-node tips (experimental, untested, shouldn't work);// (made extern int in Traversal.h)
int read_length;

#include "Bank.h"
#include "Hash16.h"
#include "Set.h"
#include "Pool.h"
#include "Bloom.h"
#include "Debloom.h"
#include "Utils.h"
#include "SortingCount.h"
#include "Terminator.h"
#include "Kmer.h"
#include "Traversal.h"
#include "rvalues.h" // for 4bloom
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int64_t nb_reads;

/*
To run the new version, type make in the directory to compile.

Then, type ./minia 1 2 3 4 5 6, where
1 = name of reads file (current format is each line has a string of characters representing the read)
2 = k
3 = read length
4 = number of reads
5 = prefix for output files ("" works for me, I'll look into what exactly this does)
6 = false positive rate (.01 is a good base rate) 

This will load a bloom filter with all the kmers from the reads, then scan through them inserting potential false positives.
No error correction corrently, and no assembly.
*/

inline int handle_arguments(int argc, char *argv[]){
if(argc <  7)
    {
        fprintf (stderr,"usage:\n");
        fprintf (stderr," %s input_file kmer_size min_abundance estimated_nb_reads prefix\n",argv[0]);
        fprintf (stderr,"hints:\n min_abundance ~ 3\n estimated_nb_reads is in bp, does not need to be accurate, only controls memory usage\n prefix is any name you want the results to start with\n");

        return 1;
    }

    //1st arg: read file name
    strcpy(solid_reads_file,argv[1]);
    printf("Reads file name: %s \n", solid_reads_file);

    // 2rd arg: kmer size.
    sizeKmer=27;
    sizeKmer = atoi(argv[2]);
    if (sizeKmer%2==0)
    {
        sizeKmer-=1;
        printf("Need odd kmer size to avoid palindromes. I've set kmer size to %d.\n",sizeKmer);
    }
    if (sizeKmer>((int)sizeof(kmer_type)*4))
    {
        printf("Max kmer size on this compiled version is %lu\n",sizeof(kmer_type)*4);
        exit(1);
    }
    if (sizeKmer == (int)(sizeof(kmer_type)*4))
        kmerMask = -1;
    else
        kmerMask=(((kmer_type)1)<<(sizeKmer*2))-1;

    //3rd arg: read length
    read_length = atoi(argv[3]);

    //4th arg: number of reads
    nb_reads  = atoll(argv[4]);
    printf("Number reads: %d .\n", nb_reads);
    
    //5th arg: output prefix
    strcpy(prefix,argv[5]);

    //6th arg: false posiive rate
    fpRate = atof(argv[6]);
}

inline void load_bloom_filter(Bloom* bloo1, const char* reads_filename){

    ifstream solidReads;
    solidReads.open(reads_filename);

    int readsProcessed = 0;
    kmer_type new_graine, kmer;
    string read;

    printf("Weight before load: %ld \n", bloo1->weight());
    while (getline(solidReads, read))
    {
        kmer_type * kmer_table_seq =  (kmer_type*) malloc (sizeof(kmer_type)*read.length()*2);
        compute_kmer_table_from_one_seq(read.length(),&read[0],kmer_table_seq);

        for (int i = 0; i <= read.length() - sizeKmer ; i++){
            kmer = kmer_table_seq[i];
            bloo1->add(get_canon(kmer));
            readsProcessed++;
            if ((readsProcessed%10000)==0) fprintf (stderr,"%c %lld",13,(long long)readsProcessed);
        }

        free(kmer_table_seq);
    }

    solidReads.close();
    printf("\n");
    printf("Weight after load: %ld \n", bloo1->weight());
}

inline Bloom* create_bloom_filter(int estimated_items, float fpRate){
     
    Bloom * bloo1;
    int bits_per_item = -log(fpRate)/log(2)/log(2); // needed to process argv[5]

    // int estimated_bloom_size = max( (int)ceilf(log2f(nb_reads * NBITS_PER_KMER )), 1);
    uint64_t estimated_bloom_size = (uint64_t) (estimated_items*bits_per_item);
    printf("Estimated bloom size: %d .\n", (int)estimated_bloom_size);
    
    max_memory =  (float)(estimated_bloom_size/8LL /1024LL)/1024;

    printf("estimated values: nbits Bloom %lli. \n",estimated_bloom_size);

    printf("Max Memory: %f MB\n", max_memory);
    bloo1 = new Bloom(estimated_bloom_size);

    printf("Bits per kmer: %d \n", bits_per_item);

    printf("Number of hash functions: %d \n", (int)floorf(0.7*NBITS_PER_KMER));
    bloo1->set_number_of_hash_func((int)floorf(0.7*NBITS_PER_KMER));

return bloo1;
}

int main(int argc, char *argv[])
{
    
    if (handle_arguments(argc, argv) == 1){
        return 1;
    }

    int estimated_kmers = nb_reads*(read_length-sizeKmer);
    Bloom* bloo1 = create_bloom_filter(estimated_kmers, fpRate);
   
    load_bloom_filter(bloo1, solid_reads_file);

    // debloom, write false positives to disk, insert them into false_positives
    debloom(order, max_memory, bloo1);

    printf("Program reached end. \n");
    return 0;
}





inline void assemble(Bloom* bloo1)
{

    //////-------------------------------------------------------------------------------------------
    fprintf (stderr,"______________________________________________________ \n");
    fprintf (stderr,"___________ Assemble from bloom filter _______________ \n");
    fprintf (stderr,"______________________________________________________ \n\n");

    //////-------------------------------------------------------------------------------------------


    long long len_left = 0;
    long long len_right = 0;
    long long contig_len =0;
    long long maxlen=10000000;

    char *left_traversal  = (char *) malloc(maxlen*sizeof(char));
    char *right_traversal = (char *) malloc(maxlen*sizeof(char));
    char *contig          = (char *) malloc(2*(maxlen+sizeKmer)*sizeof(char));
    kmer_type kmer;

    long long nbContig =0;
    long long nbSmallContig =0;
    long long totalnt=0;
    long long max_contig_len=0;
    long long mlenleft=0,mlenright=0;
    int64_t NbBranchingKmer=0;
    char kmer_seq[sizeKmer+1];
    FILE * file_assembly = fopen(return_file_name(assembly_file),"w+");

    BinaryBank *SolidKmers = new BinaryBank(return_file_name(solid_kmers_file),sizeof(kmer_type),0);

    STARTWALL(assembly);

    char *assemble_only_one_region = NULL; // debugging, set to a ASCII kmer to activate, NULL to desactivate
    bool LOAD_BRANCHING_KMERS=false; // debugging
    bool DUMP_BRANCHING_KMERS=false;
   
    BranchingTerminator *terminator;

    if (LOAD_BRANCHING_KMERS)
    {
        BinaryBank *BranchingKmers = new BinaryBank(return_file_name(branching_kmers_file),sizeof(kmer_type),false);
        terminator = new BranchingTerminator(BranchingKmers,SolidKmers, bloo1,false_positives);
        BranchingKmers->close();
    }
    else
        terminator = new BranchingTerminator(SolidKmers,nb_reads, bloo1,false_positives);

    if (DUMP_BRANCHING_KMERS)
    {
        BinaryBank *BranchingKmers = new BinaryBank(return_file_name(branching_kmers_file),sizeof(kmer_type),true);
        terminator->dump_branching_kmers(BranchingKmers);
        BranchingKmers->close();
    }

#ifdef UNITIG
    SimplePathsTraversal *traversal = new SimplePathsTraversal(bloo1,false_positives,terminator);
    fprintf (stderr,"_________________Assembling in Unitig mode ..._____________________ \n\n");
#else
    MonumentTraversal *traversal = new MonumentTraversal(bloo1,false_positives,terminator);
#endif
    //RandomBranchingTraversal *traversal = new RandomBranchingTraversal(bloo1,false_positives,terminator);
    traversal->set_maxlen(maxlen);
    traversal->set_max_depth(500);
    traversal->set_max_breadth(20);
    
    while (terminator->next(&kmer))
    {
        // keep looping while a starting kmer is available from this kmer
        // everything will be marked during the traversal()'s
        kmer_type starting_kmer;
#ifdef UNITIG
        while (traversal->get_new_starting_node_improved(kmer,starting_kmer))
#else
        while (traversal->find_starting_kmer(kmer,starting_kmer))
#endif
        {
            code2seq(starting_kmer,kmer_seq); // convert starting kmer to nucleotide seq
            traversal->revert_stats(); // set stats from the last commit (discard stats from find_starting_kmer / small contigs)

            if (assemble_only_one_region != NULL)
            {
                kmer_type dummy;
                starting_kmer = extractKmerFromRead(assemble_only_one_region,0,&kmer,&dummy,false);
            }

            // right extension
            len_right = traversal->traverse(starting_kmer,right_traversal,0);
            mlenright= max(len_right,mlenright);

            // left extension, is equivalent to right extension of the revcomp
            len_left = traversal->traverse(starting_kmer,left_traversal,1);
            mlenleft= max(len_left,mlenleft);

            // form the contig
            revcomp_sequence(left_traversal,len_left);
            strcpy(contig,left_traversal); // contig = revcomp(left_traversal)
            strcat(contig,kmer_seq);//               + starting_kmer
            strcat(contig,right_traversal);//           + right_traversal

            contig_len=len_left+len_right+sizeKmer;

            // save the contig
            if(contig_len >= MIN_CONTIG_SIZE)
            {
                max_contig_len = max(max_contig_len,contig_len);
                fprintf(file_assembly,">%lli__len__%lli \n",nbContig,contig_len);
                fprintf(file_assembly,"%s\n",contig);
                nbContig++;
                totalnt+=contig_len;
                traversal->commit_stats();
            }
            else
            {
                traversal->revert_stats();
                nbSmallContig++;
            }
            if (assemble_only_one_region != NULL)
                break;
        }
    
        NbBranchingKmer++;
        if ((NbBranchingKmer%300)==0) fprintf (stderr,"%cLooping through branching kmer n° %lld / %lld  total nt   %lld   ",13,(long long int) NbBranchingKmer,(long long int) terminator->nb_branching_kmers, (long long int)totalnt );

        if (nbContig > 0 && assemble_only_one_region != NULL)
            break;

    }
    fclose(file_assembly);

    fprintf (stderr,"\n Total nt assembled  %lli  nbContig %lli\n",totalnt,nbContig);
    fprintf (stderr," Max contig len  %lli (debug: max len left %lli, max len right %lli)\n",max_contig_len,mlenleft,mlenright);
    fprintf (stderr,"\n Debug traversal stats: %ld ends of contigs (%lld unsaved small contigs), among them:\n",traversal->final_stats.ended_traversals,nbSmallContig);
    fprintf (stderr," %ld couldn't validate consensuses\n",traversal->final_stats.couldnt_validate_consensuses);
    fprintf (stderr," %ld large bubble breadth, %ld large bubble depth, %ld marked kmer, %ld no extension\n",traversal->final_stats.couldnt_traverse_bubble_breadth,traversal->final_stats.couldnt_traverse_bubble_depth,traversal->final_stats.couldnt_because_marked_kmer,traversal->final_stats.couldnt_find_extension);
    fprintf (stderr," %ld in-branchin large depth, %ld in-branching large breadth, %ld in-branching other\n",traversal->final_stats.couldnt_inbranching_depth,traversal->final_stats.couldnt_inbranching_breadth,traversal->final_stats.couldnt_inbranching_other);
    
    STOPWALL(assembly,"Assembly");

    free(left_traversal);
    free(right_traversal);
    free(contig);
    SolidKmers->close();
}