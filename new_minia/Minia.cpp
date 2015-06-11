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
int j = 0;

char* junctions_filename = new char[100];
char* solids_file = new char[100];
bool from_kmers = true;
int read_length;
uint64_t genome_size;
bool get_kmers = false;//only to initially get the kmer files so i could use them, don't want to really use this code yet
int64_t nb_reads;
char* kmer_filename = (char*)"solid_27mers_100k";
#include "Bloom.h"
#include "ReadScanner.h"
#include "Kmer.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

set<kmer_type> all_kmers;
/*
To run the new version, type make in the directory to compile.

Then, type ./minia 1 2 3 4 5 6 7 [8], where
1 = name of reads file (current format is each line has a string of characters representing the read)
2 = k
3 = read length
4 = number of distinct kmers.  This will be directly used to size the bloom filter so try to have a good estimate.
5 = false positive rate (.01 is a good base rate) 
6 = j.  j = 0 corresponds to taking direct extensions of the reads, j = 1 is extensions of extensions, etc.
7 = from_kmers.  Put in 0, it will be from_kmers = false, do it like normal from reads.  Put in 1, it does from kmers.
8 = junctions file name.  Will print out to this file. If no filename is given, no file will be printed.

[] indicate optional arguments.

This will load a bloom filter with all the kmers from the reads, then scan through them inserting potential false positives.
No error correction corrently, and no assembly.
*/

inline int handle_arguments(int argc, char *argv[]){
if(argc <  8)
    {
        fprintf (stderr,"usage:\n");
        fprintf (stderr," %s input_file kmer_size min_abundance estimated_nb_reads prefix\n",argv[0]);
        fprintf (stderr,"hints:\n min_abundance ~ 3\n estimated_nb_reads is in bp, does not need to be accurate, only controls memory usage\n prefix is any name you want the results to start with\n");

        return 1;
    }

    //1st arg: read file name
    strcpy(solids_file,argv[1]);
    printf("Reads file name: %s \n", solids_file);

    // 2rd arg: kmer size.
    sizeKmer=27;
    sizeKmer = atoi(argv[2]);
    if (sizeKmer>((int)sizeof(kmer_type)*4))
    {
        printf("Max kmer size on this compiled version is %u\n",sizeof(kmer_type)*4);
        exit(1);
    }
    if (sizeKmer == (int)(sizeof(kmer_type)*4))
        kmerMask = -1;
    else
        kmerMask=(((kmer_type)1)<<(sizeKmer*2))-1;
    printf("k: %d: \n", sizeKmer);

    //3rd arg: read length
    read_length = atoi(argv[3]);

    //4th arg: number of reads
    genome_size = atoll(argv[4]);
    printf("Genome size: %lli .\n", genome_size);
    
    //5th arg: false posiive rate
    fpRate = atof(argv[5]);

    //6th arg: j
    j = atoi(argv[6]);
    printf("j: %d \n", j);

    //7th arg: whether to do from kmers or from reads
    from_kmers = atoi(argv[7]);
    if(from_kmers)
        printf("Working from kmers. \n");
    else
        printf("Working from reads. \n");

    if(argc > 8){
        strcpy(junctions_filename,argv[8]);
        printf("Junctions file name: %s \n", junctions_filename);
    }
}

inline void load_filter_from_reads(Bloom* bloo1, const char* reads_filename){
    ifstream solidReads;
    solidReads.open(reads_filename);

    int readsProcessed = 0;
    kmer_type new_graine, kmer;
    string read;
    time_t start;
    time_t stop;
    time(&start);
    printf("Weight before load: %ld \n", bloo1->weight());
    while (getline(solidReads, read))
    {
        getFirstKmerFromRead(&kmer,&read[0]);

        for (int i = 0; i <= read.length() - sizeKmer ; i++, 
          shift_kmer(&kmer, NT2int(read[i+sizeKmer-1]), 0)){
            bloo1->add(get_canon(kmer));
            readsProcessed++;
            if ((readsProcessed%10000)==0) fprintf (stderr,"%c %lld",13,(long long)readsProcessed);
        }
    }

    solidReads.close();
    printf("\n");
    printf("Weight after load: %ld \n", bloo1->weight()); 
    time(&stop);
    printf("Time to load: %f \n", difftime(stop,start));
}

inline void load_filter_from_kmers(Bloom* bloo1, const char* kmers_filename){
    ifstream solidKmers;
    solidKmers.open(kmers_filename);

    int kmersProcessed = 0;
    kmer_type left,right;
    string kpomer;

    //printf("Weight before load: %ld \n", bloo1->weight());
    while (getline(solidKmers, kpomer))
    {
        //printf("kpomer %s \n", &kpomer[0]);
        getFirstKmerFromRead(&left,&kpomer[0]);
        right = next_kmer(left, NT2int(kpomer[sizeKmer]),0);
        //printf("left %s \n", print_kmer(left));
        //printf("right %s \n", print_kmer(right));
        bloo1->add(get_canon(left));
        bloo1->add(get_canon(right));
        kmersProcessed++;
        if ((kmersProcessed%10000)==0) fprintf (stderr,"%c %lld",13,(long long)kmersProcessed);
    }

    solidKmers.close();
    printf("\n");
    printf("Weight after load: %ld \n", bloo1->weight()); 
}

inline Bloom* create_bloom_filter(int estimated_items, float fpRate){
     
    Bloom * bloo1;
    int bits_per_item = -log(fpRate)/log(2)/log(2); // needed to process argv[5]

    // int estimated_bloom_size = max( (int)ceilf(log2f(nb_reads * NBITS_PER_KMER )), 1);
    uint64_t estimated_bloom_size = (uint64_t) (genome_size*bits_per_item);
    printf("Estimated items: %lli \n", genome_size);
    printf("Estimated bloom size: %d .\n", (int)estimated_bloom_size);
    
    printf("BF memory: %f MB\n", (float)(estimated_bloom_size/8LL /1024LL)/1024);
    bloo1 = new Bloom(estimated_bloom_size);

    printf("Bits per kmer: %d \n", bits_per_item);

    printf("Number of hash functions: %d \n", (int)floorf(0.7*bits_per_item));
    bloo1->set_number_of_hash_func((int)floorf(0.7*bits_per_item));

return bloo1;
}

void write_kmers(const char* reads_filename){
    ifstream solidReads;
    solidReads.open(reads_filename);
    kmer_type kmer;
    char* kmerSeq = new char[sizeKmer];
    int readsProcessed = 0;
    printf("Building kmer set.\n");
    string read;
    while (getline(solidReads, read))
    {
        getFirstKmerFromRead(&kmer,&read[0]);

        for (int i = 0; i <= read.length() - sizeKmer ; i++, 
            shift_kmer(&kmer, NT2int(read[i+sizeKmer-1]), 0)){
            all_kmers.insert(kmer);
            readsProcessed++;
            if ((readsProcessed%10000)==0) fprintf (stderr,"%c %lld",13,(long long)readsProcessed);
        }
    }
    solidReads.close();
    printf("Done building kmer set.\n");
    ofstream solidKmers;
    solidKmers.open(kmer_filename);
    
    printf("Writing to kmer file\n");
    set<kmer_type>::iterator it;
    for(it = all_kmers.begin(); it != all_kmers.end(); it++){
        code2seq(*it, kmerSeq);
        solidKmers << kmerSeq;
        solidKmers << '\n';
    }
    printf("Done writing to kmer file.\n");
    solidKmers.close();

}

int main(int argc, char *argv[])
{
    if (handle_arguments(argc, argv) == 1){
        return 1;
    }

    if(get_kmers){
        write_kmers(solids_file);
    }

    int estimated_kmers = genome_size;
    Bloom* bloo1 = create_bloom_filter(estimated_kmers, fpRate);
   
   // if(from_kmers){
   //      load_filter_from_kmers(bloo1, solids_file);
   //      debloom_kpomerscan(solids_file, bloo1, j);
   //  }
    load_filter_from_reads(bloo1, solids_file);
    ReadScanner* scanner = new ReadScanner(solids_file, bloo1);
    scanner->setJ(j);

    scanner->scanReads(genome_size);
    scanner->printScanSummary();
    if(argc > 8){
        scanner->junctionMapToFile(junctions_filename);
    }   
    printf("Program reached end. \n");
    return 0;
}



/*

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
*/