#include "Debloom.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

// GUS: the false positive set can be either FPSet or FPSetCascading4,
// both inheriting from Set. The choice is made by the function called
// to load the false positives: load_false_positives() or
// load_false_positives_cascading4().
Set *false_positives; 


FILE * F_debloom_read;
FILE * F_debloom_write;
uint64_t n_false_positives=0;

Hash16 * hasht1;

void end_debloom_partition(bool last_partition)
{

    int value;
    char false_positive_kmer_char[sizeKmer+1];
    FILE *file_false_positive_kmers =NULL;
    kmer_type graine;

    /////////////////////////begin write files 
    rewind (F_debloom_read);
    rewind (F_debloom_write);

	#ifndef MINGW
	ftruncate(fileno(F_debloom_write), 0); //erase previous file 
	#else // tempfix? fileno is not accepted by mingw
	fclose(F_debloom_write);
	F_debloom_write = fopen(return_file_name("debloom2"),"wb+");
	#endif

    if (last_partition)
    {   
        // write false positive kmers to fasta file
        file_false_positive_kmers = fopen(return_file_name(false_positive_kmers_file),"wb");
    }

    n_false_positives = 0;
    while(fread(&graine, sizeof(graine),1, F_debloom_read)){

        if(hasht1->get(graine,&value)==0) //kmer not present == kmer not solid
        {
            n_false_positives ++;

            if (!fwrite(&graine, sizeof(graine), 1, F_debloom_write))
            {
                printf("error: can't fwrite (disk full?)\n");
                exit(1);
            }


            if (last_partition)
            {
                code2seq(graine,false_positive_kmer_char);
                fprintf(file_false_positive_kmers,">fp\n");
                fputs(false_positive_kmer_char,file_false_positive_kmers);
                fprintf(file_false_positive_kmers,"\n");
            }
        }
        //else kmer is a true positive, do nothing

    }

    if (last_partition)   
        fclose(file_false_positive_kmers);
} 

int debloom(int order, int max_memory, Bloom * bloo1)
{
    STARTWALL(pos);

    FILE * debloom_file = fopen(return_file_name("debloom"),"wb+");

    ifstream solidReads;
    solidReads.open(solid_reads_file);

    kmer_type new_graine, next_real, kmer;
    string read;
    int nt, nextNuc, nucPos, strand;
   
    int NbCandKmer=0, NbRawCandKmer = 0;
    uint64_t NbSolidKmer =0;

    // write all positive extensions in disk file
    printf("Weight before debloom: %d \n", bloo1->weight());
    while (getline(solidReads, read))
    {
 
      //printf("Read: %s \n", &read[0]);
      getFirstKmerFromRead(&kmer,&read[0]);

      for (int i = 0; i <= read.length() - sizeKmer ; i++, 
        shift_kmer(&kmer, NT2int(read[i+sizeKmer-1]), 0)){
        
          //printf("This kmer: %s \n", print_kmer(kmer));
         // printf("Extensions: \n");
        
        for (strand = 0; strand < 2 ; strand++){
              //printf("strand %d :\n", strand);

              nucPos = i+sizeKmer;
              if(strand == 1){
                  nucPos = i-1;
              }

              next_real = next_kmer(kmer,NT2int(read[nucPos]), strand);
            
              for(nt=0; nt<4; nt++) {
                new_graine = next_kmer(kmer,nt, strand);
                //printf("%s ", print_kmer(new_graine));
                if(bloo1->contains(get_canon(new_graine))){ 
                    //printf("is in the filter ");
                    NbRawCandKmer++;
                    if(new_graine != next_real){
                      //printf("and is not in the read. ");
                                          //printf("%s is a positive extension! \n", print_kmer(new_graine));
                  
                      NbCandKmer++;
                      if (!fwrite(&new_graine, sizeof(new_graine), 1, debloom_file))
                      {
                        printf("error: can't fwrite (disk full?)\n");
                        exit(1);
                      }
                    }
                }
                //printf("\n");
            }
        }
        NbSolidKmer++;
        if ((NbSolidKmer%10000)==0) fprintf (stderr,"%c Writing positive Bloom Kmers %lld",13,NbSolidKmer);
      }
      }
      solidReads.close();
      printf ("\n Number of candidate kmers %d \n",NbCandKmer);

      printf (" Number of raw candidate kmers %d \n",NbRawCandKmer);
      printf ("Estimated false positive rate: %f \n", float(NbCandKmer)/float(NbSolidKmer*6));
    /* nbkmers_solid =  NbSolidKmer; // GUS: it's global now

    fprintf(stderr,"\n%lli kmers written\n",cc);

    STOPWALL(pos,"Write all positive kmers");

    STARTWALL(deb);

    double bl1tai =  (double)bloo1->tai ;
    delete bloo1;

    // now that bloo1 is deleted, initialize hasht1
    int NBITS_HT = max( (int)ceilf(log2f((0.1*max_memory*1024L*1024L)/sizeof(cell_ptr_t))), 1); // set hasht1 cells to occupy 0.1 * [as much mem as poss]
    hasht1 =new Hash16(NBITS_HT); 
    
    ////////////////////////////////////////////////////////////////   --find false positive, with hash table partitioning
    uint64_t max_kmer_per_part = (uint64_t) (0.8*max_memory*1024LL*1024LL /sizeof(cell<kmer_type>));
    //adapter taille ht en fonction
    

    printf("%d partitions will be needed\n",(int)(nbkmers_solid/max_kmer_per_part));

    NbSolidKmer =0;
    int numpart = 0;
    
    return 1;

    SolidKmers->rewind_all();

    // deblooming:
    // read the list of (non-redundant) solid kmers and load it, in chunks, into a hash table
    // at each pass, check all the positive extensions and keep those which are not indicated, by the current chunk, as solid kmers
    // at the end, only the positive extensions which are not solid are kept
    while (SolidKmers->read_element(&kmer))
    {
        hasht1->add(kmer);

        NbSolidKmer++;
        if ((NbSolidKmer%10000)==0) fprintf (stderr,"%cBuild Hash table %lld",13,NbSolidKmer);

        if(hasht1->nb_elem >max_kmer_per_part) //end partition,  find false positives
        {
            fprintf(stderr,"End of debloom partition  %lli / %lld \n",hasht1->nb_elem,max_kmer_per_part);

            end_debloom_partition(false);

            //swap file pointers
            F_tmp = F_debloom_read;
            F_debloom_read = F_debloom_write;
            F_debloom_write = F_tmp;
            /////////end write files

            //reset hash table
            hasht1->empty_all();

            fprintf(stderr,"\n%lli false positives written , partition %i \n",n_false_positives,numpart);

            numpart++;
        } ///end partition


    }
    //fprintf(stderr,"Nb kmers stored in the bloom table %lld\n",nbkmers_solid);


    ///////////////////////// last partition, will write all the FP's to the good file

    end_debloom_partition(true); 

    /////////end write files


    fprintf(stderr,"Total nb false positives stored in the Debloom hashtable %lli \n",n_false_positives);

    delete hasht1;


    STOPWALL(deb,"Debloom");
 
    // GUS: will use to output summary later
    b1_size = (uint64_t) bl1tai;
  
    fclose(debloom_file);
    fclose(debloom_file_2);
    SolidKmers->close();


    return 1;
    */
}

uint64_t countFP(Bank *FalsePositives)
{
  char * rseq;
  int readlen;
  uint64_t nbFP = 0;

  while (FalsePositives->get_next_seq(&rseq,&readlen))
    nbFP++;
  
  FalsePositives->rewind_all();
  return nbFP;
}

Set *load_false_positives() 
{
    int64_t NbInsertedKmers = 0;
    char * rseq;
    int readlen;
    kmer_type kmer, graine, graine_revcomp;

    Bank *FalsePositives = new Bank(return_file_name(false_positive_kmers_file));

    // alloc false positives with the just the right estimated size

    uint64_t nbFP = countFP(FalsePositives);

    FPSet *fp = new FPSet(nbFP);
    
    while (FalsePositives->get_next_seq(&rseq,&readlen))
    {
        kmer = extractKmerFromRead(rseq,0,&graine,&graine_revcomp);
                
        fp->insert(kmer);

        NbInsertedKmers++;

        if ((NbInsertedKmers%10000)==0) fprintf (stderr,(char*)"%cInsert false positive Kmers in hash table %lld",13,NbInsertedKmers);
    }
    fp->finalize(); // always call this when finishing to create a FPSet

    fprintf (stderr,"\nInserted %lld false positive kmers in the hash structure.\n\n",NbInsertedKmers);

//    print_size_summary(fp);

    return fp;
}


Set *dummy_false_positives() 
{
    FPSet *fp = new FPSet((uint64_t)1);
    return fp;
}

Set *load_false_positives_cascading4()
{
  int64_t NbInsertedKmers;
  char * rseq;
  int readlen;
  kmer_type kmer, graine, graine_revcomp;

  
  // **** Initialize B2, B3, B4 and T4 ****
  Bank *FalsePositives = new Bank(return_file_name(false_positive_kmers_file));
  uint64_t nbFP = countFP(FalsePositives);
  
  FPSetCascading4 *fp = new FPSetCascading4;
  
  fp->bloom2 = new Bloom((uint64_t)(nbFP * NBITS_PER_KMER));
  fp->bloom2->set_number_of_hash_func((int)floorf(0.7*NBITS_PER_KMER));

  uint64_t estimated_T2_size = max((int)ceilf(nbkmers_solid * (double)powf((double)0.62, (double)NBITS_PER_KMER)), 1);
  uint64_t estimated_T3_size = max((int)ceilf(nbFP          * (double)powf((double)0.62, (double)NBITS_PER_KMER)) ,1);

  fp->bloom3 = new Bloom((uint64_t)(estimated_T2_size * NBITS_PER_KMER));
  fp->bloom3->set_number_of_hash_func((int)floorf(0.7*NBITS_PER_KMER));

  fp->bloom4 = new Bloom((uint64_t)(estimated_T3_size * NBITS_PER_KMER));
  fp->bloom4->set_number_of_hash_func((int)floorf(0.7*NBITS_PER_KMER));


  // **** Insert the false positives in B2 ****
  NbInsertedKmers = 0;
  while (FalsePositives->get_next_seq(&rseq,&readlen))
  {
    kmer = extractKmerFromRead(rseq,0,&graine,&graine_revcomp);
    
    fp->bloom2->add(kmer);
    
    NbInsertedKmers++;
   // if ((NbInsertedKmers%10000)==0)
     // fprintf (stderr,"%cInsert false positive B2 %lld",13,NbInsertedKmers);
  }
  //fprintf (stderr,"%cInsert false positive B2 %lld", 13,NbInsertedKmers);
  FalsePositives->close();

  DEBUGE(("\nInserted %lld (estimated, %lld) kmers in B2.\n", NbInsertedKmers, nbFP));


  //  **** Insert false positives in B3 and write T2 
  int addKmers = 0;
  NbInsertedKmers = 0;
  FILE *T2_file = fopen(return_file_name("t2_kmers"), "w+"); // We will read this file later, when filling T4 
  BinaryBank *SolidKmers = new BinaryBank(return_file_name(solid_kmers_file),sizeof(kmer),0);
  while(SolidKmers->read_element(&kmer))
  {
    if (fp->bloom2->contains(kmer))
    {
      if (!fwrite(&kmer, sizeof(kmer), 1, T2_file))
      {
	printf("error: can't fwrite (disk full?)\n");
	exit(1);
      }

      fp->bloom3->add(kmer);
      addKmers++;
    }

    NbInsertedKmers++;
    //if ((NbInsertedKmers%10000)==0)
      //fprintf (stderr,(char*)"%cInsert false positive B3 %lld",13,NbInsertedKmers);
  }
 // fprintf (stderr,(char*)"%cInsert false positive B3 %lld",13,NbInsertedKmers);
  SolidKmers->close();

  DEBUGE(("\nInserted %d (estimated, %llu) kmers in B3.\n", addKmers, estimated_T2_size));

  
  // **** Insert false positives in B4 (we could write T3, but it's not necessary)
  FalsePositives = new Bank(return_file_name(false_positive_kmers_file));
  NbInsertedKmers = 0;
  addKmers = 0;
  while (FalsePositives->get_next_seq(&rseq,&readlen))
  {
    kmer = extractKmerFromRead(rseq,0,&graine,&graine_revcomp);
    
    if (fp->bloom3->contains(kmer))
    {
      fp->bloom4->add(kmer);
      addKmers++;
    }

    NbInsertedKmers++;
    //if ((NbInsertedKmers%10000)==0)
      //fprintf (stderr,"%cInsert false positive B4 %lld",13,NbInsertedKmers);
  }
  //fprintf (stderr,"%cInsert false positive B4 %lld", 13,NbInsertedKmers);
  FalsePositives->close();

  DEBUGE(("\nInserted %d (estimated, %lld) kmers in B4.\n", addKmers, estimated_T3_size));
  

  // **** Count and insert false positives in T4
  rewind(T2_file);
  addKmers = 0;
  while (fread(&kmer, sizeof(kmer), 1, T2_file))
    if (fp->bloom4->contains(kmer))
      addKmers++;

  fp->false_positives = new FPSet(addKmers);
  rewind(T2_file);
  addKmers = 0;
  NbInsertedKmers = 0;
  while (fread(&kmer, sizeof(kmer), 1, T2_file))
  {  
    if (fp->bloom4->contains(kmer))
    {
      fp->false_positives->insert(kmer);
      addKmers++;
    }

    NbInsertedKmers++;
   // if ((NbInsertedKmers%10000)==0)
    //  fprintf (stderr,"%cInsert false positive T4 %lld",13,NbInsertedKmers);
  }
  fp->false_positives->finalize();
 // fprintf (stderr,"%cInsert false positive T4 %lld", 13,NbInsertedKmers);
  fclose(T2_file);

  DEBUGE(("\nInserted %d (estimated, %lld) kmers in T4.\n\n", addKmers, (uint64_t)fp->false_positives->capacity()));
  
//  print_size_summary(fp);

  return fp;
}

double toMB(double value)
{
  return value / 8LL/1024LL/1024LL;
}

void print_size_summary(FPSet *fp)
{
  int bits_per_FP_element = FPSet::bits_per_element;

  uint64_t size_B1 = b1_size,
           size_T1 = fp->capacity() * FPSet::bits_per_element;
  double total_size = (double)(size_B1 + size_T1);

  fprintf(stderr,"Size of the Bloom table  : %.2lf MB\n", toMB(size_B1) );  
  fprintf(stderr,"                           %.2lf bits / solid kmer\n", b1_size/(double)(nbkmers_solid) );  
  fprintf(stderr, "Size of the FP table     : %lli FP x %d bits =  %.2lf MB  \n", fp->capacity(), bits_per_FP_element, toMB((double)(size_T1)) );
  fprintf(stderr,"                                      actual implementation : %.2lf bits / solid kmer\n", size_T1/(double)nbkmers_solid);
  fprintf(stderr,"  assuming list of kmers, i.e. sizeof(kmer_type) bits / FP : %.2lf bits / solid kmer \n\n",(fp->capacity()*sizeof(kmer_type)*8LL)/(double)(nbkmers_solid));
  fprintf(stderr,"      Total %.2lf MB for %lld solid kmers  ==>  %.2lf bits / solid kmer\n\n", toMB(total_size), nbkmers_solid, total_size / nbkmers_solid);  
}

void print_size_summary(FPSetCascading4 *fp)
{
  uint64_t size_B1 = b1_size,
    size_B2 = fp->bloom2->tai,
    size_B3 = fp->bloom3->tai,
    size_B4 = fp->bloom4->tai,
    size_T4 = fp->false_positives->capacity() * FPSet::bits_per_element;
  double total_size = (double)(size_B1 + size_B2 + size_B3 + size_B4 + size_T4);

  DEBUGE((stderr,"Size of the Bloom table (B1)  : %.2lf MB\n", toMB((double)size_B1)));
  DEBUGE((stderr,"Size of the Bloom table (B2)  : %.2lf MB\n", toMB((double)size_B2)));
  DEBUGE((stderr,"Size of the Bloom table (B3)  : %.2lf MB\n", toMB((double)size_B3)));
  DEBUGE((stderr,"Size of the Bloom table (B4)  : %.2lf MB\n", toMB((double)size_B4)));
  DEBUGE((stderr,"Size of the FP table (T4)     : %.2lf MB\n", toMB((double)size_T4)));
  fprintf(stderr,"      Total %.2lf MB for %lld solid kmers  ==>  %.2lf bits / solid kmer\n\n", toMB(total_size), nbkmers_solid, total_size / nbkmers_solid);
}


//j-checking code from legacy.  Moved down here since we may want to refer to it but don't need it yet.
/*

  // maybe do more lax deblooming; if it's a dead-end, it's no big sdeal, don't pass it to the false positive test
                    // what would have been needed if i decided to enable order>0 (but actually this won't happen): 
                    //  - better estimate of structure size in the presence of order>0 deblooming  
                    

if (order == 1)  // this case just detects tips
                    {
                        bool is_linked = false;
                        for(int tip_nt=0; tip_nt<4; tip_nt++) 
                        {
                            int new_strand = strand;
                            kmer_type kmer_after_possible_tip = next_kmer(new_graine,tip_nt, &new_strand);
                            if(bloo1->contains(kmer_after_possible_tip))
                            {
                                is_linked = true;
                                break;
                            }
                        }
                        if (!is_linked)
                            continue; // it's a tip, because it's linked to nothing
                    }
    
                    if (order > 1) // general case. should work for order = 1, but i coded an optimized version above
                    { 
                        Frontline frontline( new_graine, strand, bloo1, NULL, NULL, NULL);
                        while (frontline.depth < order)
                        {
                            frontline.go_next_depth();
                            if (frontline.size() == 0)
                                break;
                            // don't allow a breadth too large anywqy
                            if (frontline.size()> 10)
                                break;
                        }
                        if (frontline.size() == 0)
                            continue; // it's a deadend
                    }
*/