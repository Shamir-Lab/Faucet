# import minia bf to nim
# need Bloom.h
# "LargeInt.h"
# "ttmath/ttmath.h"

{.link: "./nimcache/Bloom.o".}

# type
#     Bloom{.importcpp: "Bloom::Bloom".} = object
    
# proc cnew*[T](x: T): ptr T {.importcpp: "(new '*0#@)", nodecl.}
# # template bloom*(a: int): Bloom =
#   # constructBloom(a.int)
# proc constructBloom(a : cint): Bloom {.importcpp: "Bloom".}
# # proc add(var a: Bloom, b: string){.importcpp: "#.add(@)".} {.urh.}
# # proc contains(var a:  Bloom, b: string): bool {.importcpp: "#.contains(@)".} {.urh.}

# # proc saveXML*(sc: ptr Scene, f: UrFile) {.importcpp: "#.SaveXML(@)", header:
# #   "File.h".}

# let BF  =  cnew constructBloom(10_000)
# add(BF, "yes")
# echo(contains(BF "yes"))
# echo(contains(BF, "no"))


# create new BF as done in minia code
#### functions to use
# Bloom * bloo1;
# bloo1 = new Bloom((uint64_t)estimated_BL1_freesize);
# bloo1->set_number_of_hash_func((int)floorf(0.7*NBITS_PER_KMER));   
# bloo1->add(whatever)
# bloo1->contains(new_graine)
# bloo1->dump(filename)
# bloo1->load(filename)



# run tests used for nim bloom package on minia

######### how it's used in minia code ###########
#### functions found in Minia's Utils.cpp ####
# - wrapper for setting up first BF and loading its k-mers
# bloo1 = bloom_create_bloo1((BloomCpt *)NULL, false);

# bloom_pass_reads_binary(bloo1, bloom_counter, (char*)"%cInsert solid Kmers in Bloom %lld"); # called with NULL as bloom counter, solid_kmers as false (from main)
# BinaryBank * SolidKmers = new BinaryBank(return_file_name(solid_kmers_file),sizeof(kmer),0);

#   while(SolidKmers->read_element(&kmer))
#     {
#       // printf("kmer %lld\n",kmer);
#       bloom_to_insert->add(kmer);
#       NbInsertedKmers++;
#       NbRead++;
#       if ((NbRead%10000)==0) fprintf (stderr,stderr_message,13,(long long)NbRead);
#     }
# #  fprintf (stderr,"\nInserted %lld %s kmers in the bloom structure.\n",(long long)NbInsertedKmers,"solid");
#   SolidKmers->close();  