#!/bin/bash
#Runs Mink with k = 27, max read length = 100, estimated kmer size = 60 million, two hash functions
#Takes on parameter- the output file prefix
#Streams the input file from Roye's website at http://www.tau.ac.il/~rozovr/chr20.c50.fa.gz

READ_COMMAND=wget\ -qO-\ http://www.tau.ac.il/~rozovr/chr20.c50.fa.gz\ \|\ zcat\ -f

eval "./mink -size_kmer 27 \
-max_read_length 100 \
-estimated_kmers 60000000 \
-read_load_file <($READ_COMMAND) \
-file_prefix $1 \
--two_hash \
--just_load_bloom"

eval "./mink -size_kmer 27 \
-max_read_length 100 \
-estimated_kmers 60000000 \
-read_scan_file <($READ_COMMAND) \
-bloom_file $1.bloom \
-file_prefix $1 \
--two_hash"