#!/bin/bash
#Runs Mink with k = 27, max read length = 100, estimated kmer size = 60 million, two hash functions
#Takes on parameter- the output file prefix
#Streams the input file from Roye's website at http://www.tau.ac.il/~rozovr/chr20.c50.fa.gz

URL_FILE=wget_urls
READ_COMMAND=wget\ --read-timeout=5\ --timeout=15\ -t\ 0\ -qO-\ -i\ $URL_FILE\ \|\ bzip2\ -d\ -c\ -q

eval "./mink -read_load_file <($READ_COMMAND) -read_scan_file <($READ_COMMAND) -size_kmer 31 -max_read_length 130 -estimated_kmers 3000000000 -file_prefix $1 --two_hash --paired_ends --fastq"

# eval "./mink -size_kmer 27 \
# -max_read_length 130 \
# -estimated_kmers 3000000000 \
# -read_load_file <($READ_COMMAND) \
# -file_prefix $1 \
# --two_hash \
# --just_load_bloom \
# --fastq "

# eval "./mink -size_kmer 27 \
# -max_read_length 130 \
# -estimated_kmers 3000000000 \
# -read_scan_file <($READ_COMMAND) \
# -bloom_file $1.bloom \
# -file_prefix $1 \
# --two_hash \
# --fastq \
# --no_cleaning"