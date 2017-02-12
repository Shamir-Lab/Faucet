#!/bin/bash
#Runs Faucet
#Takes 4 parameter- 
#1) the output file prefix and a 
#2) file listing the input file URLs
#3) estimated kmers
#4) singletons

URL_FILE=$2
READ_COMMAND=wget\ --read-timeout=5\ --timeout=15\ -t\ 0\ -qO-\ -i\ $URL_FILE\ \|\ bzip2\ -d\ -c\ -q

eval "./faucet -read_load_file <($READ_COMMAND) -read_scan_file <($READ_COMMAND) -size_kmer 31 -max_read_length 130 -estimated_kmers $3 -singletons $4 -file_prefix $1 --fastq"

# eval "./faucet -size_kmer 27 \
# -max_read_length 130 \
# -estimated_kmers 3000000000 \
# -read_load_file <($READ_COMMAND) \
# -file_prefix $1 \
# --two_hash \
# --just_load_bloom \
# --fastq "

# eval "./faucet -size_kmer 27 \
# -max_read_length 130 \
# -estimated_kmers 3000000000 \
# -read_scan_file <($READ_COMMAND) \
# -bloom_file $1.bloom \
# -file_prefix $1 \
# --two_hash \
# --fastq \
# --no_cleaning"