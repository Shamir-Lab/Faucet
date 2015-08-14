#!/bin/bash
#Runs Mink with k = 27, max read length = 100, estimated kmer size = 60 million, two hash functions
#Takes on parameter- the output file prefix
#Streams the input file from Roye's website at http://www.tau.ac.il/~rozovr/chr20.c50.fa.gz

#Command: wget -qO- ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/SRA010/SRA010896/SRX016231* | bzip2 -d -c -q
ADDRESS=ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/SRA010/SRA010896/SRX016231/
READ_COMMAND=wget\ -qO-\ $ADDRESS*\ \|\ bzip2\ -d\ -c\ -q

eval "./mink -size_kmer 27 \
-max_read_length 130 \
-estimated_kmers 3000000000 \
-read_load_file <($READ_COMMAND) \
-file_prefix $1 \
--two_hash \
--just_load_bloom \
--fastq "

eval "./mink -size_kmer 27 \
-max_read_length 130 \
-estimated_kmers 3000000000 \
-read_scan_file <($READ_COMMAND) \
-bloom_file $1.bloom \
-file_prefix $1 \
--two_hash \
--fastq "