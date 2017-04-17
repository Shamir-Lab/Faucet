# Getting Faucet
    git clone https://github.com/rozovr/Faucet.git
    cd Faucet/src
	make depend
	make    

# Running Faucet (locally)
Example usage:

	./faucet -read_load_file interlaced_reads.fq -read_scan_file interlaced_reads.fq -size_kmer 31 -max_read_length 100 -estimated_kmers 1000000000 -singletons 200000000 -file_prefix faucet_outputs --fastq --paired_ends

The above command takes as input the file interlaced_reads.fq (where entries alternate between mates 1 and 2 of a paired end library), and the input format is fastq. Faucet does not accept separate mate files, but can accept fasta format and files composed of read sequences alone.

# Streaming from a remote source
A demonstration streaming reads from a remote server is provided in the script src/stream_data_from_urls_list.sh


# Requirements
Faucet was implemented in C++ 11, so requires a compiler that is not too ancient to support it, and has been tested only on Linux so far. 

# Detailed usage

Usage:
./faucet -read_load_file <filename> -read_scan_file <filename> -size_kmer <k> -max_read_length <length> -estimated_kmers <num_kmers> -singletons <num_kmers> -file_prefix <prefix>
Optional arguments: --fastq -max_spacer_dist <dist> -fp rate <rate> -j <int> --two_hash -bloom_file <filename> -junctions_file <filename> --paired_ends --no_cleaning

### required arguments:
 
	-read_load_file <filename>, a file name string 
	-read_scan_file <filename> , a file name string
	-size_kmer <k> , and odd integer <= 31
	-max_read_length <length>, the longest read length in the data (e.g., if the reads were trimmed to different sizes)
	-estimated_kmers <num_kmers> 
	-singletons <num_kmers> 
	-file_prefix <prefix>, the desired prefix string or directory path for output files 
 
we recommend applying <a href="https://github.com/bcgsc/ntCard">ntCard</a> to extract the number estimated k-mers (F0) and singletons (f1) in the dataset.


License
=======


* Low level code for dealing with binary encoded k-mers and strings, and for Bloom filters is derived from the original minia implementation, http://minia.genouest.org/; these components, mostly unmodified, are distributed under a GPL 3.0 license

* Original code is distributed under the BSD 3 clause license.
