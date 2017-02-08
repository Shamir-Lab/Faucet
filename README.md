<!-- # Getting Recycler
You can download Recycler [here](https://github.com/Shamir-Lab/Recycler/releases/download/Recycler-v0.6/Recycler-0.6.zip) or clone it via the link below. In case you download the zip, unzip the file before following the instructions below (ignoring the 'git clone' line)
 -->
# Getting Faucet

    git clone https://github.com/rozovr/Faucet.git
    cd Faucet/readScanAssembler
	make depend
	make    
    
# Running Faucet on reads that are on a local hard drive


# Introduction

Faucet is a streaming algorithm for constructing de novo assembly graphs from large metagenomes. 

# Requirements

Faucet was implemented in C++ 11, so requires a compiler that is not too ancient to support it, and has been tested only on Linux so far. 


# Detailed usage

python recycle.py -g GRAPH -k MAX_K -b BAM [-l LENGTH] [-m MAX_CV] [-i ISO]

### required arguments:
 <!--    
    -g GRAPH
    (spades 3.50+) assembly graph FASTG file to process:
    for spades 3.5, before_rr.fastg; for spades 3.6+, assembly_graph.fastg
    -k MAX_K
    integer reflecting maximum k value used by the assembler
    -b BAM
    BAM file resulting from aligning reads to contigs file, filtering for best matches
  -->
### optional arguments:
<!-- 
    -l LENGTH
    minimum length required for reporting [default: 1000]
    -m MAX_CV
    coefficient of variation used for pre-selection
    [default: 0.5, higher - less restrictive]
    -i ISO
    True or False value reflecting whether data sequenced
    was an isolated strain 
 -->
# <a name="bam-prep">Preparing the BAM input:
<!-- 
Recycler uses paired-end alignments of the reads originally assembled to the output assembly graph to filter and select amongst candidate circular sequences. In order to do so, it requires as input a BAM file containing the set of best alignment hits for each read pair. We recommend the following steps (tested on BWA 0.7.5 and samtools 1.19) to prepare the BAM file:

    python make_fasta_from_fastg.py -g assembly_graph.fastg
    
    bwa index assembly_graph.nodes.fasta
    
    bwa mem  assembly_graph.nodes.fasta R1.fastq.gz R2.fastq.gz | samtools view -buS - > reads_pe.bam
    
    samtools view -bF 0x0800 reads_pe.bam > reads_pe_primary.bam
    
    samtools sort reads_pe_primary.bam reads_pe_primary.sort
    
    samtools index reads_pe_primary.sort.bam

following these steps, we only need the files reads_pe_primary.sort.bam and reads_pe_primary.sort.bam.bai.
 -->
# Outputs:

<!-- 1. \<prefix\>.cycs.fasta  - a fasta formatted file of predicted plasmids
2. \<prefix\>.cycs.paths_w_cov.txt - a text file containing information about plasmids composed of multiple contigs.

The format for the second file is:
* *\<plasmid name\>* - e.g., RNODE_5_length_42666_cov_19.93685
* *\<node names in the original graph making up this cycle\>* - e.g., \('NODE_2801_length_42596_cov_19.8677', "NODE_2387_length_125_cov_34.7286'"\).
* *\<coverage levels of nodes at the time they are removed\>* - e.g., \[19.8677, 34.7286\]
* *\<node numbers in the original graph making up this cycle\>* - e.g., \[2801, 2387\]. This can be useful for visualizing the path in tools like [Bandage](https://rrwick.github.io/Bandage/) -->

