import random
import numpy as np
from pybloomfilter import BloomFilter

read_len = 100

# read in reads
reads = []
fname1 = "/home/nasheran/rozovr/BARCODE_test_data/chr20.c10.reads"
with open(fname1) as f:
	for line in f:
		read = line.rstrip()
		reads.append(read)

reads_set = set(reads) # for fast validation


# read in chromosome fasta, 
# optional (not done): remove all lines containing 
# something other than ACGT
fname2 = "/specific/a/home/cc/cs/rozovr/chr20.fa"
with open(fname2) as f:
    lines = f.readlines()

seq = ''.join(lines[1:])
seq = seq.replace('\n','')

print "finished reading reads and genome"

fp  = 0.01
B = BloomFilter(capacity = 6500000, error_rate=fp)
# load reads to BF
B.update(reads)

print "loaded BF"
# scan genome, query B with read_len windows 
pos = 0
TP_cnt = 0
FP_cnt = 0
FP2_cnt = 0
TN_cnt = 0
num_windows = len(seq)-read_len+1 
while (pos < num_windows):
	if pos%3000000==0:
		print pos, TN_cnt, TP_cnt, FP_cnt, FP2_cnt
	test_seq1 = seq[pos:pos+read_len]
	test_seq2 = seq[pos+1:pos+1+read_len]
	if (test_seq1 in B):
		if (test_seq1 not in reads_set):	
			FP_cnt += 1
			if (test_seq2 not in B):
				FP2_cnt += 1
		else:
			TP_cnt += 1
	else:
		TN_cnt += 1

	pos += 1

print fp, float(FP_cnt)/num_windows, float(FP2_cnt)/num_windows