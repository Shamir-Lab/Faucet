import random
import numpy as np
from pybloomfilter import BloomFilter

global kmer_len 
kmer_len = 25

def calc_expected_fps(G_len, f):
	max_dist = np.ceil(np.log(1./G_len)/ np.log(4*f))
	num = 0
	for i in range(int(max_dist)):
		num += G_len * (4*f)**(i+1)
	return num	

def get_rand_kmer(k, bases):
	return ''.join(random.choice(bases) for _ in range(k))

def sanity_check(bf, fp_param, num_trials, bases):
	count = 0	
	for i in range(num_trials):
		kmer = get_rand_kmer(kmer_len, bases)
		if kmer in bf:
			count += 1
	return float(count)/num_trials

def load_hash(load_lst, hsh):
    for l in load_lst:
        hsh.add(l)
    return hsh

fname = "NODE_39_length_14469.fa"

with open(fname) as f:
    lines = f.readlines()

seq = ''.join(lines[1:])
seq = seq.replace('\n','')

j = 0
rk_lst = []
while j + kmer_len < len(seq):
	rk_lst.append(seq[j:j+kmer_len])
	j += 1
real_kmers = load_hash(rk_lst, set([]))


bases = ['A', 'C', 'G', 'T']
orig_roots = []
for k in rk_lst:
	root = True
	for b in bases:
		if b+k[:-1] in real_kmers:
			root = False
	if root:
		orig_roots.append(k)


fp_lst = [0.0001, 0.001, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3]

# BFS traverse, record
# maintain stack of leaves of graph

for fp in fp_lst:
	exp_fp = calc_expected_fps(14443, fp)
	print "expected fps: %d" % exp_fp
	obs_fp = []
	san_fp = []
	queries = []
	fp_cnt = []
	for ind in range(10):	
		roots = orig_roots
		bf = BloomFilter(capacity = 14443, error_rate=fp)
		bf = load_hash(real_kmers, bf)
		num_queries = 0
		num_accepts = 0
		num_reals = 0
		num_fp = 0
		next = []
		level = 0
		while (len(roots) < 100000) and (level < 30000):
			for root in roots:
				for b in bases:
					test_kmer = root[1:]+b
					next_kmers = [test_kmer[1:]+x for x in bases]
					next_ins = [a in bf for a in next_kmers]
					if test_kmer in bf: #and any(next_ins):
						next.append(test_kmer)
						num_accepts += 1
						if test_kmer in real_kmers:
							num_reals += 1
						else:
							num_fp += 1
					num_queries += 1
			roots = next[:]
			next = []
			level += 1
	
		obs_fp.append(float(num_fp)/num_queries)
		san_fp.append(sanity_check(bf,fp,num_queries, bases))
		queries.append(num_queries)
		fp_cnt.append(num_fp)
	print "FP chosen: %.4f, mean FP observed: %f, std. obs: %f, mean FP sanity %f, std. san: %f, mean queries: %.1f, mean fps: %.1f" \
	% (fp, np.mean(obs_fp), np.std(obs_fp), np.mean(san_fp), np.std(san_fp), np.mean(num_queries), np.mean(fp_cnt)	)
	#level %d, accepts: %d, false pos: %d" \
	#, level, num_accepts, num_fp)


