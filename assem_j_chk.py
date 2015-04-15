from pybloomfilter import BloomFilter

read_len = 100
k = 27
fp  = 0.01
j = 10
bases = ['A', 'C', 'G', 'T']
complements = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

def get_kmers(r,k):
	return [r[i:i+k] for i in range(len(r)-k+1)]

def get_rc(dna):
	rev = reversed(dna)
	return "".join([complements[i] for i in rev])

def get_canons(kmers):
	return [min(x,get_rc(x)) for x in kmers]

def get_j_forward_buff(source,bf,depth):
	"""initialize traversal buffer
		stores source + nodes up to depth accepted by bf
		in list of sets
	"""
	buff = []
	roots = [source]
	buff.append(set(roots))
	next = []
	for level in range(depth+1):
		for root in roots:
			for b in bases:
				test_kmer = root[1:]+b
				canon = min(test_kmer, get_rc(test_kmer))
				if canon in bf:
					next.append(test_kmer)
		roots = next[:]
		buff.append(set(roots))
		next = []
	del buff[0]
	return buff


def advance_buffer(buff, bf):
	""" does BFS using Bloom filter bf from loaded buffer
		extends each front node one step, removes first level
	"""
	fronts = buff[-1]
	next = []
	for node in list(fronts):
		for b in bases:
			test_kmer = node[1:]+b
			if test_kmer in bf:
				next.append(test_kmer)
	buff.append(set(next))
	del buff[0] # modify original instead of copying


def get_buffer_level(buff, j, level):
	""" gets buffer contents from chosen level (in [0,j])
		returns dictionary containing invariant k-j as keys
		and k-mers including them as values - e.g., level = 0
		contains (k-j)-mers as suffixes, level = j+1 contains 
		them as prefixes  
	"""
	kmers = list(buff[level])
	invars = {}
	for kmer in kmers:
		if level != 0:
			inv = kmer[j - level : -level]
		else:
			inv = kmer[j - level :]
			 
		if inv in invars:
			invars[inv].append(kmer)
		else:
			invars[inv]=[kmer]
	return invars

def get_back_suffixes(buff, j):
	""" fetches dictionary including all k-j suffixes
		from first level (back) of buffer as keys and the 
		k-mers including them as values
	"""
	backs = list(buff[0])
	suffs = {}
	for b in backs:
		suff = b[j:]
		if suff in suffs:
			suffs[suff].append(b)
		else:
			suffs[suff]=[b]
	return suffs

def get_front_prefixes(buff, j):
	""" fetches dictionary including all k-j prefixes
		from last level (front) of buffer as keys and the 
		k-mers including them as values
	"""
	fronts = list(buff[-1])
	prefs = {}
	for f in fronts:
		pref = f[:-j]
		if pref in prefs:
			prefs[pref].append(f)
		else:
			prefs[pref]=[f]
			# prefs[f[:-j]]=[f]
	return prefs

def pretty_print_buffer(buff):
	for level in buff:
		print list(level)
	print "\n"

def load_bf_sources_sinks(filename,j,numreads):
	""" loads k-mers into bloom filter
		loads to lists ends of reads as potential
		sources and sinks (or nodes j away from sinks)
	"""
	sources = []
	j_sinks = []
	reals = []
	B = BloomFilter(capacity = numreads * (read_len-k), error_rate=fp)
	line_no = 0
	with open(filename) as f:
		for line in f:
			# if line_no%1000000==0:
			# 	print line_no
			read = line.rstrip()
			sources.append(read[:k])
			sources.append(get_rc(read[-k:]))
			j_sinks.append(read[-(k+j):-j])
			j_sinks.append(get_rc(read[j:j+k]))
			kmers = get_kmers(read,k)
			# insert canonical (lex-min) kmers only
			canons = get_canons(kmers) 
			B.update(canons)
			reals.extend(canons)
			line_no +=1
	reals = set(reals)
	sources = set(get_canons(sources))
	j_sinks = set(get_canons(j_sinks))
	return (B,sources,j_sinks,reals)

def clean_alt_paths_from_buff(alts, backs, fronts, buff):
	""" given (k-j)-mers of alt paths, gets their start and end
		k-mers to create paths to check (list returned), removes k-mers including
		them from buffer front (modified) and one level ahead of buffer back
	"""
	paths = []
	for a in alts:
		pref = backs[a][0]
		ends = fronts[a]
		for end in ends:
			path = pref + end[-j:]
			paths.append(path)
			buff[-1].discard(end)
	return paths


def get_candidate_false_joins(filename,bf):
	""" scan reads to find candidate false joins. 
		finds nodes having descendents at level j that differ from read sequence
		used to check against reals later, to know which are false joins, vs. 
		true branch points; returns list of candidate lists corr. to all candidate from a read
	"""
	cands = []
	line_no = 0
	with open(filename) as f:
		for line in f:
			if (line_no+1)%10000==0:
				print line_no+1, len(cands)
			read = line.rstrip()
			kmers = get_kmers(read,k)
			buff = get_j_forward_buff(kmers[0],bf,j)
			
			for ind, kmer in enumerate(kmers):
				if len(buff[-1])>1 and len(buff[0])>1:
					# print "read no: %d, position: %d, front len %d, back len %d" % (line_no, ind, len(buff[-1]), len(buff[0]))
					backs = get_back_suffixes(buff,j)
					fronts = get_front_prefixes(buff,j)
					comms = (set(backs.keys())).intersection(set(fronts.keys()))
					if ind == read_len-k:
						break
					next_real = kmers[ind+1][j:] # k-j suffix from next real k-mer						
					alts = comms - set([next_real]) 
					# print "back suffixes, back-mers: ", list(backs), backs.values()
					# print "front prefixes, front-mers: ", list(fronts), fronts.values()
					# print "commons: ", list(comms)
					# print "real next: ", next_real
					# print "kmer: ", kmer
					# print "alts: ", list(alts)
					# print "front len before cleaning: ", len(buff[-1])
					alt_paths = clean_alt_paths_from_buff(alts, backs, fronts, buff)
					if alt_paths:
						cands.append(alt_paths)
					# print "front len after cleaning: ", len(buff[-1])
					# print "paths to check: ", alt_paths
					# pretty_print_buffer(buff)
				advance_buffer(buff,bf)	
			line_no +=1
	
	print "got forward candidates"
	line_no = 0
	with open(filename) as f:
		for line in f:
			if (line_no+1)%10000==0:
				print line_no+1, len(cands)
			read = line.rstrip()
			kmers = get_kmers(get_rc(read),k)
			buff = get_j_forward_buff(kmers[0],bf,j)
			
			for ind, kmer in enumerate(kmers):
				if len(buff[-1])>1 and len(buff[0])>1:
					# print "read no: %d, position: %d, front len %d, back len %d" % (line_no, ind, len(buff[-1]), len(buff[0]))
					backs = get_back_suffixes(buff,j)
					backs2 = get_buffer_level(buff,j,0)
					fronts = get_front_prefixes(buff,j)
					fronts2 = get_buffer_level(buff,j,j)
					comms = (set(backs.keys())).intersection(set(fronts.keys()))
					if ind == read_len-k:
						break
					next_real = kmers[ind+1][j:] # k-j suffix from next real k-mer						
					alts = comms - set([next_real]) 
					print "back suffixes, back-mers: ", list(backs), backs.values()
					print "back suffixes, back-mers: ", list(backs2), backs2.values()
					print "front prefixes, front-mers: ", list(fronts), fronts.values()
					print "front prefixes, front-mers: ", list(fronts2), fronts2.values()
					
					# print "commons: ", list(comms)
					# print "real next: ", next_real
					# print "kmer: ", kmer
					# print "alts: ", list(alts)
					# print "front len before cleaning: ", len(buff[-1])
					alt_paths = clean_alt_paths_from_buff(alts, backs, fronts, buff)
					if alt_paths:
						cands.append(alt_paths)
					#### to debug junctions
					# backs = get_back_suffixes(buff,j)
					# fronts = get_front_prefixes(buff,j)
					# print "back suffixes, back-mers: ", list(backs), backs.values()
					# print "front prefixes, front-mers: ", list(fronts), fronts.values()
					####

					# print "front len after cleaning: ", len(buff[-1])
					# print "paths to check: ", alt_paths
					# pretty_print_buffer(buff)
				advance_buffer(buff,bf)	
			line_no +=1

	print "got rc candidates"
	return cands

def check_path_for_false_joins(path, bf, reals):
	""" given a path in the form of a string (incl several k-mers)
		check if any k-mers in path are false [TODO: return true only 
		if F->T transition]
	"""
	kmers = get_kmers(path, k)
	canons = get_canons(kmers)
	bools = [kmer in reals for kmer in canons]
	if all(bools): #TODO: checking paths should be cached - k-mer subpaths often repeated
		return False
	else:
		# print bools
		if bools[-1]==True: # path contains some F and ends in T
			return True 

def find_real_ends(cands, hsh, fetch_juncs = False):
	""" given candidate sources, query hash (bf or reals) on all
		those that have no extension on 5' end are sources
		those that have no extension on 3' are sinks
		if using bf, to_chk are for check against reals to make sure no real ends
		are missed due to FPs.  If using reals, to_chk are non-ends
	"""
	sources = []
	sinks = []
	to_chk = []
	juncs = []
	for c in cands:
		num_ext_r = 0
		num_ext_l = 0
		pref = c[1:]
		suff = c[:-1]		
		for b in bases:
			test_kmer_r = pref+b
			canon = min(test_kmer_r, get_rc(test_kmer_r))
			if canon in hsh:
				num_ext_r += 1
			test_kmer_l = b+suff
			canon = min(test_kmer_l, get_rc(test_kmer_l))
			if canon in hsh:
				num_ext_l += 1

		if num_ext_l == 0:
			sources.append(c)
		elif num_ext_r == 0:
			sinks.append(c)
		elif fetch_juncs and (num_ext_l>1 or num_ext_r>1):
			juncs.append(c)
		else:
			to_chk.append(c)

	sources = set(sources)
	sinks = set(sinks)
	juncs = set(juncs)
	if fetch_juncs:
		return (sources,sinks,juncs)
	else:
		return (sources,sinks,to_chk)



####### main ####### 
reads_f = "/home/nasheran/rozovr/BARCODE_test_data/chr20.c10.reads.head"
(B,src_cnd,j_sinks,reals) = load_bf_sources_sinks(reads_f,j,10000)

# get and count junctions, false joins
bf_cands = get_candidate_false_joins(reads_f,B)
fj_cnt = 0
br_cnt = 0
junc_nodes = []
for cnd_lst in bf_cands:
	for c in cnd_lst:
		if check_path_for_false_joins(c,B,reals):
			fj_cnt += 1
		else:
			br_cnt += 1
			junc_nodes.append(min(c[:k],get_rc(c[:k])))
junc_nodes = set(junc_nodes)
print "bf cands", len(bf_cands), fj_cnt, br_cnt
print "junc_nodes", len(junc_nodes)

# sanity check - for debugging counts
# count junctions using reals set
# there should be no false joins
rl_cands = get_candidate_false_joins(reads_f,reals) 
fj_cnt = 0
br_cnt = 0
real_junc_nodes = []
for cnd_lst in rl_cands:
	for c in cnd_lst:
		if check_path_for_false_joins(c,B,reals):
			fj_cnt += 1
		else:
			br_cnt += 1
			real_junc_nodes.append(min(c[:k],get_rc(c[:k])))
real_junc_nodes = set(real_junc_nodes)
print "real cands", len(bf_cands), fj_cnt, br_cnt
print "real junc_nodes", len(real_junc_nodes)
print len(junc_nodes.intersection(real_junc_nodes))

# find real ends, get candidates for cheking
sources, sinks, to_chk = find_real_ends(list(src_cnd), B)
print "sources, sinks, to_chk"
print len(sources), len(sinks), len(to_chk)

# check remaining end candidates, join with first stage ends
sources2, sinks2, to_discard = find_real_ends(to_chk, reals)
print len(sources2), len(sinks2), len(to_discard)
sinks |= sinks2
sources |= sources2
print len(sources), len(sinks)


# sanity check - for debugging counts
# sources, sinks, juncs = find_real_ends(list(reals), reals, fetch_juncs=True)
# print "real sources, sinks, juncs"
# print len(sources), len(sinks), len(juncs)


# traverse from all sources, mark junctions visited, 
# output contigs