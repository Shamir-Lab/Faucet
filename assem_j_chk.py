from pybloomfilter import BloomFilter

read_len = 100
k = 27
fp  = 0.01
j = 10
bases = ['A', 'C', 'G', 'T']


def get_kmers(r,k):
	return [r[i:i+k] for i in range(len(r)-k+1)]

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
				if test_kmer in bf:
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


def get_back_suffixes(buff, j):
	""" fetches dictionary including all k-j suffixes
		from first level (back) of buffer as keys and the 
		k-mers including them as values
	"""
	backs = list(buff[0])
	suffs = {}
	for b in backs:
		if b[j:] in suffs:
			suffs[b[j:]].append(b)
		else:
			suffs[b[j:]]=[b]
	return suffs

def get_front_prefixes(buff, j):
	""" fetches dictionary including all k-j prefixes
		from last level (front) of buffer as keys and the 
		k-mers including them as values
	"""
	fronts = list(buff[-1])
	prefs = {}
	for f in fronts:
		if f[j:] in prefs:
			prefs[f[:-j]].append(f)
		else:
			prefs[f[:-j]]=[f]
	return prefs

def pretty_print_buffer(buff):
	for level in buff:
		print list(level)
	print "\n"

def load_bf_sources_sinks(filename,j):
	""" loads k-mers into bloom filter
		loads to lists ends of reads as potential
		sources and sinks (or nodes j away from sinks)
	"""
	sources = []
	j_sinks = []
	B = BloomFilter(capacity = 500000000, error_rate=fp)
	line_no = 0
	with open(filename) as f:
		for line in f:
			if line_no%1000000==0:
				print line_no
			read = line.rstrip()
			sources.append(read[:k])
			j_sinks.append(read[-(k+j):-j])
			B.update(get_kmers(read,k))
			line_no +=1
	return (B,sources,j_sinks)

def clean_alt_paths_from_buff(alts, backs, fronts, buff):
	""" given (k-j)-mers of alt paths, gets their start and end
		k-mers to create paths to check (returned), removes k-mers including
		them from buffer front (modified)
	"""
	paths = []
	# front = buff[-1]
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
		true branch points
	"""
	cands = []
	line_no = 0
	with open(filename) as f:
		for line in f:
			if line_no%1000==0:
				print line_no
			read = line.rstrip()
			kmers = get_kmers(read,k)
			buff = get_j_forward_buff(kmers[0],B,j)
			
			for ind, kmer in enumerate(kmers):
				if len(buff[-1])>1 and len(buff[0])>1:
					print "read no: %d, position: %d" % (line_no, ind)
					backs = get_back_suffixes(buff,j)
					fronts = get_front_prefixes(buff,j)
					comms = (set(backs.keys())).intersection(set(fronts.keys()))
					if ind == read_len-k:
						break
					next_real = kmers[ind+1][j:] # k-j suffix from next real k-mer						
					alts = comms - set([next_real]) 
					print "back suffixes, back-mers: ", list(backs), backs.values()
					print "front prefixes, front-mers: ", list(fronts), fronts.values()
					print "commons: ", list(comms)
					print "real next: ", next_real
					print "kmer: ", kmer
					print "alts: ", list(alts)
					alt_paths = clean_alt_paths_from_buff(alts, backs, fronts, buff)
					print "paths to check: ", alt_paths
					# pretty_print_buffer(buff)
				advance_buffer(buff,B)
				
			line_no +=1
			if line_no==10:
				break


####### main ####### 
fname1 = "/home/nasheran/rozovr/BARCODE_test_data/chr20.c10.reads.head"
(B,sources,sinks) = load_bf_sources_sinks(fname1,j)
get_candidate_false_joins(fname1,B)




# need to also change sinks, sources from candidates to 
# real ends

# for s in sources:
# buff = get_level_j_buff(sources[0],B,j)
# # (next_step, j_step) = advance_with_buffer(s,j,B)
# print buff