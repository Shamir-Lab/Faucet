from utils import *

def extend_walk(walk, G, W, W_reps):
	""" given walk as list of nodes, 
		extends recursively until omnitig property
		is violated (next extension e has 
		alternate e' that begins path to intermed. 
		node in current walk), or walk includes next
		extension already (avoids repeated walks)
	"""
	# assumes walk is list of nodes
	if walk==None: return
	edge_set = {(walk[i],walk[i+1]) for i in range(len(walk)-1)}
	# print walk
	v_t = walk[-1]
	extended = False
	extensions = G.successors(v_t)
	for ext in extensions:
		INs = G.in_edges(walk[1:])
		X = [a[0] for a in INs] # sources
		X = set(X)
		X -= set(walk)
		G_p = G.copy()
		G_p.remove_edge(v_t,ext)
		# SP is returned as dict; keys are targets, 
		# vals are lists of nodes in shortest path
		# assume only keys for reachable targets
		SP = nx.shortest_path(G_p,source=v_t)

		if all([target not in SP for target in X]) and \
		(v_t, ext) not in edge_set:
			extend_walk(walk + (ext,), G, W, W_reps)
			extended = True
	# avoid cyclic rotationgs of same cycle - include 
	# all but re-entered last node in walk rep:
	walk_rep = get_unoriented_sorted_str(walk[:-1])
	if extended==False and walk_rep not in W_reps:
		W.add(walk)
		W_reps.add(walk_rep)
		return

def get_omnitigs(G):
	""" given assembly graph G
		calls extend_walk on each edge to get 
		omnitigs, doesn't include removal of non-maximal
		omnitigs or initial Y-to-V reduction
	"""
	W = set([])
	W_reps = set([]) # unique reps to avoid cyclic rotations
	for e in G.edges():
		extend_walk(e,G,W,W_reps)
	return W

def get_maximal_omnitigs(W):
	""" removes any omnitigs completely contained
		in some other omnitig
	"""
	max_tigs = set([])
	for w in W:
		# want to get rid of all walks having all nodes 
		# or all rcs of nodes in some other walk
		rc_path = [rc_node(a) for a in w]
		if any([(set(w) < set(x)) or (set(rc_path) < set(x)) for x in W]):
			continue
		else:
			max_tigs.add(w)
	return max_tigs	

def get_sample_graph_comp_seqs(fastg, test_node):
	# load test graph
	# fastg = "JJ1886_graph.fastg"
	# test_node = "EDGE_1243_length_1496_cov_78.6919"
	G = get_fastg_digraph(fastg)
	comps = nx.strongly_connected_component_subgraphs(G)
	COMP = nx.DiGraph()

	# choose desired SCC based on node in it
	for c in comps:
		if test_node in c.nodes():
			COMP = c.copy()
			break
	SEQS = get_fastg_seqs_dict(fastg, G)
	return G,COMP,SEQS


def peel_omnitigs(G, W, C):
	""" given component G and a set of maximal omnitigs W
		peels omnitigs by subtracting genome coverage C 
		until no covered nodes remain
	"""

# read/build graph
# G,COMP,SEQS = get_sample_graph_comp_seqs("JJ1886_graph.fastg", "EDGE_1243_length_1496_cov_78.6919") #"EDGE_216_length_72_cov_28.2941")#
# G,COMP,SEQS = get_sample_graph_comp_seqs("E2022_graph.fastg", "EDGE_286_length_92_cov_109.162")
G,COMP,SEQS = get_sample_graph_comp_seqs("KPN_graph.fastg", "EDGE_52_length_308692_cov_66.7584") #"EDGE_137_length_40951_cov_36.5289") #

# get omnitigs
W = get_omnitigs(COMP)
W2 = get_maximal_omnitigs(W)
print "maximal omnitigs are:"
for tig in W2:
	print tig
