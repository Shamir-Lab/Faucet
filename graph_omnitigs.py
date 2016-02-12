from utils import *


def extend_walk(walk, G, W):
	# assumes walk is list of nodes
	if walk==None: return
	print walk
	v_t = walk[-1]
	extended = False
	extensions = G.out_edges(v_t)
	for ext in extensions:
		INs = G.in_edges(walk)
		X = [a[0] for a in INs] # sources
		# nodes entering first node s in walk (s_ins)
		# should not be in X
		s_ins = G.predecessors(walk[0]) 
		# s_ins = [a[0] for a in s_ins]
		X = set(X)
		X -= set(walk)
		X -= set(s_ins)
		G_p = G.copy()
		G_p.remove_edge(ext[0],ext[1])
		# SP is returned as dict; keys are targets, 
		# vals are lists of nodes in shortest path
		# assume only keys for reachable targets
		SP = nx.shortest_path(G_p,source=v_t)
		# print SP
		if all([target not in SP for target in X]):
			# print "got to extend", walk + (ext[1],)
			extend_walk(walk + (ext[1],), G, W)
			extended = True
	if extended==False:
		W.add(walk)

def get_omnitigs(G):
	""" given assembly graph G
		calls extend_walk on each edge to get 
		omnitigs, doesn't include removal of non-maximal
		omnitigs or initial Y-to-V reduction
	"""
	W = set([])
	for e in G.edges():
		extend_walk(e,G,W)
	return W

def get_sample_graph_comp_seqs():
	# load test graph
	fastg = "assembly_graph.fastg"
	test_node = "EDGE_1243_length_1496_cov_78.6919"
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

# read/build graph
G,COMP,SEQS = get_sample_graph_comp_seqs()
# get & print omnitigs
W = get_omnitigs(COMP)
print "omnitigs are:"
for tig in W:
	print tig
