from utils import *

def is_bubble(start_node, G):
    neighbors = G.neighbors(start_node)
    if G.out_degree(start_node)==2 and \
    G.neighbors(neighbors[0])==G.neighbors(neighbors[1]) and \
    len(G.neighbors(neighbors[0]))==1:
        covs = [get_cov_from_spades_name(n) for n in neighbors]
        ratio = max(covs[0]/covs[1], covs[1]/covs[0])
        if get_length_from_spades_name(neighbors[0])==55: print neighbors
        print get_length_from_spades_name(neighbors[0]), ratio
        return 1
    else:
        return 0

fastg_mink = "/home/nasheran/rozovr/mink_data/JJ1886.cleaned_graph.fastg"
fastg_spades = "/home/nasheran/rozovr/recycle_paper_data/JJ1886_assem/spades3.6.2/assembly_graph.fastg"
G = get_fastg_digraph(fastg_mink)

# used_nodes = set([])
count = 0
for node in G.nodes():
    cnt = is_bubble(node, G)
    count += cnt
print count
