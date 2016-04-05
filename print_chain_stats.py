from utils import *

def get_chain_forward(start_node, G):
    # node = start_node
    # chain = [[node]]
    # pos = 0
    # while((pos%2==0 and G.out_degree(node)==2) or \
    #     (pos%2==1 and G.neighbors(chain[-1][0])==G.neighbors(chain[-1][1]) \
    #         and len(G.neighbors(chain[-1][0]))==1)):
    neighbors = G.neighbors(start_node)
    if G.out_degree(start_node)==2 and \
    G.neighbors(neighbors[0])==G.neighbors(neighbors[1]) and \
    len(G.neighbors(neighbors[0]))==1:
        covs = [get_cov_from_spades_name(n) for n in neighbors]
        ratio = max(covs[0]/covs[1], covs[1]/covs[0])
        print get_length_from_spades_name(neighbors[0]), ratio
        return 1
    else:
        return 0

    


    #     neighbors = G.neighbors(node) # list including either a single node or a pair
    #     # used_nodes.update(set(neighbors))
    #     chain.append(neighbors) 
    #     node = neighbors[0] # arbitrary choice in pair (when a pair)
    #     pos += 1
    # bubble_cnt = 0
    # if len(chain)>1:
    #     for nodes in chain[1::2]:
    #         bubble_cnt += 1
    #         covs = [get_cov_from_spades_name(n) for n in nodes]
    #         ratio = max(covs[0]/covs[1], covs[1]/covs[0])
    #         print get_length_from_spades_name(nodes[0]), ratio

    # return bubble_cnt
            # print 
    #     return None, used_nodes
    # else:
    #     print chain
    #return chain, used_nodes

# chains = []

fastg_mink = "/home/nasheran/rozovr/mink_data/JJ1886.cleaned_graph.fastg"
fastg_spades = "/home/nasheran/rozovr/recycle_paper_data/JJ1886_assem/spades3.6.2/assembly_graph.fastg"
G = get_fastg_digraph(fastg_mink)

# used_nodes = set([])
count = 0
for node in G.nodes():
    # if node in used_nodes: continue
    cnt = get_chain_forward(node, G)
    # if chain: chains.append(chain)
    count += cnt
print count
    # if count >= 100: break
# print chains
