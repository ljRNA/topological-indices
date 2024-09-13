import logging
import networkx as nx
from networkx.algorithms.community import greedy_modularity_communities, modularity

logger = logging.getLogger('modularity_coarse_graining')

def _community_tree(T, communities):
    T2 = T.copy()
    cgT = nx.Graph()

    node_communities = {}
    for i, c in enumerate(communities):
        for node in c:
            node_communities[node] = i

        cgT.add_node(i, orig_nodes=c, size=len(c))

    nx.set_node_attributes(T2, {n: {'community': i}
                           for n, i in node_communities.items()})

    communities_connecting_nodes = {i: [] for i in range(len(communities))}
    for n1, n2, data in T2.edges(data=True):
        c1, c2 = node_communities[n1], node_communities[n2]
        if c1 != c2:
            try:
                cgT.add_edge(c1, c2, length=data['length'])
            except:
                cgT.add_edge(c1, c2, length=1)

            communities_connecting_nodes[c1].append(n1)
            communities_connecting_nodes[c2].append(n2)

    return cgT


def _find_communities_N(T, target_N: int, weight=None):
    communities = greedy_modularity_communities(
        T, weight=weight, cutoff=target_N+1, best_n=target_N+1)

    return communities, modularity(T, communities)


def coarse_grain(T: nx.Graph, target_N: int):
    communities, Q = _find_communities_N(T, target_N)

    return _community_tree(T, communities)
