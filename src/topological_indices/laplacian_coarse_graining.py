import logging
import numpy as np
import networkx as nx
from scipy.linalg import expm
from scipy.optimize import brenth

logger = logging.getLogger('laplacian_coarse_graining')

def _find_between_community_edges(g, partition):
    edges = {}

    for (ni, nj) in g.edges():
        ci = partition[ni]
        cj = partition[nj]

        if ci != cj:
            # This is relevant only for graphs
            # try:
            #   edges[(ci, cj)] += [(ni, nj)]
            #   edges[(cj, ci)] += [(ni, nj)]
            # except KeyError:
            #   pass
            edges[(ci, cj)] = (ni, nj)
            edges[(cj, ci)] = (ni, nj)

    return edges


def _laplacian_coarse_grained_tree(T, tau, weight=None):
    T = nx.convert_node_labels_to_integers(T)

    # Construct Laplacian
    L = np.diag([d for _, d in T.degree()]) - \
        nx.adjacency_matrix(T, weight=weight)

    # Construct network propagator
    K = expm(-tau*L)

    # Construct partition function
    Z = np.trace(K)

    # Construct canonical density matrix
    rho = K/Z

    rhoprime = rho/np.minimum.outer(rho.diagonal(), rho.diagonal())
    zeta = np.heaviside(rhoprime-1, 1)

    meta_graph = nx.from_numpy_array(zeta)

    node_communities = {}
    communities = list(nx.connected_components(meta_graph))
    cgT = nx.Graph()
    for i, c in enumerate(communities):
        cgT.add_node(i, orig_nodes=c, size=len(c))
        for node in c:
            node_communities[node] = i
    node_communities = dict(sorted(node_communities.items()))

    edges = _find_between_community_edges(T, node_communities)

    communities_connecting_nodes = {i: set() for i in range(len(communities))}
    for (c1, c2), (n1, n2) in edges.items():
        cgT.add_edge(c1, c2, length=1)

        communities_connecting_nodes[c1].add(n1)
        communities_connecting_nodes[c2].add(n2)

    return cgT


def _laplacian_coarse_grained_tree_N(T, target_N, weight=None):
    def _lcgN(log_tau):
        tau = np.exp(log_tau)
        cgT = _laplacian_coarse_grained_tree(T, tau, weight=weight)

        return cgT.number_of_edges() - target_N

    # TODO: Better initial guess
    if target_N <= 20:
        log_tau_interval = 2, 3.5
    elif target_N <= 50:
        log_tau_interval = 1, 2.5
    else:
        log_tau_interval = 0.5, 1.5

    try:
        log_tau, res = brenth(_lcgN, *log_tau_interval, full_output=True)
    except:
        # The interval was too small, try a wider one.
        log_tau, res = brenth(_lcgN, -10, 10, full_output=True)

    logger.debug('brenth: %s', res)

    cgT = _laplacian_coarse_grained_tree(T, np.exp(log_tau), weight=weight)

    if cgT.number_of_edges() != target_N:
        logger.info('Could not coarse grain to N = %i. Returning N = %i instead.',
                    target_N, cgT.number_of_edges())

    return cgT


def coarse_grain(T: nx.Graph, target_N: int):
    """Return LRG coarse-grained tree with _approximately_ target_N edges.

    Args:
        T (nx.Graph): Tree
        target_N (int): Target number of edges. 
        weight (_type_, optional): _description_. Defaults to None.

    Returns:
        nx.Graph: Coarse-grained tree
    """
    return _laplacian_coarse_grained_tree_N(T, target_N)
