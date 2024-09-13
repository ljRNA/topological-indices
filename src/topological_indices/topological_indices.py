import numpy as np
import networkx as nx
import scipy.stats as st
from topological_indices import utilities

try:
    import pynauty
except:
    pynauty = None


INDICES = ["N", "MLD", "ALD", "NA", "λ2", "λN", "M1",
           "M2", "W", "P", "ABC", "R", "SC", "H1", "H2", "J", "EE", "B"]
INDICES_CDF = ["ABC", "H2", "W", "R"]

if pynauty is not None:
    INDICES.append("I")


dfFitRandom = utilities.load_df('random_trees/dfFitRandom-H2-ABC-W-R')
dfFitParameters = utilities.load_df('random_trees/dfFitParameters')
dfLinearStar = utilities.load_df('all_trees/dfLinearStar')

dist_names = {'H2': 'moyal', 'ABC': 'lognorm',
              'W': 'pearson3', 'R': 'pearson3'}


def calculate_TIs(T: nx.Graph, inds=None, normalization=None):
    if inds is None:
        if normalization is None or normalization == 'simple':
            inds = INDICES
        else:
            inds = INDICES_CDF + ['N']

    if type(inds) != list:
        inds = [inds]

    TIs = _topological_indices(T, inds)

    if normalization is not None:
        TIs = normalize_TIs(TIs, normalization)

    return TIs


def normalize_TIs(TIs: dict, normalization, N=None):
    if N is None:
        N = int(TIs['N'])

    if normalization == 'simple':
        TIs_normalized = {
            ind + "'": _normalize_simple(TI, ind, N) for ind, TI in TIs.items() if ind != 'N'}

    elif normalization == 'cdf':
        TIs_normalized = {
            ind + "*": _normalize_cdf(TI, ind, N) for ind, TI in TIs.items() if ind in INDICES_CDF}

    return TIs_normalized


def _topological_indices(T, indices=None):
    N = T.number_of_edges()

    assert N != 0, "Topological indices are ill-defined for single-node trees (N = 0)."

    if set(['MLD', 'ALD', 'W', 'J', 'P']) & set(indices):
        distance_matrix = nx.floyd_warshall_numpy(T)
        distances = distance_matrix[np.triu_indices_from(
            distance_matrix, k=1)].flatten()

    if 'J' in indices:
        DS = {node: DSi for node, DSi in zip(
            T.nodes(), np.sum(distance_matrix, axis=0))}

    if set(['λ2', 'λN', 'H1', 'H2', 'logλ2', 'logλN']) & set(indices):
        laplacian_spectrum = nx.laplacian_spectrum(T)

    if 'EE' in indices:
        adjacency_spectrum = nx.adjacency_spectrum(T)

    res = {}
    for ind in indices:
        if ind == 'N':
            val = T.number_of_edges()
        elif ind == 'MLD':
            val = int(np.max(distances))
        elif ind == 'ALD':
            val = np.mean(distances)
        elif ind == 'λ2':
            val = laplacian_spectrum[1]
        elif ind == 'λN':
            val = laplacian_spectrum[N]
        elif ind == 'logλ2':
            val = np.log(laplacian_spectrum[1])
        elif ind == 'logλN':
            val = np.log(laplacian_spectrum[N])
        elif ind == 'NA':
            val = int((1/2)*np.sum([d*(d-1) for (_, d) in T.degree()]))
        elif ind == 'M1':
            val = np.sum(np.array([d for (_, d) in T.degree()])**2, dtype=int)
        elif ind == 'M2':
            val = np.sum([T.degree(n1)*T.degree(n2)
                         for n1, n2 in T.edges], dtype=int)
        elif ind == 'W':
            val = np.sum(distances)
        elif ind == 'P':
            val = np.count_nonzero(distances == 3)
        elif ind == 'ABC':
            val = np.sum([np.sqrt((T.degree(n1)+T.degree(n2)-2) /
                         (T.degree(n1)*T.degree(n2))) for n1, n2 in T.edges])
        elif ind == 'R':
            val = np.sum([1/np.sqrt(T.degree(n1)*T.degree(n2))
                         for n1, n2 in T.edges])
        elif ind == 'SC':
            val = np.sum([1/np.sqrt(T.degree(n1)+T.degree(n2))
                         for n1, n2 in T.edges])
        elif ind == 'H1':
            val = (1/(2*N))*np.sum(1/laplacian_spectrum[1:])
        elif ind == 'H2':
            val = (1/(2*N))*np.sum(1/laplacian_spectrum[1:]**2)
        elif ind == 'J':
            val = (N/(N-T.number_of_nodes()+2)) * \
                np.sum([1/np.sqrt(DS[n1]*DS[n2]) for n1, n2 in T.edges])
        elif ind == 'EE':
            val = np.sum(np.exp(adjacency_spectrum)).real
        elif ind == 'B':
            B = 0
            T2 = T.copy()

            while len(T2.nodes()) > 1:
                c = 0
                g = T2.copy()
                for node in g.nodes():
                    if nx.degree(g, node) == 1:
                        c += 1
                        T2.remove_node(node)
                B += c**2
            val = B
        elif ind == 'I':
            if pynauty is not None:
                T2 = nx.convert_node_labels_to_integers(T)
                g = pynauty.Graph(T2.number_of_nodes())
                for u, v in T2.edges():
                    g.connect_vertex(u, v)

                # generators, grpsize1, grpsize2, orbits, numorbits
                _, _, _, orbits, _ = pynauty.autgrp(g)

                _, p = np.unique(orbits, return_counts=True)
                assert np.sum(p) == T2.number_of_nodes()
                p = p/np.sum(p)

                val = -np.sum(p*np.log2(p))
            else:
                val = None

        elif ind == 'D2':
            val = np.count_nonzero(np.array([d for _, d in T.degree]) == 2)
        elif ind == 'min(D)':
            val = min(d for _, d in T.degree if d > 1)
        elif ind == 'max(D)':
            val = max(d for _, d in T.degree)

        res[ind] = val

    return res


def _TI_linear_star(N, ind):
    pts = dfLinearStar.loc[[f'linear-{N}', f'star-{N}'], ind].values

    return pts


def _interpolate_fit(N, ind, param_i, alpha, beta):
    if ind == 'W' and param_i == 0:
        return alpha/N + beta
    elif ind == 'R' and param_i == 0:
        return -np.exp(alpha)*(N**beta)
    else:
        return np.exp(alpha)*(N**beta)


def _normalize_simple(TI, ind, N):
    TI_linear, TI_star = _TI_linear_star(N, ind)

    return (TI - TI_linear) / (TI_star - TI_linear)


def _normalize_cdf(TI, ind, N):
    dist_name = dist_names[ind]
    distribution = getattr(st, dist_name)

    params = dfFitRandom[(dfFitRandom['ind'] == ind) & (
        dfFitRandom['distribution'] == dist_name) & (dfFitRandom['N'] == N)]['parameters'].values

    if len(params) == 0:
        # Interpolation of parameters
        row = dfFitParameters[(dfFitParameters['ind'] == ind) & (
            dfFitParameters['distribution'] == dist_name)].sort_values('param')
        params = []
        for _, row in row.iterrows():
            params.append(_interpolate_fit(
                N, ind, row['param'], row['alpha'], row['beta']))
    else:
        params = params[0]

    TI_linear, TI_star = _TI_linear_star(N, ind)

    args = params[:-2]
    loc = params[-2]
    scale = params[-1]

    def cdf(TI):
        return distribution.cdf(TI, loc=loc, scale=scale, *args)

    TIn_linear = cdf(TI_linear)
    TIn_star = cdf(TI_star)
    TIn = cdf(TI)

    return (TIn - TIn_linear) / (TIn_star - TIn_linear)
