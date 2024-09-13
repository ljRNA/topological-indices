import subprocess
from collections import deque
import numpy as np
import networkx as nx
from topological_indices.config import nauty_path


def linear_tree(N):
    return nx.path_graph(N + 1)


def star_tree(N):
    return nx.star_graph(N)


def second_star_tree(N):
    p = int(np.ceil(N / 2))
    return nx.Graph(
        [[0, i] for i in range(1, p + 1)] + [[i, i + p] for i in range(1, N - p + 1)]
    )


def dumbbell_tree(N, k):
    assert k <= N - 4, "k must be <= N - 4"
    s1 = int(np.floor((N - k) / 2))
    return nx.Graph(
        [[0, i] for i in range(1, s1 + 1)]
        + [[s1 + i, s1 + i + 1] for i in range(0, k + 2)]
        + [[s1 + k + 1, s1 + k + 1 + i] for i in range(2, N - s1 - k)]
    )


def starlike_tree(N, p):
    p = N - p
    return nx.Graph(
        [[0, i] for i in range(1, p + 1)] + [[i, i + p] for i in range(1, N - p + 1)]
    )


def spider_tree(N):
    return starlike_tree(N, int(N / 2))


def dandellion_tree(N, k):
    return nx.Graph(
        [[0, i] for i in range(1, N - k + 1)] + [[i, i + 1] for i in range(N - k, N)]
    )


def kary_tree(N, k):
    return nx.full_rary_tree(k, N + 1)


def bethe_tree(N, k):
    edges = []
    i = 1
    queue = deque([i])
    newqueue = deque()
    first = True

    while i <= N:
        if len(queue) != 0:
            curr = queue.pop()
            for _ in range(k + 1 if first else k):
                i += 1
                edges.append((curr, i))
                newqueue.append(i)
                if i == N + 1:
                    break
            first = False
        else:
            queue = newqueue.copy()
            newqueue.clear()

    return nx.Graph(edges)


def cayley_tree(N, k):
    edges = []
    left = N
    leaves = [0]
    j = 0
    level = 1
    while left > 0:
        N_previous = k ** (level - 1)
        new_leaves = []
        for i in range(min(left, k**level)):
            j += 1
            edges.append((leaves[i % N_previous], j))
            new_leaves.append(j)
        left -= len(new_leaves)
        level += 1
        leaves = new_leaves

    return nx.Graph(edges)


def flower_tree(N, k):
    assert N % k == 0, "N has to be divisible by k"

    m = int(N / k)
    edges = [
        (max((i - 2) * m + j, 0), (i - 1) * m + j)
        for i in range(1, k + 1)
        for j in range(1, m + 1)
    ]

    return nx.Graph(edges)


def is_chemical(T):
    return max(d for _, d in T.degree) <= 4


def is_caterpillar(T):
    deg = dict(T.degree)
    for n, d in deg.items():
        if d >= 3:
            i = 0
            for n_ in T.adj[n]:
                if deg[n_] > 1:
                    i += 1
                if i > 2:
                    return False
    return True


def is_lobster(T):
    T2 = T.copy()
    T2.remove_nodes_from([n for n in T if T.degree(n) == 1])
    return is_caterpillar(T2)


def all_trees(N):
    return list(nx.nonisomorphic_trees(N + 1))


def random_tree(N, seed=None):
    return random_trees(N, 1, seed)[0]


def random_BA_tree(N, seed=None):
    return random_BA_trees(N, 1, seed)[0]


def random_trees(N, number_of_trees=1, seed=None):
    Ts = nx.random_unlabeled_tree(N + 1, number_of_trees=number_of_trees, seed=seed)

    for T in Ts:
        nx.set_edge_attributes(T, {e: {"length": 1} for e in T.edges()})

    return Ts


def random_BA_trees(N, number_of_trees=1, seed=None):
    Ts = [nx.barabasi_albert_graph(N + 1, 1, seed=seed) for i in range(number_of_trees)]

    for T in Ts:
        nx.set_edge_attributes(T, {e: {"length": 1} for e in T.edges()})

    return Ts


def all_ncyclic_graphs(N, n):
    with subprocess.Popen(
        [nauty_path / "geng", "-c", f"{N-n+1}", f"{N}:{N}"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
    ) as geng:
        output, error = geng.communicate()

        assert error is None, error

        output = output.decode().strip("\n").split("\n")

        graphs = []
        for out in output:
            G = nx.from_graph6_bytes(out.encode())
            assert G.number_of_edges() == N

            graphs.append(G)

        return graphs
