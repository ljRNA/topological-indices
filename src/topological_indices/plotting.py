import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from topological_indices import topological_indices, tree_generation, topological_indices

def latex_on():
    custom_preamble = {
        "text.usetex": True,
        "text.latex.preamble": r"\usepackage{amsmath}",
        "font.family": "serif",
        'font.size': 10
    }
    plt.rcParams.update(custom_preamble)


def latex_off():
    custom_preamble = {
        "text.usetex": False,
        "text.latex.preamble": r"\usepackage{amsmath}",
    }
    plt.rcParams.update(custom_preamble)


def plot_tree(T, iterations=500, seed=None, **nx_draw_args):
    pos = nx.spring_layout(T, iterations=iterations, seed=seed)

    nx.draw(T, pos, **nx_draw_args)

    plt.gca().set_aspect('equal', 'datalim')

    

def plot_phase_space_boundary(N, inds, ax=None, normalization=False):
    if ax is None:
        ax = plt.gca()

    boundary_coarseness = 1 if N < 50 else int(N/10)
    Tss = [
        [tree_generation.linear_tree(N)], 
        [tree_generation.linear_tree(N)] + list(reversed([tree_generation.dumbbell_tree(N, i) for i in range(0, N-4, boundary_coarseness)])) + [tree_generation.star_tree(N)], 
        [tree_generation.star_tree(N)], 
        [tree_generation.star_tree(N)] + list(reversed([tree_generation.cayley_tree(N, i) for i in range(int(N/2), N, boundary_coarseness)])),
        ]
    Tss.append([tree_generation.starlike_tree(N, int(N/2)), Tss[-1][-1]])

    if 'H2' in inds and 'ABC' in inds:
        Tss.append([Tss[-2][-1], tree_generation.linear_tree(N)])

    Tss_colors = ['tab:blue', 'tab:green', 'tab:red', 'tab:orange', 'tab:brown', 'tab:purple']
    plot_markers = {'tab:red': (8,1,0), 'tab:blue': 's', 'tab:orange': 'X', 'tab:green': 'h', 'tab:purple': '>', 'tab:brown': 'D'}

    for Ts, color in zip(Tss, Tss_colors):
        res = np.array([list(topological_indices.calculate_TIs(T, inds).values()) for T in Ts])
        
        zorder = 10000 if color in ['tab:red', 'tab:blue', 'tab:brown'] else None
        ms = 8 if color in ['tab:red', 'tab:blue', 'tab:brown'] else 5
        ms = 12 if color == 'tab:red' else ms

        if normalization == 'simple':
            res_orig = res.copy()
            res[:,0] = [topological_indices._normalize_simple(TI, inds[0], N) for TI in res_orig[:,0]]
            res[:,1] = [topological_indices._normalize_simple(TI, inds[1], N) for TI in res_orig[:,1]]
        elif normalization == 'cdf':
            res_orig = res.copy()
            res[:,0] = [topological_indices._normalize_cdf(TI, inds[0], N) for TI in res_orig[:,0]]
            res[:,1] = [topological_indices._normalize_cdf(TI, inds[1], N) for TI in res_orig[:,1]]

        ax.plot(res[:,0], res[:,1], '-', color=color, ms=ms, marker=plot_markers[color], zorder=zorder)

    if normalization:
        ax.set_aspect('equal')
