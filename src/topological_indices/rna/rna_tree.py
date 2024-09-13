from collections import defaultdict
import logging
import forgi
import forgi.graph.bulge_graph as fgb
import networkx as nx

logger = logging.getLogger('rna_graphs')


# Generate tree with nucleotide numbers
def dotbracket_to_tree(db):
    bg = fgb.BulgeGraph.from_dotbracket(db)
    T = _generate_tree(bg, True)

    return T


# Generate Secondary Structure Tree
def _generate_tree(bg, include_nucleotides=False):
    """Generate networkx tree of RNA secondary structure from forgi BulgeGraph

    Args:
        bg (forgi.graph.bulge_graph): forgi BulgeGraph
        include_nucleotides (bool): Add nucleotide data to each node in a tree

    Returns:
        networkx.Graph: Secondary structure tree
    """
    G = nx.Graph()

    connect_to_virtual_node = False

    # 1. Add all hairpin and interior loops as nodes.
    for el in bg.defines.keys():
        if el[0] in ['h', 'i']:
            # For interior loops, stem_length returns length of the shorter strand
            if include_nucleotides:
                # nucs=unpaired nucs, defines=nucleotides defining this loop
                G.add_node(el, type=el[0], unpaired=bg.stem_length(
                    el), nucs=bg.elements_to_nucleotides([el]), defines=bg.define_a(el))
            else:
                G.add_node(el, type=el[0], unpaired=bg.stem_length(el))

    # 2. Add all stems, not connected to a multiloop, as edges.
    for el, adj in bg.edges.items():
        if el[0] == 's':
            if len(adj) == 2:
                a, b = adj
                if a[0] != 'm' and b[0] != 'm':
                    G.add_edge(a, b, length=bg.stem_length(el))
                elif a[0] == 'm' and b[0] == 'm':   # Special base: see below
                    G.add_node('o', type='n', unpaired=0)
                    connect_to_virtual_node = (a, bg.stem_length(el))
            # Special case: stem at the beginning/end of the sequence with no t0 or f0; add virtual node
            elif len(adj) == 1:
                G.add_node('o', type='n', unpaired=0)
                G.add_edge('o', list(adj)[0], length=bg.stem_length(el))

    curr_m_i = 0
    connected_ms = []
    mstars = {}                 # To which m* was m merged to
    mstars_ms = defaultdict(lambda: [])  # Which m-s contain a m*
    tf_set = set(['t0', 'f0'])

    # 2a. Find zero-length bulges
    zerolen = set([e for e in bg.defines if len(bg.defines[e]) == 0])
    logger.debug('Zerolen elements: %s', zerolen)

    # 3. Iterate over all true multiloops. (ms contains only m-nodes)
    for ms in bg.junctions:
        curr_m = f'm*{curr_m_i}'

        ms = set(ms)  # - tf_set

        # 4. Generate new multiloop node.
        if include_nucleotides:
            G.add_node(curr_m, type='m*',
                       unpaired=sum(bg.stem_length(m) for m in ms), nucs=bg.elements_to_nucleotides(ms), defines=bg.define_a(el))
        else:
            G.add_node(curr_m, type='m*',
                       unpaired=sum(bg.stem_length(m) for m in ms))

        # 5. Get all s-nodes in current multiloop.
        ss = []
        for m in ms:
            for adj in (bg.edges[m]):
                ss.append(adj)

            mstars[m] = curr_m
            mstars_ms[curr_m].append(m)

        # 6. Check if current multiloop is exterior.
        exterior_loop = 't0' in ms or 'f0' in ms or len(ms) == 1

        # 7. Connect current multiloop node to all adjacent elements, and save adjacent multiloop components.
        for s in set(ss):
            adj = bg.edges[s]
            adj_other = adj - ms        # Elements, adjacent to s, but not from current multiloop
            adj_other_types = [a[0]
                               for a in adj_other if a[0] in ['m', 't', 'f']]
            stem_length = bg.stem_length(s)

            if len(adj_other_types):     # Multiloop is connected to an adjacent multiloop
                other_ms = adj - ms  # - tf_set
                connected_ms.append([curr_m, stem_length, other_ms.pop()])
            else:                   # Multiloop is NOT connected to any adjacent multiloop
                for a in adj_other:
                    G.add_edge(curr_m, a, length=stem_length)

        curr_m_i += 1

    # 8. Remove both overhangs.
    if 't0' in G:
        G.remove_node('t0')
    if 'f0' in G:
        G.remove_node('f0')

    # 9. Connect adjacent multiloops.
    logger.debug('Connected ms: %s', connected_ms)
    logger.debug('mstars: %s', mstars)

    for curr_m, stem_length, other_m in connected_ms:
        logger.debug('Adding edge {curr_m}-{mstars[other_m]}')
        G.add_edge(curr_m, mstars[other_m], length=stem_length)

    if connect_to_virtual_node:
        G.add_edge('o', mstars[connect_to_virtual_node[0]],
                   length=connect_to_virtual_node[1])

    # 10. Remove zerolen bulges
    logger.debug('mstars %s', mstars)
    logger.debug('mstars_ms %s', mstars_ms)
    logger.debug('zerolen %s', zerolen)
    for z in zerolen:
        # It can happen that two (or more) zerolen bulges get merged into one m*
        zero_mstar = mstars[z]
        nonzero_ms = set(mstars_ms[zero_mstar]) - zerolen
        logger.debug('z: %s zero_mstar: %s G[zero_mstar]:  %s', z, zero_mstar, G[zero_mstar])

        if len(nonzero_ms) == 0:
            if len(G[zero_mstar]) == 2:     # Common case
                e1, e2 = G[zero_mstar]
                l = G.edges[zero_mstar, e1]['length'] + \
                    G.edges[zero_mstar, e2]['length']
                nx.contracted_nodes(G, e1, zero_mstar, self_loops=False, copy=False)                
                G.edges[e1, e2]['length'] = l

                logger.debug('Contracting %s -> %s / %s', e1, zero_mstar, e2)
            elif len(G[zero_mstar]) > 2:    # Rare case: 3(or more)-starred junction
                # Node should not be removed.
                logger.debug('Node not removed.')
            else:
                assert False, 'Error! len(G[zero_mstar]) < 2 ?!'

    for i, j, _ in G.edges(data=True):
        G[i][j]['1/length'] = 1/G[i][j]['length']

    # Validate generated tree
    try:
        nx.algorithms.cycles.find_cycle(G)
        assert False, 'Error! Cycle found!'
    except nx.NetworkXNoCycle:
        pass

    for node in G.nodes:
        if node[0] == 'm' and node[0:2] != 'm*':
            assert False, 'Error! m-node left in tree!'
            
    assert nx.number_connected_components(G) == 1, 'Generated tree is not 1-component!'

    tree_basepairs = sum(nx.get_edge_attributes(G, 'length').values())
    seq_basepairs = sum(bg.stem_length(s) for s in bg.stem_iterator())
    assert tree_basepairs == seq_basepairs, f'Error! Number of basepairs not equal! Tree: {tree_basepairs} vs. sequence: {seq_basepairs}'

    return G


def RAGify_tree(T):
    """Convert simple mapping tree to RAG tree.

    Args:
        T (nx.Graph): Simple mapping tree

    Returns:
        nx.Graph: RAG mapping tree
    """
    T_new = T.copy()
    merges = {}

    order = {'i':0, 'h':1, 'm*':2}

    for n1, n2, data in T.edges(data=True):
        if data['length'] == 1:
            if n1 not in T_new.nodes and n2 not in T_new.nodes:
                logger.debug('%s and %s already gone', n1, n2)
                continue

            if n1 not in T_new.nodes:
                n1 = merges[n1]

                while n1 not in T_new.nodes:
                    n1 = merges[n1]

            if n2 not in T_new.nodes:
                n2 = merges[n2]

                while n2 not in T_new.nodes:
                    n2 = merges[n2]

            t1, t2 = T_new.nodes[n1]['type'], T_new.nodes[n2]['type']
            o1, o2 = order[t1], order[t2]

            # Check priority
            if o1 > o2:
                n1, n2 = n2, n1
            
            logger.debug('Merging %s %s', n1, n2)
            T_new = _merge_nodes(T_new, n1, n2)
            merges[n1] = n2

    T_copy = T_new.copy()
    for node in T_copy.nodes:
        if T.nodes[node]['type'] != 'i':
            continue

        if T.degree[node] != 2:
            logger.error('Internal loop with degree %s?!', T.degree[node])

        nucs = T_new.nodes[node]['nucs']

        if (len(nucs) == 1) or (len(nucs) == 2 and nucs[1] - nucs[0] > 1):
            logger.debug('Deleting %s - %s neighbors %s %s', node, nucs, list(T_new.neighbors(node)), T_new.nodes[node])
            T_new = _delete_node(T_new, node)

    return T_new

def _merge_nodes(T: nx.Graph, node1, node2):
    """Merge two nodes into one node.

    Args:
        T (nx.Graph): Tree
        node1: First node
        node2: Second node

    Returns:
        nx.Graph: Tree with merged nodes
    """
    if node1 not in T.nodes:
        logger.debug(node1, 'already gone')
        return False

    if T.degree[node1] > 2 and T.degree[node1] > 2:
        # Merging two multiloops
        # TODO: transfer unpaired nucleotides?
        T = nx.contracted_nodes(T, node1, node2, False)

    else: 
        if T.degree[node1] > 1:
            node1_neighbors = set(T.neighbors(node1))
            connect_to = (node1_neighbors - set([node2])).pop()

            edge_remove, other_edge = T.edges[(node1, node2)], T.edges[(node1, connect_to)]

            total_length = edge_remove['length'] + other_edge['length']

            T.add_edge(connect_to, node2, length=total_length)
            T[connect_to][node2]['1/length'] = 1/total_length

            new_unpaired_nucs = T.nodes[node1]['nucs'] + \
                list(set(T.nodes[node1]['defines']) & set(T.nodes[node2]['defines']))
        else:
            new_unpaired_nucs = T.nodes[node1]['nucs'] + T.nodes[node1]['defines']

        T.nodes[node2]['nucs'] += new_unpaired_nucs

        logger.debug('Adding %s to %s', new_unpaired_nucs, T.nodes[node2]['nucs'])

        T.remove_node(node1)

    return T

def _delete_node(T: nx.Graph, node):
    """Delete degree 2 node and merge adjacent edges.

    Args:
        T (nx.Graph): Tree  
        node: Node

    Returns:
        nx.Graph: Tree without deleted node
    """
    n1, n2 = T.neighbors(node)
    e1, e2 = T.edges[(node, n1)], T.edges[(node, n2)]

    total_length = e1['length'] + e2['length']

    T.add_edge(n1, n2, length=total_length)
    T[n1][n2]['1/length'] = 1/total_length

    T.remove_node(node)

    return T
