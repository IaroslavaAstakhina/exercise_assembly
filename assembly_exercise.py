from collections import Counter

# Functions carried over from Exercise 04

def build_adjacency_list(edges, directed=False):
    graph = {}
    for u, v in edges:
        graph.setdefault(u, [])
        graph.setdefault(v, [])
        graph[u].append(v)
        if not directed:
            graph[v].append(u)
    return graph


def degrees(graph):
    return {v: len(neighbours) for v, neighbours in graph.items()}


def in_out_degree(graph):
    out_degree = {v: len(neighbours) for v, neighbours in graph.items()}
    in_degree = {v: 0 for v in graph}
    for neighbours in graph.values():
        for w in neighbours:
            in_degree[w] = in_degree.get(w, 0) + 1
    for v in in_degree:
        out_degree.setdefault(v, 0)
    return in_degree, out_degree


def dfs(graph, start):
    """Vertices reachable from start, in order of first visit.
    Unvisited neighbours are visited in sorted (alphabetical) order."""
    visited = set()
    order = []
    stack = [start]
    while stack:
        v = stack.pop()
        if v in visited:
            continue
        visited.add(v)
        order.append(v)
        for w in sorted(graph.get(v, []), reverse=True):
            if w not in visited:
                stack.append(w)
    return order


def connected_components(graph):
    components = []
    seen = set()
    for v in sorted(graph):
        if v not in seen:
            comp = set(dfs(graph, v))
            components.append(comp)
            seen |= comp
    return components


def is_weakly_connected(graph):
    shadow = {v: [] for v in graph}
    for u, neighbours in graph.items():
        for v in neighbours:
            shadow.setdefault(v, [])
            shadow[u].append(v)
            shadow[v].append(u)
    return len(connected_components(shadow)) == 1


def eulerian_type(graph, directed=False):
    """Returns 'circuit', 'path' or 'none'."""
    if not directed:
        deg = degrees(graph)
        with_edges = {v: nb for v, nb in graph.items() if nb}
        if len(connected_components(with_edges)) > 1:
            return "none"
        odd = sum(1 for d_ in deg.values() if d_ % 2 == 1)
        if odd == 0:
            return "circuit"
        if odd == 2:
            return "path"
        return "none"

    in_deg, out_deg = in_out_degree(graph)
    with_edges = {v: nb for v, nb in graph.items() if nb or in_deg[v] > 0}
    if with_edges and not is_weakly_connected(with_edges):
        return "none"
    excess = [out_deg[v] - in_deg[v] for v in graph]
    if any(e not in (-1, 0, 1) for e in excess):
        return "none"
    plus = excess.count(1)
    minus = excess.count(-1)
    if plus == 0 and minus == 0:
        return "circuit"
    if plus == 1 and minus == 1:
        return "path"
    return "none"


# Task 2: De Bruijn graph

def kmer_composition(text, k):
    return [text[i:i+k] for i in range(len(text) - k + 1)]


def debruijn_graph(kmers):
    graph = {}
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        graph.setdefault(prefix, [])
        graph.setdefault(suffix, [])
        graph[prefix].append(suffix)
    return graph


def spell_path(nodes):
    genome = nodes[0]
    for node in nodes[1:]:
        genome += node[-1]
    return genome


def test_task2():
    assert kmer_composition('GATTACA', 4) == ['GATT', 'ATTA', 'TTAC', 'TACA']

    toy_a = ['ACT', 'AGA', 'CAG', 'CCA', 'CTG', 'GAA',
             'GAC', 'GGT', 'GTC', 'TCC', 'TGA']
    graph = debruijn_graph(toy_a)
    expected_nodes = {km[:-1] for km in toy_a} | {km[1:] for km in toy_a}
    assert set(graph) == expected_nodes
    assert len(graph) == 11
    assert sum(len(neighbors) for neighbors in graph.values()) == len(toy_a)
    assert Counter((u, v) for u, ns in graph.items() for v in ns) == Counter(
        (km[:-1], km[1:]) for km in toy_a)

    incoming, outgoing = in_out_degree(graph)
    assert sum(incoming.values()) == sum(outgoing.values()) == len(toy_a)
    assert (incoming['GG'], outgoing['GG']) == (0, 1)
    assert (incoming['AA'], outgoing['AA']) == (1, 0)
    assert (incoming['GA'], outgoing['GA']) == (2, 2)
    assert is_weakly_connected(graph)
    assert len(connected_components({u: list(set(ns)) for u, ns in _shadow(graph).items()})) == 1
    assert eulerian_type(graph, directed=True) == 'path'
    assert spell_path(['GAT', 'ATT', 'TTA']) == 'GATTA'

    for text in ['GATTACA', 'CATCATCAT', 'AAAAAA', 'ATAT', 'ACGT']:
        for k in range(1, len(text) + 1):
            kmers = kmer_composition(text, k)
            assert len(kmers) == len(text) - k + 1
            assert spell_path(kmers) == text
            if k >= 2:
                g = debruijn_graph(kmers)
                assert sum(len(ns) for ns in g.values()) == len(kmers)
                assert set(g) == ({km[:-1] for km in kmers} | {km[1:] for km in kmers})
                assert eulerian_type(g, directed=True) in ('path', 'circuit')

    repeats = debruijn_graph(kmer_composition('CATCATCAT', 3))
    assert repeats['CA'].count('AT') == 3
    assert repeats['AT'].count('TC') == 2
    assert repeats['TC'].count('CA') == 2
    assert eulerian_type(repeats, directed=True) == 'path'
    assert eulerian_type(debruijn_graph(kmer_composition('ATATA', 2)), directed=True) == 'circuit'

    print('Toy A adjacency list:')
    for node in sorted(graph):
        print(f'{node} -> {graph[node]}')
    print('\nNode | in | out')
    for node in sorted(graph):
        print(f'{node:>4} | {incoming[node]:>2} | {outgoing[node]:>3}')
    print('\nToy A: 11 nodes, 11 edges, Eulerian path')
    print('Repeated k-mers and string spelling: passed')
    print('All Task 2 tests passed!')


def _shadow(graph):
    shadow = {v: [] for v in graph}
    for u, ns in graph.items():
        for v in ns:
            shadow[u].append(v)
            shadow[v].append(u)
    return shadow




# Task 3: Hierholzer's algorithm (walks and splicing)

def hierholzer(graph, start, verbose=False):
    # Work on a copy: keep the original adjacency lists unchanged.
    unused = {node: sorted(neighbors) for node, neighbors in graph.items()}

    def make_walk(node):
        walk = [node]
        while unused[node]:
            next_node = unused[node].pop(0)  # Alphabetically first unused edge
            walk.append(next_node)
            node = next_node
        return walk

    result = make_walk(start)
    if verbose:
        print('Walk 1:', ' -> '.join(result))

    walk_number = 1
    while any(unused.values()):
        # Find the first POSITION in the current spliced walk with unused edges.
        for position, node in enumerate(result):
            if unused[node]:
                break
        else:
            raise ValueError('Unused edges are not reachable from the walk')

        new_walk = make_walk(node)
        walk_number += 1
        if verbose:
            print(f'Walk {walk_number} (restart at {node}):',
                  ' -> '.join(new_walk))

        # Replace only the selected occurrence of the restart node.
        result = result[:position] + new_walk + result[position + 1:]
        if verbose:
            print('After splice:', ' -> '.join(result))

    return result


def eulerian_path(graph, verbose=False):
    kind = eulerian_type(graph, directed=True)
    if kind == 'none':
        raise ValueError('Graph has no Eulerian path')

    if kind == 'path':
        incoming, outgoing = in_out_degree(graph)
        start = next(node for node in graph
                     if outgoing[node] - incoming[node] == 1)
    else:
        start = min(graph)

    return hierholzer(graph, start, verbose)


def check_eulerian_walk(graph, kmers):
    original = {node: neighbors.copy() for node, neighbors in graph.items()}
    walk = eulerian_path(graph)
    edge_count = sum(len(neighbors) for neighbors in graph.values())

    # 1. Exactly E + 1 nodes.
    assert len(walk) == edge_count + 1

    # 2. Every edge, including parallel edges, used exactly once.
    expected_edges = Counter((u, v) for u, ns in graph.items() for v in ns)
    actual_edges = Counter(zip(walk, walk[1:]))
    assert actual_edges == expected_edges

    # 3. Reconstructed string contains exactly the input k-mers.
    assert Counter(kmer_composition(spell_path(walk), len(kmers[0]))) == Counter(kmers)

    # 4. Original graph is unchanged.
    assert graph == original

    # 6. Reversing input k-mers must not change the walk.
    reversed_graph = debruijn_graph(kmers[::-1])
    assert eulerian_path(reversed_graph) == walk
    return walk


def test_task3():
    toy_a = ['ACT', 'AGA', 'CAG', 'CCA', 'CTG', 'GAA',
             'GAC', 'GGT', 'GTC', 'TCC', 'TGA']
    toy_b = ['ACT', 'AGC', 'ATT', 'CCG', 'CGA', 'CGT', 'CTG',
             'GAC', 'GAG', 'GAT', 'GTC', 'TCC', 'TCG', 'TGA', 'TTC']

    walk_a = check_eulerian_walk(debruijn_graph(toy_a), toy_a)
    assert spell_path(walk_a) == 'GGTCCAGACTGAA'

    graph_b = debruijn_graph(toy_b)
    print('_________________TASK 3_____________________')
    print('\nToy B - Hierholzer trace:')
    walk_b = eulerian_path(graph_b, verbose=True)
    expected_b = ['GA', 'AT', 'TT', 'TC', 'CG', 'GT', 'TC', 'CC',
                  'CG', 'GA', 'AC', 'CT', 'TG', 'GA', 'AG', 'GC']
    assert walk_b == expected_b
    assert check_eulerian_walk(graph_b, toy_b) == expected_b
    assert spell_path(walk_b) == 'GATTCGTCCGACTGAGC'

    # Extra examples: parallel edges, cycles, and repeated k-mers.
    for text, k in [('CATCATCAT', 3), ('ATATA', 2),
                    ('AAAAAA', 3), ('GATTACA', 4)]:
        kmers = kmer_composition(text, k)
        check_eulerian_walk(debruijn_graph(kmers), kmers)

    print('Toy B genome:', spell_path(walk_b))
    print('Task 3 tests passed!')


if __name__ == '__main__':
    test_task2()
    test_task3()
