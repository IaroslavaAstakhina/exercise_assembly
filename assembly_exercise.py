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
    print('All Task 3 tests passed!')


# Task 4: Repeats and contigs
#A node is 1-in-1-out if its in-degree and its out-degree are both 1. A path is non-branching if every node inside it (not its first or last node) is 1-in-1-out.
#It is maximal if it cannot be extended at either end without becoming branching.
def contigs(graph):
    incoming, outgoing = in_out_degree(graph)
    unused = {node: sorted(neighbors) for node, neighbors in graph.items()}
    result = []

    def one_in_one_out(node):
        return incoming[node] == 1 and outgoing[node] == 1

    # Start a path at each node that is not 1-in-1-out.
    for start in sorted(graph):
        if one_in_one_out(start):
            continue
        while unused[start]:
            node = unused[start].pop(0)
            path = [start, node]
            while one_in_one_out(node):
                next_node = unused[node].pop(0)
                path.append(next_node)
                node = next_node
            result.append(spell_path(path))


    for start in sorted(graph):
        while unused[start]:
            path = [start]
            node = start
            while True:
                node = unused[node].pop(0)
                path.append(node)
                if node == start:
                    break
            result.append(spell_path(path))

    return sorted(result)


def check_contigs(graph, kmers):
    original = {node: neighbors.copy() for node, neighbors in graph.items()}
    pieces = contigs(graph)
    k = len(kmers[0])

    # Every edge must appear in exactly one contig (including repeats).
    expected = Counter((u, v) for u, ns in graph.items() for v in ns)
    actual = Counter()
    for piece in pieces:
        piece_graph = debruijn_graph(kmer_composition(piece, k))
        actual.update((u, v) for u, ns in piece_graph.items() for v in ns)
    assert actual == expected
    assert sum(len(piece) - k + 1 for piece in pieces) == len(kmers)
    assert graph == original
    assert pieces == contigs(debruijn_graph(kmers[::-1]))
    return pieces


def test_task4():
    print('\n_________________TASK 4_____________________')
    toy_c = ['ATC', 'CCT', 'CTT', 'GGT', 'GTA', 'GTG',
             'TAT', 'TCT', 'TGG', 'TGT', 'TTG']
    graph_c = debruijn_graph(toy_c)

    print('\nToy C - Hierholzer trace:')
    walk_c = eulerian_path(graph_c, verbose=True)
    genome = spell_path(walk_c)
    assert genome == 'CCTTGTGGTATCT'
    assert check_eulerian_walk(graph_c, toy_c) == walk_c
    print('Toy C genome:', genome)

    print('\nToy C contigs:')
    pieces = check_contigs(graph_c, toy_c)
    expected = ['CCT', 'CTTG', 'GTATCT', 'GTG', 'TGGT', 'TGT']
    assert pieces == expected
    for piece in pieces:
        print(' ', piece)
    print('Toy C:', len(pieces), 'contigs,', len(toy_c), 'edges covered')

    # A string with no repeated (k-1)-mer: one contig, the string itself.
    text = 'GATTACA'
    kmers = kmer_composition(text, 3)
    assert check_contigs(debruijn_graph(kmers), kmers) == [text]

    # A single 1-in-1-out cycle: one contig, starting at min node.
    text = 'ACGTTGCAAC'
    kmers = kmer_composition(text, 3)
    cycle_pieces = check_contigs(debruijn_graph(kmers), kmers)
    assert len(cycle_pieces) == 1
    assert len(cycle_pieces[0]) == len(text)
    assert cycle_pieces[0][:2] == min(debruijn_graph(kmers))


    kmers = ['AAB', 'AAB']
    assert check_contigs(debruijn_graph(kmers), kmers) == ['AAB', 'AAB']

    kmers = ['ABC', 'BCD', 'XYX', 'YXY']
    assert check_contigs(debruijn_graph(kmers), kmers) == ['ABCD', 'XYXY']

    print('Task 4 tests passed!')


if __name__ == '__main__':
    test_task2()
    test_task3()
    test_task4()

# Task 5: A mystery genome

def read_fasta_kmers(filename, k):
    """Collect distinct k-mers from all FASTA reads."""
    from Bio import SeqIO

    kmers = set()
    read_count = 0
    for record in SeqIO.parse(filename, 'fasta'):
        read_count += 1
        sequence = str(record.seq).upper()
        kmers.update(kmer_composition(sequence, k))
    if read_count == 0 or not kmers:
        raise ValueError('No reads or no k-mers found in ' + filename)
    return sorted(kmers)


def trim_cycle(sequence, k):

    overlap = k - 1
    if sequence[:overlap] == sequence[-overlap:]:
        return sequence[:-overlap]
    return sequence


def same_circle(s, t):
    return len(s) == len(t) and s in (t + t)


def check_task5_contigs(graph, pieces, k):

    expected = Counter((u, v) for u, ns in graph.items() for v in ns)
    actual = Counter()
    for piece in pieces:
        for kmer in kmer_composition(piece, k):
            actual[(kmer[:-1], kmer[1:])] += 1
    assert actual == expected


def task5(filename='reads.fasta'):
    from pathlib import Path

    if not Path(filename).is_file():
        raise FileNotFoundError(
            f'Cannot find {filename}. Put reads.fasta next to assembly_exercise.py '
            'or pass the full path to task5().'
        )

    print('\n_________________TASK 5_____________________')
    results = []

    for k in [9, 12, 13, 14, 21]:
        kmers = read_fasta_kmers(filename, k)
        kmer_set = set(kmers)
        graph = debruijn_graph(kmers)
        node_count = len(graph)
        edge_count = sum(len(ns) for ns in graph.values())
        kind = eulerian_type(graph, directed=True)
        assert edge_count == len(kmer_set)

        pieces = contigs(graph)
        check_task5_contigs(graph, pieces, k)
        trimmed_pieces = [trim_cycle(piece, k) for piece in pieces]
        longest = max(map(len, trimmed_pieces), default=0)

        spelled = None
        if kind != 'none':
            original = {u: ns.copy() for u, ns in graph.items()}
            walk = eulerian_path(graph)
            full_string = spell_path(walk)

            # Self-checks BEFORE trimming (as requested by the README).
            assert len(walk) == edge_count + 1
            assert len(full_string) == edge_count + k - 1
            assert set(kmer_composition(full_string, k)) == kmer_set
            assert Counter(zip(walk, walk[1:])) == Counter(
                (u, v) for u, ns in graph.items() for v in ns
            )
            assert graph == original

            if kind == 'circuit':
                assert walk[0] == walk[-1]
                spelled = trim_cycle(full_string, k)
            else:
                spelled = full_string

        results.append({
            'k': k, 'V': node_count, 'E': edge_count,
            'type': kind, 'contigs': len(pieces),
            'longest': longest, 'spelled': spelled,
            'pieces': trimmed_pieces,
        })

    print('k   |V|    |E|    type      contigs  longest  spelled length')
    print('-' * 68)
    for r in results:
        spelled_len = str(len(r['spelled'])) if r['spelled'] is not None else '-'
        print(f"{r['k']:<3} {r['V']:<6} {r['E']:<6} {r['type']:<9} "
              f"{r['contigs']:<8} {r['longest']:<8} {spelled_len}")

    print('\nCircular sequence comparison:')
    spelled_results = [r for r in results if r['spelled'] is not None]
    for i, first in enumerate(spelled_results):
        for second in spelled_results[i + 1:]:
            equal = same_circle(first['spelled'], second['spelled'])
            print(f"k={first['k']} vs k={second['k']}: "
                  f"{'same circle' if equal else 'DIFFERENT circles'}")


    first_single = next((r for r in results if r['contigs'] == 1), None)
    blast_queries = []
    if first_single is not None:
        blast_queries.append((
            f"single_contig_k{first_single['k']}", first_single['pieces'][0]
        ))
        for r in results:
            if r['k'] < first_single['k'] and r['spelled'] is not None:
                blast_queries.append((f"spelled_k{r['k']}", r['spelled']))

    if blast_queries:
        output = Path(filename).with_name('blast_queries.fasta')
        with output.open('w') as f:
            for name, sequence in blast_queries:
                f.write(f'>{name}\n')
                for i in range(0, len(sequence), 80):
                    f.write(sequence[i:i+80] + '\n')
        print(f'\nBLAST queries saved to: {output}')
    else:
        print('\nNo single contig found at the tested k values.')

    print('Task 5 checks passed!')
    return results


if __name__ == '__main__':
    from pathlib import Path
    fasta = Path(__file__).with_name('reads.fasta')
    if fasta.exists():
        task5(str(fasta))
    else:
        print('\nTask 5 not run: reads.fasta is not next to this script.')
