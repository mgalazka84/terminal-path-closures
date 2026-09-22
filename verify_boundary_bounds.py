#!/usr/bin/env python3
"""Finite falsification checks for the proposed K_{2,t} boundary lemma.

The minor oracle enumerates connected set partitions and inspects their
quotient graphs. It is independent of the BFS-leaf proof of the lemma.
Uses only the Python standard library. This is evidence, not a proof.
"""
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import json
import random
import time


def bits(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


@lru_cache(None)
def partitions(n):
    if n == 0:
        return ((),)
    bit = 1 << (n - 1)
    ans = []
    for part in partitions(n - 1):
        ans.append(part + (bit,))
        for i in range(len(part)):
            ans.append(part[:i] + (part[i] | bit,) + part[i + 1:])
    return tuple(ans)


def graph_data(adj):
    n = len(adj)
    neighbors = [0] * (1 << n)
    connected = [False] * (1 << n)
    radii = [None] * (1 << n)
    for mask in range(1, 1 << n):
        bit = mask & -mask
        neighbors[mask] = neighbors[mask ^ bit] | adj[bit.bit_length() - 1]
        seen = bit
        while True:
            nxt = seen | (neighbors[seen] & mask)
            if nxt == seen:
                break
            seen = nxt
        if seen != mask:
            continue
        connected[mask] = True
        best = n
        for v in bits(mask):
            seen, eccentricity = 1 << v, 0
            while seen != mask:
                seen |= neighbors[seen] & mask
                eccentricity += 1
            best = min(best, eccentricity)
        radii[mask] = best
    return neighbors, connected, radii


def largest_k2_minor(adj, neighbors, connected):
    """Largest s with K_{2,s} a minor; 0 if none with s >= 1.

    Every minor model occurs in one enumerated partition: use its branch
    sets as blocks and all unused vertices as singleton blocks. Conversely,
    every quotient is a minor, and s common neighbors certify K_{2,s}.
    """
    best = 0
    for part in partitions(len(adj)):
        if len(part) <= best + 2 or any(not connected[b] for b in part):
            continue
        quotient = [0] * len(part)
        for i, j in combinations(range(len(part)), 2):
            if neighbors[part[i]] & part[j]:
                quotient[i] |= 1 << j
                quotient[j] |= 1 << i
        for i, j in combinations(range(len(part)), 2):
            best = max(best, (quotient[i] & quotient[j]).bit_count())
    return best


def check_graph(adj, counts):
    neighbors, connected, radii = graph_data(adj)
    s = largest_k2_minor(adj, neighbors, connected)
    t = max(2, s + 1)  # strongest applicable parameter; bounds grow in t
    subsets = [a for a in range(1, 1 << len(adj)) if connected[a]]
    for a, b in combinations(subsets, 2):
        if a & b or not neighbors[a] & b:
            continue
        wa = (a & neighbors[b]).bit_count()
        wb = (b & neighbors[a]).bit_count()
        m = sum((adj[v] & b).bit_count() for v in bits(a))
        assert wa <= 1 + radii[a] * (t - 1), (adj, a, b, t, wa)
        assert wb <= 1 + radii[b] * (t - 1), (adj, a, b, t, wb)
        assert 2 * m <= t + 1 + (radii[a] + radii[b]) * (t - 1)**2, (
            adj, a, b, t, m, radii[a], radii[b])
        counts['connected_set_pairs'] += 1
        counts['scalar_inequalities'] += 3
    counts['graphs'] += 1
    counts['minor_free_parameter_histogram'][str(t)] = (
        counts['minor_free_parameter_histogram'].get(str(t), 0) + 1)


def from_edges(n, edges):
    adj = [0] * n
    for u, v in edges:
        adj[u] |= 1 << v
        adj[v] |= 1 << u
    return adj


def oracle_sanity():
    cases = [
        ('path', from_edges(7, [(i, i+1) for i in range(6)]), 1),
        ('triangle', from_edges(3, combinations(range(3), 2)), 1),
        ('cycle', from_edges(6, [(i, (i+1) % 6) for i in range(6)]), 2),
        ('clique', from_edges(7, combinations(range(7), 2)), 5),
        ('K2,4', from_edges(6, [(i, j) for i in range(2) for j in range(2, 6)]), 4),
    ]
    for name, adj, expected in cases:
        neigh, conn, _ = graph_data(adj)
        actual = largest_k2_minor(adj, neigh, conn)
        assert actual == expected, (name, actual, expected)
    return [name for name, _, _ in cases]


def ladder_sanity():
    cases = []
    for k in (1, 2, 3, 5, 10):
        m = 2*k+1
        adj = from_edges(2*m, [(i, i+1) for i in range(m-1)]
                         + [(m+i, m+i+1) for i in range(m-1)]
                         + [(i, m+i) for i in range(m)])
        def distances(root):
            dist = {root: 0}
            queue = [root]
            for u in queue:
                for v in bits(adj[u]):
                    if v not in dist:
                        dist[v] = dist[u] + 1
                        queue.append(v)
            return dist
        da, db = distances(k), distances(m+k)
        assert {v for v in range(2*m) if da[v] < db[v]} == set(range(m))
        assert {v for v in range(2*m) if db[v] < da[v]} == set(range(m, 2*m))
        assert max(min(da[v], db[v]) for v in range(2*m)) == k
        assert m == 1 + k*(3-1)
        cases.append({'k': k, 'vertices': 2*m, 'boundary_per_cell': m})
    return cases


def main():
    start = time.monotonic()
    counts = {'graphs': 0, 'connected_set_pairs': 0, 'scalar_inequalities': 0,
              'minor_free_parameter_histogram': {}}
    oracle = oracle_sanity()
    for n in range(1, 6):
        possible = list(combinations(range(n), 2))
        for mask in range(1 << len(possible)):
            adj = from_edges(n, [e for i, e in enumerate(possible) if mask & (1 << i)])
            check_graph(adj, counts)
    exhaustive = counts['graphs']
    rng = random.Random(22092026)
    for n in (6, 7, 8):
        for _ in range(100):
            p = rng.choice((0.15, 0.3, 0.5, 0.7, 0.9))
            adj = from_edges(n, [e for e in combinations(range(n), 2) if rng.random() < p])
            check_graph(adj, counts)
    result = {
        'status': 'all assertions passed',
        'scope': 'All labelled simple graphs on 1..5 vertices; 100 seeded random graphs for each of n=6,7,8.',
        'exhaustive_graphs': exhaustive,
        'random_graphs': 300,
        'seed': 22092026,
        'oracle_sanity_cases': oracle,
        'ladder_cases': ladder_sanity(),
        **counts,
        'seconds': round(time.monotonic() - start, 3),
        'limitation': 'Finite falsification checks; no claim of proof or exhaustive literature novelty.',
    }
    path = Path(__file__).with_name('verification_results.json')
    path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
