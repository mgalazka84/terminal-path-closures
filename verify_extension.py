#!/usr/bin/env python3
"""Exact small-graph and constructed-family checks for the second extension.

This is falsification evidence, not a proof or a novelty certificate.
The independent minor oracle is imported from the first verifier.
"""
from collections import deque
from itertools import combinations
from math import prod
from pathlib import Path
import json
import random
import time

from verify_boundary_bounds import (
    bits, from_edges, graph_data, largest_k2_minor, partitions,
)


def closure(adj, terminals, h):
    result = set(terminals)
    for root in terminals:
        stack = [(root, (root,))]
        while stack:
            u, path = stack.pop()
            if len(path) - 1 == h:
                continue
            for v in bits(adj[u]):
                if v in path:
                    continue
                newpath = path + (v,)
                if v in terminals:
                    result.update(newpath)
                else:
                    stack.append((v, newpath))
    return result


def distances(adj, root, allowed=None):
    allowed = set(range(len(adj))) if allowed is None else set(allowed)
    d = {root: 0}
    queue = deque([root])
    while queue:
        u = queue.popleft()
        for v in bits(adj[u]):
            if v in allowed and v not in d:
                d[v] = d[u] + 1
                queue.append(v)
    return d


def small_graph_checks():
    counts = dict(graphs=0, path_closure_inequalities=0,
                  partition_boundary_inequalities=0)
    for n in range(1, 6):
        possible = list(combinations(range(n), 2))
        for mask in range(1 << len(possible)):
            adj = from_edges(n, [e for i, e in enumerate(possible)
                                 if mask & (1 << i)])
            neigh, connected, radii = graph_data(adj)
            t = max(2, largest_k2_minor(adj, neigh, connected) + 1)
            for terminal_mask in range(1, 1 << n):
                q = set(bits(terminal_mask))
                for h in range(1, 5):
                    a = prod(1 + i*(t-1) for i in range(1, h))
                    actual = len(closure(adj, q, h))
                    assert actual <= 1+a*(len(q)-1), (adj, q, h, t, actual, a)
                    counts['path_closure_inequalities'] += 1
            for part in partitions(n):
                if any(not connected[b] for b in part):
                    continue
                boundary = set()
                active = set()
                for i, j in combinations(range(len(part)), 2):
                    if neigh[part[i]] & part[j]:
                        active.update((i, j))
                        boundary.update(bits(part[i] & neigh[part[j]]))
                        boundary.update(bits(part[j] & neigh[part[i]]))
                k = max((radii[part[i]] for i in active), default=0)
                bound = len(active) + 2*k*(t-1)*max(0, len(active)-1)
                assert len(boundary) <= bound
                counts['partition_boundary_inequalities'] += 1
            counts['graphs'] += 1
    return counts


def outerplanar_graph(n, rng):
    # An edge is active exactly when it is on the outer boundary.
    edges = [(0, 1)]
    active = [(0, 1)]
    for v in range(2, n):
        u, w = active.pop(rng.randrange(len(active)))
        edges.extend(((u, v), (v, w)))
        active.extend(((u, v), (v, w)))
    return from_edges(n, edges)


def check_partition(adj, centres, groups):
    ds = [distances(adj, c) for c in centres]
    owner = [min(range(len(centres)), key=lambda i: (ds[i][v], i))
             for v in range(len(adj))]
    cells = [{v for v in range(len(adj)) if owner[v] == i}
             for i in range(len(centres))]
    radius = []
    for i, c in enumerate(centres):
        local = distances(adj, c, cells[i])
        assert set(local) == cells[i]
        radius.append(max(local.values()))
    boundary = set()
    active = set()
    for u in range(len(adj)):
        for v in bits(adj[u]):
            i, j = owner[u], owner[v]
            if groups[i] != groups[j]:
                boundary.update((u, v))
                active.update((i, j))
    k = max((radius[i] for i in active), default=0)
    bound = len(active) + 4*k*max(0, len(active)-1)  # t=3, outerplanar
    assert len(boundary) <= bound, (len(adj), centres, k, len(boundary), bound)
    return len(boundary), bound, len(active), max(radius), owner


def random_outerplanar_checks():
    rng = random.Random(2209202602)
    cases = nonvacuous = 0
    for n in (20, 50, 100, 200):
        for _ in range(15):
            adj = outerplanar_graph(n, rng)
            for m in (2, 3, 5, 10):
                centres = rng.sample(range(n), m)
                for groups in (list(range(m)), [rng.randrange(2) for _ in range(m)]):
                    b, bound, _, _, _ = check_partition(adj, centres, groups)
                    cases += 1
                    nonvacuous += int(0 < bound < n)
    return {'cases': cases, 'cases_with_positive_bound_below_order': nonvacuous,
            'seed': 2209202602, 'orders': [20, 50, 100, 200]}


def ladder_chain(k, m):
    # Centres 0,...,m-1. Consecutive ladders share just their common centre.
    edges = []
    n = m
    expected = list(range(m))
    for block in range(m-1):
        rails = []
        for centre in (block, block+1):
            rail = []
            for pos in range(-k, k+1):
                if pos == 0:
                    rail.append(centre)
                else:
                    rail.append(n)
                    expected.append(centre)
                    n += 1
            rails.append(rail)
            edges.extend(zip(rail, rail[1:]))
        edges.extend(zip(*rails))
    return from_edges(n, edges), expected


def sharpness_checks():
    out = []
    for k in (1, 2, 5):
        for m in (2, 3, 10, 20):
            adj, expected = ladder_chain(k, m)
            b, bound, active, radius, owner = check_partition(
                adj, list(range(m)), list(range(m)))
            assert owner == expected
            assert active == m and radius == k
            assert b == len(adj) == (4*k+1)*m - 4*k
            assert bound == b == (4*k+1)*m - 4*k
            out.append({'k': k, 'cells': m, 'boundary': b, 'upper_bound': bound})
    return out


def binary_triangulation(depth):
    edges = [(0, 1)]
    frontier = [(0, 1)]
    n = 2
    for _ in range(depth):
        nxt = []
        for u, v in frontier:
            w = n
            n += 1
            edges.extend(((u, w), (w, v)))
            nxt.extend(((u, w), (w, v)))
        frontier = nxt
    assert n == 2**depth + 1
    return from_edges(n, edges)


def iterated_closure_checks():
    out = []
    for h, j in ((2, 1), (2, 5), (3, 1), (3, 2), (3, 3), (4, 2)):
        depth = j*(h-1)
        adj = binary_triangulation(depth)
        q = {0, 1}
        sizes = [len(q)]
        for _ in range(j):
            q = closure(adj, q, h)
            sizes.append(len(q))
        assert len(q) == len(adj) == 2**depth + 1
        out.append({'h': h, 'iterations': j, 'depth': depth,
                    'vertices': len(adj), 'closure_sizes': sizes})
    return out


def main():
    start = time.monotonic()
    result = {'status': 'all assertions passed',
              'small_graphs': small_graph_checks(),
              'outerplanar_boundary_checks': random_outerplanar_checks(),
              'sharp_ladder_chains': sharpness_checks(),
              'iterated_closure_examples': iterated_closure_checks(),
              'limitations': ['Finite checks supplement the mathematical proofs.',
                              'No test of the full distributed algorithm.',
                              'No assertion of literature priority.']}
    result['seconds'] = round(time.monotonic() - start, 3)
    Path(__file__).with_name('extension_verification_results.json').write_text(
        json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
