#!/usr/bin/env python3
"""Targeted falsification checks for the random-colouring closure proof.

Tests only the new argument: retained-white contractions, exact witness
probabilities, and nonvacuous closure bounds on explicit outerplanar graphs.
The independent exact minor oracle is inherited from the first verifier.
"""
from functools import lru_cache
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
import json
import random
import time

from verify_boundary_bounds import bits, from_edges, graph_data, largest_k2_minor
from verify_extension import binary_triangulation, closure

RED, BLUE, WHITE = 0, 1, 2


@lru_cache(maxsize=None)
def max_k2_minor(adj):
    neigh, connected, _ = graph_data(adj)
    return largest_k2_minor(adj, neigh, connected)


def retained_minor(adj, terminals, colours):
    """Build the minor directly from monochromatic connected branch sets."""
    n = len(adj)
    seen, components = set(), []
    for root in range(n):
        if root in seen or colours[root] == WHITE:
            continue
        component = {root}
        seen.add(root)
        stack = [root]
        while stack:
            u = stack.pop()
            for v in bits(adj[u]):
                if v not in seen and colours[v] == colours[root]:
                    seen.add(v)
                    component.add(v)
                    stack.append(v)
        if component & terminals:
            components.append((component, colours[root]))
    owner = {v: i for i, (component, _) in enumerate(components)
             for v in component}
    retained = []
    for v in range(n):
        if colours[v] != WHITE:
            continue
        neighbour_colours = {components[owner[u]][1] for u in bits(adj[v])
                             if u in owner}
        if neighbour_colours == {RED, BLUE}:
            retained.append(v)
    y = len(components)
    owner.update({v: y+i for i, v in enumerate(retained)})
    quotient_edges = set()
    for u in range(n):
        if u not in owner:
            continue
        for v in bits(adj[u]):
            if v not in owner:
                continue
            a, b = owner[u], owner[v]
            if a != b and not (a >= y and b >= y):
                quotient_edges.add(tuple(sorted((a, b))))
    quotient = tuple(from_edges(y+len(retained), quotient_edges))
    return retained, y, quotient


def check_outcome(adj, terminals, colours, t):
    retained, y, quotient = retained_minor(adj, terminals, colours)
    assert 1 <= y <= len(terminals)
    for i in range(y, len(quotient)):
        neighbours = set(bits(quotient[i]))
        assert len(neighbours) >= 2 and all(j < y for j in neighbours)
    assert max_k2_minor(quotient) < t, (adj, terminals, colours, quotient, t)
    assert len(retained) <= (t-1)*(y-1)
    return bool(retained)


def contraction_checks():
    graphs = outcomes = positive = 0
    # Every labelled graph on 1--4 vertices; every nonempty terminal set;
    # every allowed colouring (terminals red/blue; other vertices 3 colours).
    for n in range(1, 5):
        possible = list(combinations(range(n), 2))
        for edge_mask in range(1 << len(possible)):
            adj = tuple(from_edges(n, [e for i, e in enumerate(possible)
                                      if edge_mask & (1 << i)]))
            t = max(2, max_k2_minor(adj)+1)
            graphs += 1
            for terminal_mask in range(1, 1 << n):
                terminals = set(bits(terminal_mask))
                allowed = [(RED, BLUE) if v in terminals else (RED, BLUE, WHITE)
                           for v in range(n)]
                for colours in product(*allowed):
                    positive += check_outcome(adj, terminals, colours, t)
                    outcomes += 1
    rng = random.Random(2209202609)
    sampled_graphs = sampled_outcomes = sampled_positive = 0
    for n in (5, 6, 7):
        possible = list(combinations(range(n), 2))
        for _ in range(12):
            adj = tuple(from_edges(n, [e for e in possible if rng.random() < .35]))
            t = max(2, max_k2_minor(adj)+1)
            sampled_graphs += 1
            for _ in range(80):
                terminals = set(rng.sample(range(n), rng.randrange(1, n+1)))
                colours = tuple(rng.randrange(2 if v in terminals else 3)
                                for v in range(n))
                sampled_positive += check_outcome(adj, terminals, colours, t)
                sampled_outcomes += 1
    return dict(exhaustive_graphs=graphs, exhaustive_orders=[1, 2, 3, 4],
                exhaustive_colouring_outcomes=outcomes,
                exhaustive_outcomes_with_retained_white_vertices=positive,
                sampled_graphs=sampled_graphs, sampled_orders=[5, 6, 7],
                sampled_colouring_outcomes=sampled_outcomes,
                sampled_outcomes_with_retained_white_vertices=sampled_positive,
                seed=2209202609,
                independent_oracle='Connected-set partitions; largest K2,s minor')


def exact_probability_checks():
    records, total_outcomes = [], 0
    for h in range(3, 7):
        p = Fraction(h-2, 2*(h-1))
        white = 1-2*p
        for ell in range(2, h+1):
            adj = tuple(from_edges(ell+1, [(v, v+1) for v in range(ell)]))
            terminals = {0, ell}
            masses = [Fraction(0) for _ in range(ell+1)]
            total_mass = Fraction(0)
            for endpoints in product((RED, BLUE), repeat=2):
                for internal in product((RED, BLUE, WHITE), repeat=ell-1):
                    colours = (endpoints[0],) + internal + (endpoints[1],)
                    probability = Fraction(1, 4)
                    for colour in internal:
                        probability *= white if colour == WHITE else p
                    total_mass += probability
                    retained, _, _ = retained_minor(adj, terminals, colours)
                    for v in retained:
                        masses[v] += probability
                    total_outcomes += 1
            assert total_mass == 1
            expected = white*p**(ell-2)/2
            uniform_lower = white*p**(h-2)/2
            for x in range(1, ell):
                assert masses[x] == expected >= uniform_lower
                records.append(dict(h=h, path_length=ell, internal_vertex=x,
                                    exact_probability=str(masses[x]),
                                    required_lower_probability=str(uniform_lower)))
    return dict(exact_probability_equalities=len(records),
                colourings_enumerated=total_outcomes,
                arithmetic='fractions.Fraction; no floating-point tolerances',
                records=records)


def factor(h, t):
    if h == 1:
        return Fraction(1)
    if h == 2:
        return Fraction(t)
    p = Fraction(h-2, 2*(h-1))
    white = 1-2*p
    return 1 + Fraction(2*(t-1), 1)/(white*p**(h-2))


def constructed_closure_checks():
    records = []
    for h in range(3, 7):
        for depth in (h+4, h+5):
            adj = binary_triangulation(depth)
            actual = len(closure(adj, {0, 1}, h))
            bound = 1+factor(h, 3)
            assert actual <= bound < len(adj)
            # Independent direct short-path enumeration above; the ambient
            # graph extends beyond the levels guaranteed by the lower example.
            assert actual >= 2**(h-1)+1
            records.append(dict(h=h, t=3, ambient_depth=depth,
                                graph_order=len(adj), terminals=2,
                                enumerated_closure_size=actual,
                                upper_bound_exact=str(bound),
                                upper_bound_below_ambient_order=True))
    return dict(cases=len(records), nonvacuous_cases=len(records), records=records)


def main():
    started = time.monotonic()
    result = dict(status='passed', scope='Targeted tests of the new colouring proof only',
                  contractions=contraction_checks(),
                  exact_probabilities=exact_probability_checks(),
                  constructed_closures=constructed_closure_checks())
    result['elapsed_seconds'] = round(time.monotonic()-started, 3)
    target = Path(__file__).with_name('colouring_verification_results.json')
    target.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ('exact_probabilities', 'constructed_closures')}, indent=2))
    print('Exact probability equalities:', result['exact_probabilities']['exact_probability_equalities'])
    print('Nonvacuous constructed closures:', result['constructed_closures']['nonvacuous_cases'])


if __name__ == '__main__':
    main()
