#!/usr/bin/env python3
# By Sebastian Raaphorst, 2025.

from copy import deepcopy
from itertools import product
from sys import argv
from typing import Optional, TypeAlias, TypeVar

# Support graphs of arbitrary type.
T = TypeVar('T')
adjacency_list: TypeAlias = dict[T, set[T]]


def _build_de_bruijn_digraph(k: int, n: int) -> adjacency_list[tuple[int, ...]]:
    """
    Constructs a De Bruijn digraph for a given alphabet size `k` and length `n`.
    Note that this means that the vertices correspond to tuples over `range(k)` pf size `n-1`.

    :param k: alphabet size
    :param n: length of substrings represented by the digraph
    :return: A dictionary representing the de Bruijn digraph where keys are vertices (n-1 tuples)
             and values are lists of edges (n-mers).
    :rtype: Adjacency list of the digraph.
    """
    alphabet = frozenset(range(k))
    vertices = frozenset(product(alphabet, repeat=n-1))
    return {v: {v[1:] + (s,) for s in alphabet} for v in vertices}


def _reverse_digraph(orig_digraph: adjacency_list[T]) -> adjacency_list[T]:
    """
    Given a digraph, reverse all the edges and return a mutable digraph.
    :param orig_digraph: the adjacency list of the original digraph
    :return: the reversed digraph
    """
    reversed_digraph: adjacency_list[T] = {}

    for v, nbrs in orig_digraph.items():
        if v not in reversed_digraph:
            reversed_digraph[v] = set()
        for w in nbrs:
            if w not in reversed_digraph:
                reversed_digraph[w] = set()
            reversed_digraph[w].add(v)

    return reversed_digraph


def _degree_condition(digraph: adjacency_list[T]) -> bool:
    """
    Check that the in-degree count of every vertex is equal to the out-degree count.
    :param digraph: the adjacency list of the graph or digraph
    :type digraph: the type of the graph's vertices
    :return: true if the condition is satisfied, false otherwise
    :rtype: bool
    """
    in_degree = {v: 0 for v in digraph}
    out_degree = {v: 0 for v in digraph}

    for v, successors in digraph.items():
        out_degree[v] = len(successors)
        for vp in successors:
            in_degree[vp] += 1

    return all(in_degree[v] == out_degree[v] for v in digraph)


def _dfs_iter(start: T, digraph: adjacency_list[T]) -> set[T]:
    """
    Perform an iterative DFS (i.e using a stack) from start to determine
    the set of vertices that can be reached from in the digraph.
    :param start: A starting vertex in the digraph.
    :param digraph: The adjacency list of the digraph.
    :return: The set of vertices that can be reached from start in the digraph.
    """
    stack = [start]
    visited = set(start)

    while stack:
        current = stack.pop()
        for nxt in digraph[current]:
            if nxt not in visited:
                visited.add(nxt)
                stack.append(nxt)
    return visited


def _is_strongly_connected(digraph: adjacency_list[T]) -> bool:
    """
    Check if digraph is strongly connected using two stack-based DFS passes (original & reversed).
    """
    if not digraph:
        return True

    # Pick a start node.
    start = next(iter(digraph))

    # 1. Visit all nodes in original graph.
    visited_original = _dfs_iter(start, digraph)
    if len(visited_original) < len(digraph):
        return False

    # 2. Reverse the graph
    graph_rev = _reverse_digraph(digraph)

    # 3. Visit all nodes in the reversed graph.
    visited_reversed = _dfs_iter(start, graph_rev)
    if len(visited_reversed) < len(digraph):
        return False

    return True


def _eulerian_cycle(graph: adjacency_list[T]) -> Optional[list[T]]:
    """
    Using Hierholzer's Algorithm, find an Eulerian cycle in the given graph (simple or directed)..
    :param graph: the adjacency list of the graph or digraph
    :type graph: the type of the graph's vertices
    :return: a list of vertices representing the Eulerian cycle, or None if no such cycle exists
    :rtype: the type of the graph's vertices
    """
    # First determine if an Eulerian circuit exists.
    if not _degree_condition(graph):
        return None
    if not _is_strongly_connected(graph):
        return None

    # Find the Eulerian circuit using Hierholzer's Algorithm.
    # We need a version of the adjacency list we can mutate, so we start by copying.
    graph_cpy = deepcopy(graph)
    stack = [next(iter(graph_cpy))]
    cycle = []
    while stack:
        v = stack[-1]
        if graph_cpy[v]:
            w = graph_cpy[v].pop()
            stack.append(w)
        else:
            cycle.append(stack.pop())
    cycle.reverse()
    return cycle


def de_bruijn_sequence(k: int, n: int) -> Optional[list[int]]:
    assert k >= 2, "k must be at least 2."
    assert n >= 3, "n must be at least 3."

    graph = _build_de_bruijn_digraph(k, n)
    e_cycle = _eulerian_cycle(graph)
    sequence = list(e_cycle[0]) + [e[-1] for e in e_cycle[1:]]

    # The sequence should be of size k^n + n - 1:
    # This is because there are k^n keywords, and every window of size n defines one.
    # We have to have the initial window, which must be of size n, so we subtract n-1
    # as the first n characters define the first covering, so n-1 of those are simply
    # part of the sequence cycle. The last n-1 characters should be the same as the first n-1.
    full_cycle_length = k**n + n - 1
    repeated_cycle_length = n - 1
    assert len(sequence) == full_cycle_length, f"Incorrect cycle length: Expected {full_cycle_length}, got {len(sequence)}."
    assert sequence[:repeated_cycle_length] == sequence[-repeated_cycle_length:], f"Does not cycle: {sequence}."

    # Check the cycle to make sure all elements are covered.
    uncovered_elements = set(product(range(k), repeat=n)) - set(tuple(sequence[i:i+n]) for i in range(k**n))
    assert len(uncovered_elements) == 0, f"Uncovered elements: {uncovered_elements}."

    # Chop the unnecessary end off the cycle and output it.
    return sequence[:-repeated_cycle_length]


def _main(k: int, n: int):
    print(de_bruijn_sequence(k, n))


if __name__ == '__main__':
    if len(argv) != 3:
        print(f'Usage: {argv[0]} k n, where')
        print('\tk is the alphabet size, where k ≥ 2')
        print('\tn is the length of substrings, where n ≥ 3')
        exit(1)
    k = int(argv[1])
    n = int(argv[2])
    _main(k, n)
