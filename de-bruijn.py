#!/usr/bin/env python3
# By Sebastian Raaphorst, 2025.

from itertools import product
from sys import argv
from typing import Optional, TypeAlias, TypeVar

# Support graphs of arbitrary type.
T = TypeVar('T')
adjacency_list: TypeAlias = dict[T, frozenset[T]]


def build_de_bruijn_graph(k: int, n: int) -> adjacency_list[tuple[int, ...]]:
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
    adjacencies = {v: [] for v in vertices}

    for v in vertices:
        # If v = (t_0, t_1, ..., t_{n-2}),, make out edges with all vertices of the form
        # w = (t_1, ..., t_{n-2}, s).
        for s in alphabet:
            w = v[1:] + (s,)
            adjacencies[v].append(w)

    return {v: frozenset(adjacencies[v]) for v in vertices}


def _degree_condition(graph: adjacency_list[T]) -> bool:
    """
    Check that the in-degree count of every vertex is equal to the out-degree count.
    :param graph: the adjacency list of the graph or digraph
    :type graph: the type of the graph's vertices
    :return: true if the condition is satisfied, false otherwise
    :rtype: bool
    """
    in_degree = {v: 0 for v in graph}
    out_degree = {v: 0 for v in graph}

    for v, successors in graph.items():
        out_degree[v] = len(successors)
        for vp in successors:
            in_degree[vp] += 1

    return all(in_degree[v] == out_degree[v] for v in graph)


def _strongly_connected(graph: adjacency_list[T]) -> bool:
    """
    Determines if a directed graph is strongly connected.

    We do so by using Kosaraju's algorithm to find strongly connected components
    and make sure that there is only one such component.

    Link to the pseudocode:
    https://en.wikipedia.org/wiki/Kosaraju%27s_algorithm

    :param graph: The directed graph represented as an adjacency list where
        each vertex maps to a list of adjacent vertices.
    :type graph: the type of the graph's vertices
    :return: A boolean indicating whether the graph is strongly connected.
    :rtype: bool
    """
    # Step 1 of Kosaraju's algorithm.
    visited = {v: False for v in graph}
    block = []

    # Step 2 of Kosaraju's algorithm.
    def visit(u: T) -> None:
        if not visited[u]:
            visited[u] = True
            for v in graph[u]:
                visit(v)
            block.insert(0, u)

    for u in graph:
        visit(u)

    # Step 3 of Kosaraju's algorithm.
    components: dict[T, T] = {}
    def assign(u: T, root: T):
        if u not in components:
            components[u] = root
            for v in graph:
                if u in graph.get(v, frozenset()):
                    assign(v, root)

    for u in block:
        assign(u, u)

    # Make sure there is only one strongly connected component.
    return len(set(components.values())) == 1


def _hierholzer(graph: adjacency_list[T]) -> list[T]:
    # Make a copy of graph, but we need it to be mutable.
    graph_mut = {v: set(nbrs) for v, nbrs in graph.items()}

    # We want to start with the vertex that is all 0s.
    stack = [(0,) * len(next(iter(graph)))]
    # stack = [next(iter(graph_mut))]
    cycle = []
    while stack:
        v = stack[-1]
        if graph_mut[v]:
            w = graph_mut[v].pop()
            stack.append(w)
        else:
            cycle.append(stack.pop())
    cycle.reverse()
    return cycle


def eulerian_cycle(graph: adjacency_list[T]) -> Optional[list[T]]:
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
    if not _strongly_connected(graph):
        return None
    return _hierholzer(graph)


def main(k: int, n: int):
    graph = build_de_bruijn_graph(k, n)
    e_cycle = eulerian_cycle(graph)
    de_bruijn_sequence = list(e_cycle[0]) + [e[-1] for e in e_cycle[1:]]

    # If there is a zero at the end of the de Bruijn sequence, move it to the beginning.
    if de_bruijn_sequence[-1] == 0:
        de_bruijn_sequence = [0] + de_bruijn_sequence[:-1]

    # The sequence should be of size k^n + n - 1:
    # This is because there are k^n keywords, and every window of size n defines one.
    # We have to have the initial window, which must be of size n, so we subtract n-1
    # as the first n characters define the first covering, so n-1 of those are simply
    # part of the sequence cycle. The last n-1 characters should be the same as the first n-1.
    full_cycle_length = k**n + n - 1
    repeated_cycle_length = n - 1
    assert len(de_bruijn_sequence) == full_cycle_length, f"Incorrect cycle length: Expected {full_cycle_length}, got {len(de_bruijn_sequence)}."
    assert de_bruijn_sequence[:repeated_cycle_length] == de_bruijn_sequence[-repeated_cycle_length:], f"Does not cycle: {de_bruijn_sequence}."

    # Check the cycle to make sure all elements are covered.
    uncovered_elements = set(product(range(k), repeat=n)) - set(tuple(de_bruijn_sequence[i:i+n]) for i in range(k**n))
    assert len(uncovered_elements) == 0, f"Uncovered elements: {uncovered_elements}."

    # Chop the unnecessary end off the cycle and output it.
    de_bruijn_sequence = de_bruijn_sequence[:-repeated_cycle_length]
    print(de_bruijn_sequence)


if __name__ == '__main__':
    if len(argv) != 3:
        print(f'Use: {argv[0]} k n, where')
        print('\tk is the alphabet size')
        print('\tn is the length of substrings')
        exit(1)

    k = int(argv[1])
    n = int(argv[2])
    main(k, n)
