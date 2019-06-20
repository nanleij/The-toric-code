from collections import defaultdict
from random import random
from Toric_code import *


def peeling_decoder(erasure, syndrome):
    """(TESTED) Given an erasure and syndrome, construct a Pauli Z error P
    such that P is a subset of the erasure and such that the syndrome of P equals the syndrome
    1) construct spanning forest F
    2) initialize A (empty set/list)
    3) while F nonempty: pick leaf e = {u,v}. u at the outside
        -remove e from F
        - if u in syndrome: add e to A, remove u from syndrome, flip v in syndrome
        -else: do nothing
    return P = Z errors on all edges in A"""
    # syndrome is defaultdict(int)
    F = spanning_forest_dict(erasure)
    A = []
    vertex_count = defaultdict(int)
    for el in F:
        vertex_count[el[0]] += 1
        vertex_count[el[1]] += 1

    while F:
        leaf_edge = F.pop()
        if vertex_count[leaf_edge[0]] == 0:
            u = leaf_edge[0]
            v = leaf_edge[1]
        else:
            u = leaf_edge[1]
            v = leaf_edge[0]

        vertex_count[u] -= 1
        vertex_count[v] -= 1
        if syndrome[u] == 1:
            A.append(leaf_edge)
            syndrome[u] -= 1  # not necessary? # can be used to check at the end
            if syndrome[v] == 1:
                syndrome[v] -= 1
            else:
                syndrome[v] += 1
    return A


def edge_list_to_graph(e):
    """Creates a graph in the form of adjacency lists from an edge list
    Input:  e:      list of all edges of the graph
    Output: graph:  a dict with vertices as keys, lists of adjacent vertices as values """
    graph = defaultdict(list)
    for e1, e2 in e:
        graph[e1].append(e2)
        graph[e2].append(e1)
    return graph


def spanning_tree_dict(graph, visited, tree, vertex):
    """Creates a spanning tree of the connected component of vertex in the graph
    Input:  graph: dict with vertices as keys, for each vertex the value is a list of neighbouring vertices
            visited: dict with vertices as keys, 0 or 1 as value. 0 =  not visited, 1 = visited
            tree: the tree which is being build
            vertex: the current vertex from where the tree is built further
    Output: tree: the final spanning tree
            visited: dict with vertices as keys, 0 or 1 as value. 0 =  not visited, 1 = visited. Not useful but needed for recursion"""
    visited[vertex] = 1

    for neighbour in graph[vertex]:
        if not visited[neighbour]:
            tree.append((vertex, neighbour))
            tree, visited = spanning_tree_dict(graph, visited, tree, neighbour)
    return tree, visited


def spanning_forest_dict(graph):
    """Creates a spanning forest of the graph:
    Input:  graph: dict with vertices as keys, for each vertex the value is a list of neighbouring vertices
                    or an edge list (list of all edges in the graph)
    Output: tree: a spanning forest of the graph, as a list of edges"""
    if type(graph) != dict:
        # assume its an edge list
        graph = edge_list_to_graph(graph)
    tree = []
    visited = defaultdict(int)
    vertices = graph.keys()
    for v in vertices:
        if not visited[v]:
            tree, visited = spanning_tree_dict(graph, visited, tree, v)
    return tree


def make_erasure(grid_s, grid_q, p_e):
    """"Make an erasure on grid_q, and add the syndrome to grid_s. Erasure with probability p_e,
    if a qubit is erased, 50% chance on X/Y error(detectable)
    Input:
        grid_q: qubit grid
        grid_s: stabilizer grid
        p_e: probability on erasure
    Output:
        erasures: the qubits which are erased (between which stabs)
        grid_s: stabilizer grid with errors added
        grid_q: qubit grid with errors added
        """
    erasures = []
    for row_idx in range(len(grid_q)):
        for col_idx in range(len(grid_q[0])):
            error = random() <= p_e
            if error:
                if row_idx % 2 == 0:
                    # above/under stabilizers -> same column
                    stab_row = int(row_idx / 2) % len(grid_s)
                    erasures.append(((stab_row, col_idx), ((stab_row - 1) % len(grid_s), col_idx)))
                    if random() <= 0.5:  # 50% chance on Y or Z error
                        grid_s[stab_row][col_idx] += 1  # stabilizer under qubit
                        grid_s[(stab_row - 1) % len(grid_s)][col_idx] += 1  # stabilizer above qubit
                        grid_q[row_idx][col_idx] += 1
                else:
                    # left/right of stabilizers -> same row
                    stab_row = int((row_idx - 1) / 2) % len(grid_s)
                    erasures.append(((stab_row, col_idx), (stab_row, (col_idx - 1) % len(grid_s[0]))))
                    if random() <= 0.5:  # 50% chance on Y or Z error
                        grid_s[stab_row][col_idx] += 1  # stabilizer right of qubit
                        grid_s[stab_row][(col_idx - 1) % len(grid_s[0])] += 1  # stabilizer right of qubit
                        grid_q[row_idx][col_idx] += 1
    return erasures, grid_s, grid_q


def get_syndrome(grid_s):
    """Determines the syndrome based on grid_s
    Input: grid_s: stabilizer grid of which the syndrome should be determined
    Output: syndrome_res: dict with 1 as value for the coords of the stabilizers with -1 as outcome"""
    syndrome_res = defaultdict(int)
    for row_idx in range(len(grid_s)):
        for col_idx in range(len(grid_s[0])):
            if grid_s[row_idx][col_idx] % 2 == 1:
                syndrome_res[(row_idx, col_idx)] += 1
    return syndrome_res


def apply_peeling_correction(grid_q, correction):
    """Apply the correction on the grid.
    Input:
        grid_q: qubit grid on which the correction must be applied
        correction: the qubits which need to be flipped, given as (stab1, stab2)
    Output:
        grid_q: corrected qubit grid
    """
    for bit in correction:
        if bit[0][0] == bit[1][0]:  # same row
            row = bit[0][0] * 2 + 1
            left_col = min(bit[0][1], bit[1][1])
            right_col = max(bit[0][1], bit[1][1])
            if right_col - left_col > 1:  # edge case
                col = 0
            else:
                col = right_col
        else:  # same col
            col = bit[0][1]
            up_row = min(bit[0][0], bit[1][0])
            low_row = max(bit[0][0], bit[1][0])
            if low_row - up_row > 1:  # edge case
                row = 0
            else:
                row = low_row * 2
        grid_q[row][col] += 1
    return grid_q


def simulate_peeling(L, pe):
    """Simulate the toric code with erasure errors and the peeling decoder.
    Input:
        L: gridsize
        pe: probability on erasure error
    Output:
        True if correction correct
        False if correction gives logical error
    """
    # make grids, generate erasure error and get the syndrome
    g, q = make_grids(L)
    erasure, error_g, error_q = make_erasure(g, q, pe)
    syndrome = get_syndrome(error_g)

    # apply the peeling decoder
    correction = peeling_decoder(erasure, syndrome)
    corrected_grid = apply_peeling_correction(error_q, correction)

    # check if the correction is correct
    check = check_correction(corrected_grid)
    return check[0]
