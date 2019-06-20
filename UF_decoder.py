from Peeling_decoder import peeling_decoder, get_syndrome, \
    apply_peeling_correction, print_grid_stab, print_grid_qubits
from Toric_code import make_grids, generate_error, check_correction
from collections import defaultdict
import heapq


def find(el, clusters):  # path splitting
    """Returns the root of the cluster of the element with path splitting
    Input:
        el: the element of which the root must be found
        clusters: the clusters in defaultdict form
    Output:
        el: the root of the element"""
    if clusters[el][0] == 0:  # new element, so [element, size, parity] = [el,1,0]
        clusters[el][0] = el
        return el
    while clusters[el][0] != el:
        el, clusters[el][0] = clusters[el][0], clusters[clusters[el][0]][0]
    return el


def union(x, y, clusters):  # by size
    """Merges clusters of x and y
    Input:  x,y : roots of clusters which must be merged
            clusters: dict where all elements link to their parents
    Output: None (clusters are merged in place)"""
    if x == y:  # same root, so no merge
        return
    if clusters[x][1] < clusters[y][1]:
        x, y = y, x  # y is smallest cluster

    clusters[y][0] = x  # make y_root a child of x_root
    clusters[x][1] += clusters[y][1]  # update size
    clusters[x][2] += clusters[y][2]  # update parity


def grow(cluster_roots, boundaries, support, clusters,L):
    """Grows the clusters, merges the touching clusters.
    Input:  cluster_roots:  list of roots of the clusters which must be grown
            boundaries:     dict with the boundary vertices for each cluster_root
            support:        dict with edges as keys, 0, 1 or 2 as value (not grown, half grown, fully grown)
            clusters:       dict where all elements link to their parents
    Output: final_roots:    list of the roots of the odd clusters which have been changed (grown/merged)"""
    # growth
    next_nodes = [(1, 0), (0, 1), (-1, 0), (0, -1)]  # where to go for the neighbouring nodes
    fusion_edges = []
    for u in cluster_roots:  # loop over all given roots
        for v in boundaries[u]:  # each boundary vertex
            for n in next_nodes:  # each possible neighbouring vertex
                other = ((v[0] + n[0]) % L, (v[1] + n[1]) % L)  # define neighbouring vertex
                edge = (min(v, other), max(v, other))  # define edge from boundary vertex to neighbouring vertex
                if support[edge] == 0:  # grow edge
                    support[edge] += 1
                elif support[edge] == 1:
                    support[edge] += 1
                    fusion_edges.append(edge)  # if fully grown: possible merge of clusters

    # print('fusion_edges', fusion_edges)

    # fusion of clusters + fusion of boundary lists + updating roots
    new_roots = defaultdict(int)
    found_roots = defaultdict(int)  # to see which roots we've already seen
    while len(fusion_edges) > 0:  # loop over all fusion edges
        u, v = fusion_edges.pop()  # ensure that the loop will terminate
        # print(u, v)
        u_root = find(u, clusters)
        v_root = find(v, clusters)
        if u_root != v_root:
            found_roots[u_root] += 1  # first check if they are not of the same cluster, then add to found roots.
            found_roots[v_root] += 1
            if clusters[u_root][1] < clusters[v_root][1]:
                u_root, v_root = v_root, u_root  # smallest is v_root
            if len(boundaries[v_root]) == 0:  # v_root is a single generator, no real cluster
                boundaries[u_root].append(v_root)
            else:
                boundaries[u_root].extend(boundaries[v_root])  # v_root has a cluster
            union(u_root, v_root,
                  clusters)  # merge u and v #distinct clusters, thus no need to check in function anymore (2x less find(u) each time)
            if new_roots[
                v_root]:  # v_root is the new root of a merged cluster, but now merging again, to u_root. So it is not a new root anymore
                new_roots[
                    v_root] = 0  # (does not make a difference, only speeds it up a bit, as less things are wrongly returned
            if clusters[u_root][2] % 2 == 1:  # only roots with odd parity added
                new_roots[u_root] += 1  # if two times, then >1, but doesn't matter. keys are roots with odd parity.
                # maybe not true when multiple clusters merge into one?
                # -> doesn't matter, we check if parity is odd before growing another round

    # print('new_roots', new_roots.keys())

    # update boundary lists
    for u in new_roots:  # loop over keys (= roots)
        for v in boundaries[u]:  # for each boundary vertex
            for n in next_nodes:  # if all incident edges are fully grown, remove v from boundary
                other = ((v[0] + n[0]) % L, (v[1] + n[1]) % L)
                edge = (min(v, other), max(v, other))
                if not support[edge] == 2:
                    break
            else:
                boundaries[u].remove(v)

    for x in cluster_roots:  # add original root if not merged
        if not found_roots[x]:
            new_roots[x] += 1  # add roots of non merged clusters (non merged-> same parity -> always add).
            # Note that these boundary lists will stay the same, as they are not merged, which would change the boundary
    final_roots = [x for x in new_roots.keys() if new_roots[x]]  # O(len(new_roots))
    return final_roots


def union_find_decoder(syndrome,L):
    """Initializes all data structures and runs the algorithm
    Input:
        syndrome: the syndrome of an error
    Output:
        the edges which need to be corrected"""
    support = defaultdict(int)  # 0 = unoccupied, 1 = half grown from node, 2 = grown
    boundaries = defaultdict(list)
    clusters = defaultdict(lambda: [0, 1, 0])
    cluster_roots = []
    grow_order = []
    entry_num = 1
    for g in syndrome:
        clusters[g][0] = g  # generator with syndrome
        clusters[g][2] = 1  # parity
        cluster_roots.append(g)  # make g a cluster root
        boundaries[g] = [g]  # boundaries of g
        heapq.heappush(grow_order, [1, entry_num, g])  # add g to the grow order, as parity is odd
        entry_num += 1  # picks lowest in order of entry -> do not grow the same cluster two times in a row if more clusters ¿of same size
    while grow_order:  # until there are no elements to grow
        grow_root = heapq.heappop(grow_order)  # get element with smallest size
        if clusters[grow_root[2]][0] != grow_root[2]:  # not a root anymore, so skip to the next root in grow order
            # print('not a root', grow_root)
            continue
        if len(boundaries[grow_root[2]]) != grow_root[0]:  # boundary has grown, so not smallest anymore
            # print('larger boundary',grow_root)
            continue
        if clusters[grow_root[2]][2] % 2 == 0:  # don't grow if parity is even
            # print('even parity',grow_root)
            continue

        # print('old', grow_root)
        # print('boundary', boundaries[grow_root[1]])
        new_odd_cluster_roots = grow([grow_root[2]], boundaries, support, clusters,L)  # grow smallest cluster
        # print('new:', new_odd_cluster_roots)
        if len(new_odd_cluster_roots) > 0:  # if there are roots to add, add element to grow_order
            for el in new_odd_cluster_roots:
                heapq.heappush(grow_order, [len(boundaries[el]), entry_num,
                                            el])  # possible problem: duplicates. Thus no longer in right order. Fixed by checking at start of loop
                entry_num += 1
        # print('order:', grow_order)
    erasure = []
    for el in support.keys():
        if support[el] == 2:
            erasure.append(el)  # add fully grown edges to erasure

    return peeling_decoder(erasure, syndrome)

def simulate_UF(L,px):
    """Simulate the toric code with the UF decoder.
    Input:
        L: gridsize
        px: the probability on an error
    Output:
        True if the correction is correct
        False if the correction gives a logical error"""
    #make the grids, generate the error and get the syndrome
    grid_g, grid_q = make_grids(L)
    g_errors, q_errors = generate_error(grid_g, grid_q, px)
    syndrome = get_syndrome(g_errors)
    #apply the decoder
    correction = union_find_decoder(syndrome,L)
    grid_corrected = apply_peeling_correction(q_errors, correction)

    correct = check_correction(grid_corrected)
    return correct[0]
