from Toric_code import *
import networkx as nx


def calc_path_lengths(grid_s):
    """TESTED: calculate the path lenghts between the stabilizers which  have an error
    Input:
            grid_s: LxL grid with errors
    Output:
            path_lengths: array with for each element:
             [stabilizer coord 1(tuple), stabilizer coord 2(tuple), path_length]
    """
    # get all stabilizer coords with errors:
    L = len(grid_s)
    stab_errors = []  # which stabilizer measure the errors
    path_lengths = []

    for row_idx in range(len(grid_s)):
        for col_idx in range(len(grid_s[0])):
            if grid_s[row_idx][col_idx] % 2 == 1:  # we only see the error for 1 or 3 qubit errors
                stab_errors.append((row_idx, col_idx))

    # for each stabilizer pair, calculate the minimum path lenght
    for stab1_idx in range(len(stab_errors) - 1):
        for stab2_idx in range(stab1_idx + 1, len(stab_errors)):
            # we have now that stab1_idx< stab2_idx, so each pair occurs once
            # |row_1-row_2| or L-|row_1-row_2|
            min_row_dif = min(abs(stab_errors[stab1_idx][0] - stab_errors[stab2_idx][0]),
                              L - abs(stab_errors[stab1_idx][0] - stab_errors[stab2_idx][0]))

            # |col_1-col_2| or L-|col_1-col_2|:
            min_col_dif = min(abs(stab_errors[stab1_idx][1] - stab_errors[stab2_idx][1]),
                              L - abs(stab_errors[stab1_idx][1] - stab_errors[stab2_idx][1]))
            path_lengths.append([stab_errors[stab1_idx], stab_errors[stab2_idx], min_row_dif + min_col_dif])
    return path_lengths


def matching_to_path(matchings, grid_q):
    """TESTED(for 1 matching):Add path of matchings to qubit grid
    input:
        matchings: array with tuples of two matched stabilizers as elements(stabilizer = tuple of coords)
        grid_q: grid of qubits with errors before correction
    output:
        grid_q: grid of qubits with all errors(correction=adding errors)
    """
    L = len(grid_q[0])
    for stab1, stab2 in matchings:
        error_path = [0, 0]
        row_dif = abs(stab1[0] - stab2[0])
        if row_dif > L - row_dif:
            # path through edge
            error_path[0] += 1
        col_dif = abs(stab1[1] - stab2[1])
        if col_dif > L - col_dif:
            # path through edge
            error_path[1] += 1
        last_row = stab1[0]
        if stab1[0] != stab2[0]:  # not the same row
            up_stab = min(stab1, stab2)
            down_stab = max(stab1, stab2)
            q_col = up_stab[1]  # column of the upper stabilizer
            last_row = down_stab[0]
            if error_path[0]:  # through edge
                for s_row in range(down_stab[0] - L, up_stab[0]):
                    q_row = (s_row + 1) * 2  # row under current stabilizer
                    grid_q[q_row][q_col] += 1  # make error = flip bit
            else:
                for s_row in range(up_stab[0], down_stab[0]):
                    q_row = (s_row + 1) * 2  # row under current stabilizer
                    grid_q[q_row][q_col] += 1

        if stab1[1] != stab2[1]:  # not the same col
            left_stab = min(stab1, stab2, key=lambda x: x[1])
            right_stab = max(stab1, stab2, key=lambda x: x[1])
            q_row = 2 * last_row + 1
            if error_path[1]:  # through edge
                for s_col in range(right_stab[1] - L, left_stab[1]):
                    q_col = s_col + 1  # col right of stabilizer
                    grid_q[q_row][q_col] += 1  # make error = flip bit
            else:
                for s_col in range(left_stab[1], right_stab[1]):
                    q_col = s_col + 1  # col right of stabilizer
                    grid_q[q_row][q_col] += 1
    return grid_q

def simulate_MWPM(L,px):
    """Simulate the toric code with the MWPM decoder, and return the result
    Input:
        L: grid size
        px: probability on an X error (0<=px<=1)
    Output:
        True if correction correct
        False if correction gives logical error
    """
    #make the grids and generate the error
    g, q = make_grids(L)
    g_errors, q_errors = generate_error(g, q, px)

    # get the graph from the stabilizer grid
    path_lengths = calc_path_lengths(g_errors)
    G = nx.Graph()
    for edge in path_lengths:
        G.add_edge(edge[0], edge[1], weight=-edge[2])
    # decode
    matching = nx.algorithms.matching.max_weight_matching(G, maxcardinality=True)
    matched_error_grid = matching_to_path(matching, q_errors)

    # check if decoding worked
    check = check_correction(matched_error_grid)
    return check[0]


