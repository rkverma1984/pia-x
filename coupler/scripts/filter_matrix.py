import numpy as np

def filter_main(matrix, xlabel, ylabel, cutoff):
    not_empty_pairs = []

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if not matrix[i][j] < cutoff:
                not_empty_pairs.append((i,j)) 

    not_empty_x = []
    not_empty_y = []
    for n in not_empty_pairs:
        not_empty_x.append(n[0])
        not_empty_y.append(n[1])
    not_empty_x = sorted(list(set(not_empty_x)))
    not_empty_y = sorted(list(set(not_empty_y)))

    filtered_matrix = np.zeros((len(not_empty_x), len(not_empty_y)))
    for i in range(len(not_empty_x)):
        for j in range(len(not_empty_y)):
            xi = not_empty_x[i]
            yj = not_empty_y[j]
            if (xi, yj) in not_empty_pairs:
                filtered_matrix[i][j] = matrix[xi][yj]
    new_xlabel = [xlabel[n] for n in not_empty_x]
    new_ylabel = [ylabel[n] for n in not_empty_y]

    return filtered_matrix, new_xlabel, new_ylabel


