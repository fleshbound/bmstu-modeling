import numpy as np
import scipy as sp


def get_kolmogorov_eq(matrix):
    count = matrix.shape[0]

    right = np.zeros([count] * 2)

    for state in range(count):
        right[state][state] = -sum(matrix[state, :])
        # t = matrix[:, state]
        # t[state] = 0
        right[state] += matrix[:, state]

    return right


def normalization(statesNumber):
    return 1, np.zeros(statesNumber) + 1


def calculate_p(matrix):
    count = matrix.shape[0]

    left = np.zeros(count)
    right = get_kolmogorov_eq(matrix)

    left[-1], right[-1] = normalization(count)

    return sp.linalg.solve(right, left)


def calculate_time(prob, matrix):
    count = matrix.shape[0]

    time = np.zeros(count)

    for i in range(count):
        time[i] = prob[i] / abs(sum(matrix[: , i]) - sum(matrix[i, :]))

    return time
