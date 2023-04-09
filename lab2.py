from copy import deepcopy

import numpy as np


def err(n, x, x_accurate):
    return [abs(x[row] - x_accurate[row]) for row in range(n)]


def gauss(n, a, b, x_accurate):
    print("\ngauss")

    steps = 0
    mrow = -1

    for start in range(n - 1):
        steps += 1
        for row in range(start, n):
            steps += 1
            if abs(a[row][start]) > a[mrow][start]:
                mrow = row

        a[start], a[mrow] = a[mrow], a[start]
        b[start], b[mrow] = b[mrow], b[start]

        for row in range(start + 1, n):
            steps += 1
            r = a[row][start] / a[start][start]
            for col in range(start, n):
                steps += 1
                a[row][col] -= r * a[start][col]
            b[row] -= r * b[start]

    x = b

    for row in reversed(range(n)):
        steps += 1
        for col in range(row + 1, n):
            steps += 1
            x[row] -= a[row][col] * x[col]
        x[row] /= a[row][row]

    print(f"Solution: {x}")
    print(f"Error: {err(n, x, x_accurate)}")
    print(f"Steps: {steps}")


def seidel(n, a, b, x_accurate):
    print("\nseidel")

    m = np.concatenate((a, np.transpose(np.asarray([b]))), axis=1)
    print(f"Input (bad) matrix):\n{m}")

    m[1] -= 0.67 * m[0]
    m[2] += 1 * m[0]
    m[3] -= 1.35 * m[0]
    m[2] += 0.27 * m[1]
    m[3] += 0.22 * m[2]

    m[2] += 5.48 * m[3]
    m[0] += -0.46 * m[2] + 2.04 * m[3] + 0.2 * m[1]

    print(f"Good matrix:\n{m}")

    x = np.asarray(n * [0.0])  # IMPORTANT: `0.0`, NOT `0`!
    steps = 0
    eps = 10 ** -6

    while True:
        steps += 1
        x_last = x.copy()
        for i in range(n):
            steps += 1
            x[i] = 0.0
            for j in range(n):
                steps += 1
                if i != j:
                    x[i] += x[j] * m[i][j]
            x[i] = (x[i] - m[i][-1]) / -m[i][i]
        if np.linalg.norm(x - x_last) < eps:
            break

    print(f"Solution: {list(x)}")
    print(f"Error: {err(n, x, x_accurate)}")
    print(f"Steps: {steps}")


###################################################################################################


def main():
    np.set_printoptions(suppress=True, precision=4)
    n = 4

    a = [
        [ 0.74, -0.62,  2.11,  0.55],
        [ 0.50,  0.98,  1.79,  0.09],
        [-0.73,  0.25,  2.07,  1.00],
        [ 1.00, -0.85,  1.95,  0.15],
    ]

    b = [
         3.18,
         0.56,
        -2.89,
         5.20,
    ]

    x_accurate = [
         2.0,
        -2.0,
         1.0,
        -3.0,
    ]

    gauss(n, deepcopy(a), deepcopy(b), deepcopy(x_accurate))
    seidel(n, deepcopy(a), deepcopy(b), deepcopy(x_accurate))


if __name__ == '__main__':
    main()
