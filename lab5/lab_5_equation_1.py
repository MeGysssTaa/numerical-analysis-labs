from math import cos

import matplotlib.pyplot as plt
import numpy as np


def fmt(x: float) -> str:
    result = str(round(x, precision))
    if "." not in result:
        result += ".0"
    return result if "e" in result else result + "0" * max(0, (precision - len(result.split(".")[1])))


def f(x, y):
    return cos(x) - 0.5 * (y ** 2)


def heun_rule(n: int, h: float) -> tuple[list[float], list[float]]:
    x = [a + i * h for i in range(n)]
    y = n * [0.0]
    y[0] = 0
    dy = n * [0.0]
    dy[0] = f(x[0], y[0])

    for i in range(n - 1):
        y[i + 1] = y[i] + h * dy[i]
        dy[i + 1] = f(x[i], y[i + 1])
        y[i + 1] = y[i] + h * (dy[i] + dy[i + 1]) / 2

    return x, y


def runge_kutta_rule(n: int, h: float) -> tuple[list[float], list[float]]:
    x = [a + i * h for i in range(n)]
    y = n * [0.0]
    y[0] = 0
    dy = n * [0.0]
    dy[0] = f(x[0], y[0])

    for i in range(n - 1):
        k1 = h * f(x[i], y[i])
        k2 = h * f(x[i] + h / 2, y[i] + k1 / 2)
        k3 = h * f(x[i] + h / 2, y[i] + k2 / 2)
        k4 = h * f(x[i] + h, y[i] + k3)
        delta = (k1 + 2 * k2 + 2 * k3 + k4) / 6
        y[i + 1] = y[i] + delta

    return x, y


def max_dist(x1: list[float], y1: list[float], x2: list[float], y2: list[float]) -> float:
    d = 0.0
    for i in range(len(x1)):
        j = None
        for k in range(len(x2)):
            if abs(x1[i] - x2[k]) < 10 ** -6:
                j = k
                break
        if j:
            d = max(d, abs(y1[i] - y2[j]))
    return d


def solve(rule, n=4, steps=1, x_last=None, y_last=None) -> tuple[list[float], list[float], int, int, float]:
    h = (b - a) / (n - 1)
    x, y = rule(n, h)
    d = None if x_last is None else max_dist(x, y, x_last, y_last)
    if d is not None and d < eps:
        return x, y, steps, n, d
    return solve(rule, 2 * n, steps + 1, x, y)


###################################################################################################


# Для всех вариантов и уравнений y(a) = 0, [a, b] = [0; 0,5], для уравнения 2 – y’(a) = 1. Погрешность решения 0,001.
a = 0.0
b = 0.5
precision = 3
eps = 10 ** -precision


def main_heun_rule():
    print("heun_rule")
    x, y, steps, n, d = solve(heun_rule)
    print(f"\tn = {n}")
    print(f"\td = {d}")
    print(f"\tx        |    y")
    for i in range(len(x)):
        print(f"\t{fmt(x[i])}    |    {fmt(y[i])}")
    fig, ax = plt.subplots()
    ax.plot(x, y, linewidth=2.0)
    plt.show()


def main_runge_kutta_rule():
    print("runge_kutta_rule")
    x, y, steps, n, d = solve(runge_kutta_rule)
    print(f"\tn = {n}")
    print(f"\td = {d}")
    print(f"\tx        |    y")
    for i in range(len(x)):
        print(f"\t{fmt(x[i])}    |    {fmt(y[i])}")
    fig, ax = plt.subplots()
    ax.plot(x, y, linewidth=2.0)
    plt.show()


if __name__ == '__main__':
    main_runge_kutta_rule()
