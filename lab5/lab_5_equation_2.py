from math import cos, sin
from typing import Callable

import matplotlib.pyplot as plt
import numpy as np


def fmt(x: float) -> str:
    result = str(round(x, precision))
    if "." not in result:
        result += ".0"
    return result if "e" in result else result + "0" * max(0, (precision - len(result.split(".")[1])))


def f(_x, _y, z):
    return z


def g(x, y, _z):
    return 0.6 * sin(x) - 1.25 * (y ** 2)


def adams3_rule_ode2(n: int, h: float) -> tuple[list[float], list[float]]:
    # Решаем ОДУ второго порядка
    #     y'' = g(x,y)
    # как систему двух ОДУ первого порядка
    #     { y' = f(x,y,z),
    #     { z' = g(x,y,z),
    # где в нашем случае
    #     f(x,y,z) = z,
    #     g(x,y,z) = 0.6 * sin(x) - 1.25 * (y ** 2).

    # Значения y[0], z[0] даны по условию.
    x = [a + i * h for i in range(n)]
    y = n * [0.0]
    y[0] = 0  # y(a) = 0
    z = n * [0.0]
    z[0] = 1  # y'(a) = 1
    dz = n * [0.0]

    # Находим y[1], z[1], y[2], z[2] методом РК4 для систем двух ОДУ первого порядка.
    for i in range(2):
        k1 = h * f(x[i], y[i], z[i])
        l1 = h * g(x[i], y[i], z[i])
        k2 = h * f(x[i] + h / 2, y[i] + k1 / 2, z[i] + l1 / 2)
        l2 = h * g(x[i] + h / 2, y[i] + k1 / 2, z[i] + l1 / 2)
        k3 = h * f(x[i] + h / 2, y[i] + k2 / 2, z[i] + l2 / 2)
        l3 = h * g(x[i] + h / 2, y[i] + k2 / 2, z[i] + l2 / 2)
        k4 = h * f(x[i] + h, y[i] + k3, z[i] + l3)
        l4 = h * g(x[i] + h, y[i] + k3, z[i] + l3)
        y[i + 1] = y[i] + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        z[i + 1] = z[i] + (l1 + 2 * l2 + 2 * l3 + l4) / 6

    for i in range(3):
        dz[i] = g(x[i], y[i], z[i])

    # Находим оставшиеся значения (y[3], z[3], ..., y[n - 1], z[n - 1]) методом Адамса 3-го порядка.
    for i in range(2, n - 1):
        z[i + 1] = z[i] + h * (23 * dz[i] - 16 * dz[i - 1] + 5 * dz[i - 2]) / 12
        dz[i + 1] = g(x[i + 1], y[i + 1], z[i + 1])
        y[i + 1] = y[i] + h * (23 * z[i] - 16 * z[i - 1] + 5 * z[i - 2]) / 12

    return x, y


def adams4_rule_ode2(n: int, h: float) -> tuple[list[float], list[float]]:
    # Решаем ОДУ второго порядка
    #     y'' = g(x,y)
    # как систему двух ОДУ первого порядка
    #     { y' = f(x,y,z),
    #     { z' = g(x,y,z),
    # где в нашем случае
    #     f(x,y,z) = z,
    #     g(x,y,z) = 0.6 * sin(x) - 1.25 * (y ** 2).

    # Значения y[0], z[0] даны по условию.
    x = [a + i * h for i in range(n)]
    y = n * [0.0]
    y[0] = 0  # y(a) = 0
    z = n * [0.0]
    z[0] = 1  # y'(a) = 1
    dz = n * [0.0]

    # Находим y[1], z[1], ..., y[3], z[3] методом РК4 для систем двух ОДУ первого порядка.
    for i in range(3):
        k1 = h * f(x[i], y[i], z[i])
        l1 = h * g(x[i], y[i], z[i])
        k2 = h * f(x[i] + h / 2, y[i] + k1 / 2, z[i] + l1 / 2)
        l2 = h * g(x[i] + h / 2, y[i] + k1 / 2, z[i] + l1 / 2)
        k3 = h * f(x[i] + h / 2, y[i] + k2 / 2, z[i] + l2 / 2)
        l3 = h * g(x[i] + h / 2, y[i] + k2 / 2, z[i] + l2 / 2)
        k4 = h * f(x[i] + h, y[i] + k3, z[i] + l3)
        l4 = h * g(x[i] + h, y[i] + k3, z[i] + l3)
        y[i + 1] = y[i] + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        z[i + 1] = z[i] + (l1 + 2 * l2 + 2 * l3 + l4) / 6

    for i in range(3):
        dz[i] = g(x[i], y[i], z[i])

    # Находим оставшиеся значения (y[4], z[4], ..., y[n - 1], z[n - 1]) методом Адамса 4-го порядка.
    for i in range(3, n - 1):
        z[i + 1] = z[i] + h * (55 * dz[i] - 59 * dz[i - 1] + 37 * dz[i - 2] - 9 * dz[i - 3]) / 12
        dz[i + 1] = g(x[i + 1], y[i + 1], z[i + 1])
        y[i + 1] = y[i] + h * (55 * z[i] - 59 * z[i - 1] + 37 * z[i - 2] - 9 * z[i - 3]) / 12

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


def main_adams3_rule():
    print("adams3_rule_ode2")
    x, y, steps, n, d = solve(adams3_rule_ode2)
    print(f"\tn = {n}")
    print(f"\td = {d}")
    print(f"\tx        |    y")
    for i in range(len(x)):
        print(f"\t{fmt(x[i])}    |    {fmt(y[i])}")
    fig, ax = plt.subplots()
    ax.plot(x, y, linewidth=2.0)
    plt.show()


def main_adams4_rule():
    print("adams4_rule_ode2")
    x, y, steps, n, d = solve(adams4_rule_ode2)
    print(f"\tn = {n}")
    print(f"\td = {d}")
    print(f"\tx        |    y")
    for i in range(len(x)):
        if i == 0 or i == len(x) - 1 or i % 200 == 0:
            print(f"\t{fmt(x[i])}    |    {fmt(y[i])}")
            if i != len(x) - 1:
                print(f"\t{19 * '.'}")

    fig, ax = plt.subplots()
    ax.plot(x, y, linewidth=2.0)
    plt.show()


if __name__ == '__main__':
    main_adams4_rule()
