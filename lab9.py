from typing import Optional, Callable

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def u_x_0(x: float) -> float:
    return 1 - x ** 2


def u_0_t(_t: float) -> float:
    return 1


def u_1_t(_t: float) -> float:
    return 0


def solve(
        title: str,  # заголовок с описанием случая для вывода
        lam: float,
        I: int,
        x_start: float,
        x_end: float,
        method: Callable  # схема решения -- одна из функций: method_1, method_2
):
    # Подготовка (дано)
    x = np.linspace(x_start, x_end, I, dtype=float)

    u = np.zeros(shape=(I, J), dtype=float)
    for i in range(I):
        u[i][0] = u_x_0(x[i])
    for j in range(J):
        u[0][j] = u_0_t(t[j])
        u[-1][j] = u_1_t(t[j])

    # Решение (применение схемы)
    method(u, I, J, lam)

    # Вывод ответа (таблица, график)
    i_plot_start = None
    i_plot_end = None
    for i in range(I):  # ищем, где начинается область, которую нужно вывести (т.к. x мог быть смещён - рисуем не всё)
        if abs(x[i] - a) < 10 ** -6:
            i_plot_start = i
        if abs(x[i] - b) < 10 ** -6:
            i_plot_end = i

    x_plot = [x[i] for i in range(i_plot_start, i_plot_end + 1)]
    t_plot = t
    u_plot = np.array([u[i] for i in range(i_plot_start, i_plot_end + 1)])

    print(title)
    print(u_plot)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    X, T = np.meshgrid(x_plot, t_plot)
    ax.plot_surface(X.T, T.T, u_plot, cmap=plots_color_theme)
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("t")
    ax.set_zlabel("u(x,t)")
    ax.set_xlim(a, b)
    ax.set_ylim(c, d)
    plt.show()


# Схема 1
def method_1(u, I, J, lam):
    for j in range(J - 1):
        for i in range(1, I - 1):
            u[i][j + 1] = lam * u[i + 1][j] + (1 - 2 * lam) * u[i][j] + lam * u[i - 1][j]


# Схема 2
def method_2(u, I, J, lam):
    ai = lam
    bi = -(1 + 2 * lam)
    ci = lam

    for j in range(1, J - 1):
        A: list[Optional[float]] = [None for _ in range(I)]
        B: list[Optional[float]] = [None for _ in range(I)]
        e: list[Optional[float]] = [None for _ in range(I)]
        d: list[Optional[float]] = [None for _ in range(I)]

        d[1] = -u[1][j] - lam * u[0][j + 1]
        d[I - 2] = -u[I - 2][j] - lam * u[I - 1][j + 1]
        A[1] = -ci / bi
        B[1] = d[1] / bi

        for i in range(2, I - 2):
            d[i] = -u[i][j]
            e[i] = ai * A[i - 1] + bi
            A[i] = -ci / e[i]
            B[i] = (d[i] - A[i] * B[i - 1]) / e[i]

        u[I - 2][j + 1] = (d[I - 2] - ai * B[I - 3]) / (bi + ai * A[I - 3])
        for i in range(I - 3, 0, -1):
            u[i][j + 1] = A[i] * u[i + 1][j + 1] + B[i]


###################################################################################################


def solve_with_method_1():
    x_start = a
    x_end = b
    I = int((x_end - x_start) / h) + 1
    lam = r / (h ** 2)  # лямбда
    solve(
        title="Схема 1",
        lam=lam,
        I=I,
        x_start=x_start,
        x_end=x_end,
        method=method_1,
    )


def solve_with_method_2():
    x_start = a
    x_end = b
    I = int((x_end - x_start) / h) + 1
    lam = r / (h ** 2)  # лямбда
    solve(
        title="Схема 2",
        lam=lam,
        I=I,
        x_start=x_start,
        x_end=x_end,
        method=method_2,
    )


###################################################################################################


precision = 6
a = 0  # def_x_start
b = 1  # def_x_end
c = 0  # def_t_start = t_start
d = 10  # def_t_end = t_end
# a^2 (D^2) = 1
# f(x,t) = 0
h = 0.1  # шаг по x
r = 0.00001  # тау - шаг по t
plots_color_theme = "plasma"

J = int((d - c) / r) + 1
t = np.linspace(c, d, J, dtype=float)


def main():
    # Для красивого вывода чисел в массивах numpy.
    np.set_printoptions(linewidth=100, precision=precision, suppress=True, floatmode="fixed")

    # Для интерактивных графиков в matplotlib.
    # Ещё в настройках PyCharm (Settings -> Tools -> Python Scientific) нужно отключить Show plots in tool window.
    matplotlib.use("TkAgg")

    print("...")

    # Решение выбранным способом: solve_with_method_1(), solve_with_method_2()
    solve_with_method_2()


if __name__ == '__main__':
    main()
