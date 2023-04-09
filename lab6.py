from typing import Optional, Callable

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def solve(
    title: str,
    I: int,
    x_start: float,
    x_end: float,
    J: int,
    t_start: float,
    t_end: float,
    a: float,
    r: float,
    fun_f: Callable[[float, float], float],
    fun_u_x_0: Callable[[float], float],
    fun_u_0_x: Optional[Callable[[float], float]],
    scheme: Callable,
):
    # Подготовка (дано)
    lam = a * r / h
    x: list[float] = [x_start + i * h for i in range(I)]
    t: list[float] = [t_start + j * r for j in range(J)]
    u: list[list[Optional[float]]] = [[None for _ in range(J)] for _ in range(I)]

    for i in range(I):  # u(x,0) - нужно как для полуплоскости, так и для прямоугольной области
        u[i][0] = fun_u_x_0(x[i])

    if fun_u_0_x is not None:  # u(0,x) - нужно только для прямоугольной области
        for j in range(J):
            u[0][j] = fun_u_0_x(t[j])

    f: list[list[float]] = [[fun_f(x[i], t[j]) for j in range(J)] for i in range(I)]  # правая часть

    # Решение (применение схемы)
    scheme(u, I, J, h, r, lam, f)

    # Вывод ответа (таблица, график)
    i_plot_start = None
    i_plot_end = None
    for i in range(I):  # ищем, где начинается область, которую нужно вывести (т.к. x мог быть смещён - рисуем не всё)
        if abs(x[i] - def_x_start) < 10 ** -6:
            i_plot_start = i
        if abs(x[i] - def_x_end) < 10 ** -6:
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
    ax.set_xlim(def_x_start, def_x_end)
    ax.set_ylim(def_t_start, def_t_end)
    plt.show()


# Схема 1 для полуплоскости
def scheme_1_halfplane(u, I, J, h, r, lam, f):
    for j in range(J - 1):
        for i in range(j + 1, I):
            u[i][j + 1] = lam * u[i - 1][j] + (1 - lam) * u[i][j] + r * f[i][j]


# Схема 1 для прямоугольной области
def scheme_1_rect(u, I, J, h, r, lam, f):
    for j in range(J - 1):
        for i in range(1, I):
            u[i][j + 1] = lam * u[i - 1][j] + (1 - lam) * u[i][j] + r * f[i][j]


###################################################################################################


def scheme_1_case_1():
    r = 0.04
    J = int((def_t_end - def_t_start) / r) + 1
    x_start = def_x_start - J * h
    I = int((def_x_end - x_start) / h) + 1
    solve(
        title="Схема 1. 1 случай: полуплоскость, a > 0, однородное",
        I=I,
        x_start=x_start,
        x_end=def_x_end,
        J=J,
        t_start=def_t_start,
        t_end=def_t_end,
        a=a_pos,
        r=r,
        fun_f=lambda _x, _t: 0,
        fun_u_x_0=lambda _x: _x ** 2 + 2 * _x,
        fun_u_0_x=None,
        scheme=scheme_1_halfplane,
    )


def scheme_1_case_2():
    r = 0.04
    J = int((def_t_end - def_t_start) / r) + 1
    x_start = def_x_start - J * h
    I = int((def_x_end - x_start) / h) + 1
    solve(
        title="Схема 1. 2 случай: полуплоскость, a > 0, неоднородное",
        I=I,
        x_start=x_start,
        x_end=def_x_end,
        J=J,
        t_start=def_t_start,
        t_end=def_t_end,
        a=a_pos,
        r=r,
        fun_f=lambda _x, _t: _x,
        fun_u_x_0=lambda _x: _x ** 2 + 2 * _x,
        fun_u_0_x=None,
        scheme=scheme_1_halfplane,
    )


def scheme_1_case_3():
    r = 0.04
    J = int((def_t_end - def_t_start) / r) + 1
    x_start = def_x_start
    I = int((def_x_end - x_start) / h) + 1
    solve(
        title="Схема 1. 3 случай: прямоугольная область, a > 0, однородное",
        I=I,
        x_start=x_start,
        x_end=def_x_end,
        J=J,
        t_start=def_t_start,
        t_end=def_t_end,
        a=a_pos,
        r=r,
        fun_f=lambda _x, _t: 0,
        fun_u_x_0=lambda _x: _x ** 2 + 2 * _x,
        fun_u_0_x=lambda _t: _t ** 2 + 2 * _t,
        scheme=scheme_1_rect,
    )


def scheme_1_case_4():
    r = 0.04
    J = int((def_t_end - def_t_start) / r) + 1
    x_start = def_x_start
    I = int((def_x_end - x_start) / h) + 1
    solve(
        title="Схема 1. 4 случай: прямоугольная область, a > 0, неоднородное",
        I=I,
        x_start=x_start,
        x_end=def_x_end,
        J=J,
        t_start=def_t_start,
        t_end=def_t_end,
        a=a_pos,
        r=r,
        fun_f=lambda _x, _t: _x,
        fun_u_x_0=lambda _x: _x ** 2 + 2 * _x,
        fun_u_0_x=lambda _t: _t ** 2 + 2 * _t,
        scheme=scheme_1_rect,
    )


###################################################################################################


precision = 6
def_x_start = 0
def_x_end = 1
def_t_start = 0
def_t_end = 10
h = 0.1
a_pos = 2
a_neg = -a_pos
plots_color_theme = "plasma"


def main():
    # Для красивого вывода чисел в массивах numpy.
    np.set_printoptions(linewidth=100, precision=precision, suppress=True, floatmode="fixed")

    # Для интерактивных графиков в matplotlib.
    # Ещё в настройках PyCharm (Settings -> Tools -> Python Scientific) нужно отключить Show plots in tool window.
    matplotlib.use("TkAgg")

    # Запуск решения одного выбранного случая.
    scheme_1_case_4()


if __name__ == "__main__":
    main()
