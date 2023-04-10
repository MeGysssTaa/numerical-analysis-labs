from typing import Optional, Callable

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


# Общая функция для решения указанной задачи заданным методом и вывода результата.
def solve(
    title: str,  # заголовок с описанием случая для вывода
    I: int,
    x_start: float,
    x_end: float,
    J: int,
    t_start: float,
    t_end: float,
    a: float,
    fun_f: Callable[[float, float], float],
    fun_u_x_0: Callable[[float], float],  # u(x,0) (нужно во всех случаях)
    fun_u_0_x: Optional[Callable[[float], float]],  # u(0,x) (нужно для прямоугольной области при a > 0)
    fun_u_1_x: Optional[Callable[[float], float]],  # u(1,x) (нужно для прямоугольной области при a < 0)
    scheme: Callable,  # схема решения -- одна из функций: scheme_1_halfplane_a_pos, ...
):
    # Подготовка (дано)
    lam = a * r / h
    x: list[float] = [x_start + i * h for i in range(I)]
    t: list[float] = [t_start + j * r for j in range(J)]
    u: list[list[Optional[float]]] = [[None for _ in range(J)] for _ in range(I)]

    for i in range(I):
        u[i][0] = fun_u_x_0(x[i])

    for j in range(J):
        if fun_u_0_x is not None:
            u[0][j] = fun_u_0_x(t[j])
        if fun_u_1_x is not None:
            u[-1][j] = fun_u_1_x(t[j])

    # Решение (применение схемы)
    scheme(u, a, I, x, J, t, lam, fun_f)

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


# Схема 1 для полуплоскости. a > 0
# noinspection PyUnusedLocal
def scheme_1_halfplane_a_pos(u, a, I, x, J, t, lam, fun_f):
    for j in range(J - 1):
        for i in range(j + 1, I):
            u[i][j + 1] = lam * u[i - 1][j] + (1 - lam) * u[i][j] + r * fun_f(x[i], t[j])


# Схема 1 для прямоугольной области. a > 0
# noinspection PyUnusedLocal
def scheme_1_rect_a_pos(u, a, I, x, J, t, lam, fun_f):
    for j in range(J - 1):
        for i in range(0 + 1, I):
            u[i][j + 1] = lam * u[i - 1][j] + (1 - lam) * u[i][j] + r * fun_f(x[i], t[j])


# Схема 2 для полуплоскости. a < 0
# noinspection PyUnusedLocal
def scheme_2_halfplane_a_neg(u, a, I, x, J, t, lam, fun_f):
    for j in range(J - 1):
        for i in range(I - j - 1):
            u[i][j + 1] = -lam * u[i + 1][j] + (1 + lam) * u[i][j] + r * fun_f(x[i], t[j])


# Схема 2 для прямоугольной области. a < 0
# noinspection PyUnusedLocal
def scheme_2_rect_a_neg(u, a, I, x, J, t, lam, fun_f):
    for j in range(J - 1):
        for i in range(I - 0 - 1):
            u[i][j + 1] = -lam * u[i + 1][j] + (1 + lam) * u[i][j] + r * fun_f(x[i], t[j])


# Схема 3 для прямоугольной области. a > 0
# noinspection PyUnusedLocal
def scheme_3_rect_a_pos(u, a, I, x, J, t, lam, fun_f):
    for j in range(J - 1):
        for i in range(0 + 1, I):
            u[i][j + 1] = (u[i][j] + lam * u[i - 1][j + 1] + r * fun_f(x[i], t[j])) / (1 + lam)


# Схема 4 для прямоугольной области. a > 0
# noinspection PyUnusedLocal
def scheme_4_rect_a_pos(u, a, I, x, J, t, lam, fun_f):
    for j in range(J - 1):
        for i in range(0 + 1, I):
            fc = fun_f(x[i] + h / 2, t[j] + r / 2)
            u[i][j + 1] = \
                (u[i - 1][j] * (1 + lam) + (u[i][j] - u[i - 1][j + 1]) * (1 - lam) + 2 * r * fc) \
                / (1 + lam)


# Схема 4 для прямоугольной области. a < 0
# noinspection PyUnusedLocal
def scheme_4_rect_a_neg(u, a, I, x, J, t, lam, fun_f):
    for j in range(J - 1):
        for i in range(I - 1, 0, -1):  # здесь идём по координатам в обратном порядке (от правой границе к левой)
            fc = fun_f(x[i] + h / 2, t[j] + r / 2)
            u[i - 1][j + 1] = \
                (
                        2 * r * h * fc
                        + h * (u[i - 1][j] + u[i][j] - u[i][j + 1])
                        - a * r * (u[i][j] + u[i][j + 1] - u[i - 1][j])
                ) \
                / (h - a * r)


###################################################################################################


def scheme_1_case_1():
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
        fun_f=lambda x, t: 0,
        fun_u_x_0=lambda x: x ** 2 + 2 * x,
        fun_u_0_x=None,
        fun_u_1_x=None,
        scheme=scheme_1_halfplane_a_pos,
    )


def scheme_1_case_2():
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
        fun_f=lambda x, t: x,
        fun_u_x_0=lambda x: x ** 2 + 2 * x,
        fun_u_0_x=None,
        fun_u_1_x=None,
        scheme=scheme_1_halfplane_a_pos,
    )


def scheme_1_case_3():
    J = int((def_t_end - def_t_start) / r) + 1
    I = int((def_x_end - def_x_start) / h) + 1
    solve(
        title="Схема 1. 3 случай: прямоугольная область, a > 0, однородное",
        I=I,
        x_start=def_x_start,
        x_end=def_x_end,
        J=J,
        t_start=def_t_start,
        t_end=def_t_end,
        a=a_pos,
        fun_f=lambda x, t: 0,
        fun_u_x_0=lambda x: x ** 2 + 2 * x,
        fun_u_0_x=lambda t: t ** 2 + 2 * t,
        fun_u_1_x=None,
        scheme=scheme_1_rect_a_pos,
    )


def scheme_1_case_4():
    J = int((def_t_end - def_t_start) / r) + 1
    I = int((def_x_end - def_x_start) / h) + 1
    solve(
        title="Схема 1. 4 случай: прямоугольная область, a > 0, неоднородное",
        I=I,
        x_start=def_x_start,
        x_end=def_x_end,
        J=J,
        t_start=def_t_start,
        t_end=def_t_end,
        a=a_pos,
        fun_f=lambda x, t: x,
        fun_u_x_0=lambda x: x ** 2 + 2 * x,
        fun_u_0_x=lambda t: t ** 2 + 2 * t,
        fun_u_1_x=None,
        scheme=scheme_1_rect_a_pos,
    )


###################################################################################################


def scheme_2_case_1():
    J = int((def_t_end - def_t_start) / r) + 1
    x_end = def_x_end + J * h
    I = int((x_end - def_x_start) / h) + 1
    solve(
        title="Схема 2. 1 случай: полуплоскость, a < 0, однородное",
        I=I,
        x_start=def_x_start,
        x_end=x_end,
        J=J,
        t_start=def_t_start,
        t_end=def_t_end,
        a=a_neg,
        fun_f=lambda x, t: 0,
        fun_u_x_0=lambda x: x ** 2 + 2 * x,
        fun_u_0_x=None,
        fun_u_1_x=None,
        scheme=scheme_2_halfplane_a_neg,
    )


def scheme_2_case_2():
    J = int((def_t_end - def_t_start) / r) + 1
    x_end = def_x_end + J * h
    I = int((x_end - def_x_start) / h) + 1
    solve(
        title="Схема 2. 2 случай: полуплоскость, a < 0, неоднородное",
        I=I,
        x_start=def_x_start,
        x_end=x_end,
        J=J,
        t_start=def_t_start,
        t_end=def_t_end,
        a=a_neg,
        fun_f=lambda x, t: x,
        fun_u_x_0=lambda x: x ** 2 + 2 * x,
        fun_u_0_x=None,
        fun_u_1_x=None,
        scheme=scheme_2_halfplane_a_neg,
    )


def scheme_2_case_3():
    J = int((def_t_end - def_t_start) / r) + 1
    I = int((def_x_end - def_x_start) / h) + 1
    solve(
        title="Схема 2. 3 случай: прямоугольная область, a < 0, однородное",
        I=I,
        x_start=def_x_start,
        x_end=def_x_end,
        J=J,
        t_start=def_t_start,
        t_end=def_t_end,
        a=a_neg,
        fun_f=lambda x, t: 0,
        fun_u_x_0=lambda x: x ** 2 + 2 * x,
        fun_u_0_x=None,
        fun_u_1_x=lambda t: t ** 2 + 2 * t + 3,
        scheme=scheme_2_rect_a_neg,
    )


def scheme_2_case_4():
    J = int((def_t_end - def_t_start) / r) + 1
    I = int((def_x_end - def_x_start) / h) + 1
    solve(
        title="Схема 2. 4 случай: прямоугольная область, a < 0, неоднородное",
        I=I,
        x_start=def_x_start,
        x_end=def_x_end,
        J=J,
        t_start=def_t_start,
        t_end=def_t_end,
        a=a_neg,
        fun_f=lambda x, t: x,
        fun_u_x_0=lambda x: x ** 2 + 2 * x,
        fun_u_0_x=None,
        fun_u_1_x=lambda t: t ** 2 + 2 * t + 3,
        scheme=scheme_2_rect_a_neg,
    )


###################################################################################################


def scheme_3_case_1():
    J = int((def_t_end - def_t_start) / r) + 1
    I = int((def_x_end - def_x_start) / h) + 1
    solve(
        title="Схема 3. 1 случай: прямоугольная область, a > 0, однородное",
        I=I,
        x_start=def_x_start,
        x_end=def_x_end,
        J=J,
        t_start=def_t_start,
        t_end=def_t_end,
        a=a_pos,
        fun_f=lambda x, t: 0,
        fun_u_x_0=lambda x: x ** 2 + 2 * x,
        fun_u_0_x=lambda t: t ** 2 + 2 * t,
        fun_u_1_x=None,
        scheme=scheme_3_rect_a_pos,
    )


def scheme_3_case_2():
    J = int((def_t_end - def_t_start) / r) + 1
    I = int((def_x_end - def_x_start) / h) + 1
    solve(
        title="Схема 3. 2 случай: прямоугольная область, a > 0, неоднородное",
        I=I,
        x_start=def_x_start,
        x_end=def_x_end,
        J=J,
        t_start=def_t_start,
        t_end=def_t_end,
        a=a_pos,
        fun_f=lambda x, t: x,
        fun_u_x_0=lambda x: x ** 2 + 2 * x,
        fun_u_0_x=lambda t: t ** 2 + 2 * t,
        fun_u_1_x=None,
        scheme=scheme_3_rect_a_pos,
    )


###################################################################################################


def scheme_4_case_1():
    J = int((def_t_end - def_t_start) / r) + 1
    I = int((def_x_end - def_x_start) / h) + 1
    solve(
        title="Схема 4. 1 случай: прямоугольная область, a > 0, однородное",
        I=I,
        x_start=def_x_start,
        x_end=def_x_end,
        J=J,
        t_start=def_t_start,
        t_end=def_t_end,
        a=a_pos,
        fun_f=lambda x, t: 0,
        fun_u_x_0=lambda x: x ** 2 + 2 * x,
        fun_u_0_x=lambda t: t ** 2 + 2 * t,
        fun_u_1_x=None,
        scheme=scheme_4_rect_a_pos,
    )


def scheme_4_case_2():
    J = int((def_t_end - def_t_start) / r) + 1
    I = int((def_x_end - def_x_start) / h) + 1
    solve(
        title="Схема 4. 1 случай: прямоугольная область, a > 0, неоднородное",
        I=I,
        x_start=def_x_start,
        x_end=def_x_end,
        J=J,
        t_start=def_t_start,
        t_end=def_t_end,
        a=a_pos,
        fun_f=lambda x, t: x,
        fun_u_x_0=lambda x: x ** 2 + 2 * x,
        fun_u_0_x=lambda t: t ** 2 + 2 * t,
        fun_u_1_x=None,
        scheme=scheme_4_rect_a_pos,
    )


def scheme_4_case_3():
    J = int((def_t_end - def_t_start) / r) + 1
    I = int((def_x_end - def_x_start) / h) + 1
    solve(
        title="Схема 4. 3 случай: прямоугольная область, a < 0, однородное",
        I=I,
        x_start=def_x_start,
        x_end=def_x_end,
        J=J,
        t_start=def_t_start,
        t_end=def_t_end,
        a=a_neg,
        fun_f=lambda x, t: 0,
        fun_u_x_0=lambda x: x ** 2 + 2 * x,
        fun_u_0_x=None,
        fun_u_1_x=lambda t: t ** 2 + 2 * t + 3,
        scheme=scheme_4_rect_a_neg,
    )


def scheme_4_case_4():
    J = int((def_t_end - def_t_start) / r) + 1
    I = int((def_x_end - def_x_start) / h) + 1
    solve(
        title="Схема 4. 4 случай: прямоугольная область, a < 0, неоднородное",
        I=I,
        x_start=def_x_start,
        x_end=def_x_end,
        J=J,
        t_start=def_t_start,
        t_end=def_t_end,
        a=a_neg,
        fun_f=lambda x, t: x,
        fun_u_x_0=lambda x: x ** 2 + 2 * x,
        fun_u_0_x=None,
        fun_u_1_x=lambda t: t ** 2 + 2 * t + 3,
        scheme=scheme_4_rect_a_neg,
    )


###################################################################################################


precision = 6
def_x_start = 0
def_x_end = 1
def_t_start = 0
def_t_end = 10
h = 0.1
r = 0.04
plots_color_theme = "plasma"

a_pos = 2
a_neg = -a_pos


def main():
    # Для красивого вывода чисел в массивах numpy.
    np.set_printoptions(linewidth=100, precision=precision, suppress=True, floatmode="fixed")

    # Для интерактивных графиков в matplotlib.
    # Ещё в настройках PyCharm (Settings -> Tools -> Python Scientific) нужно отключить Show plots in tool window.
    matplotlib.use("TkAgg")

    # Запуск решения одного выбранного случая. Выбрать одну из функций: scheme_1_case_1(), ...
    scheme_4_case_4()


if __name__ == "__main__":
    main()
