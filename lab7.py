from typing import Optional, Callable

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def u_x_0(x: float) -> float:
    return 3 if x < 0.5 else 2


def solve(
    title: str,  # заголовок с описанием случая для вывода
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

    # Решение (применение схемы)
    method(u, I, J)

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


# Метод введения искусственной вязкости
def method_1(u, I, J):
    for j in range(J - 1):
        for i in range(j + 1, I - j - 1):
            u[i][j + 1] = u[i][j] - r / h * u[i][j] * (u[i][j] - u[i - 1][j]) \
                          - (
                                  (z ** 2) * r / (2 * (h ** 3))
                                  * (u[i + 1][j] - u[i - 1][j])
                                  * (u[i + 1][j] - u[i][j] + u[i - 1][j])
                          )


# Консервативная схема
def method_2(u, I, J):
    for j in range(J - 1):
        for i in range(j + 1, I):
            u[i][j + 1] = r * (u[i - 1][j] ** 2 - u[i][j] ** 2) / (2 * h) + u[i][j]


###################################################################################################


def solve_with_method_1():
    x_start = a - J * h
    x_end = b + J * h
    I = int((x_end - x_start) / h) + 1
    solve(
        title="Метод введения искусственной вязкости",
        I=I,
        x_start=x_start,
        x_end=x_end,
        method=method_1,
    )


def solve_with_method_2():
    x_start = a - J * h
    x_end = b
    I = int((x_end - x_start) / h) + 1
    solve(
        title="Консервативная схема",
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
d = 1  # def_t_end = t_end
h = 0.01  # шаг по x
r = 0.001  # тау - шаг по t
z = 0.01  # эпсилон - для метода введения искусственной вязкости
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
    solve_with_method_1()


if __name__ == '__main__':
    main()
