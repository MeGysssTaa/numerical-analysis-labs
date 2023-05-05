import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def f(x, y):
    return x - 3 * y


def u_x_c(x):
    return x + c


def u_x_d(x):
    return x + d


def u_a_y(y):
    return a + y


def u_b_y(y):
    return b + y


def _solve(seidel: bool, grid_size: int, h: float, x, y, u):
    if not seidel:
        u_last = u.copy()  # в методе простой итерации используем только значения, известные с предыдущей итерации
    for i in range(1, grid_size - 1):
        for j in range(1, grid_size - 1):
            if seidel:
                u[i, j] = (u[i + 1, j] + u[i - 1, j] + u[i, j + 1] + u[i, j - 1]) / 4 \
                          + (h ** 2) * f(x[i], y[j])
            else:
                u[i, j] = (u_last[i + 1, j] + u_last[i - 1, j] + u_last[i, j + 1] + u_last[i, j - 1]) / 4 \
                          + (h ** 2) * f(x[i], y[j])
    return u


def solve(seidel: bool, grid_size: int, x, y) -> tuple[object, int]:
    h = (b - a) / (grid_size - 1)
    u = np.zeros(shape=(grid_size, grid_size), dtype=float)
    u[0, :] = u_a_y(y)
    u[-1, :] = u_b_y(y)
    u[:, 0] = u_x_c(x)
    u[:, -1] = u_x_d(x)

    u1 = _solve(seidel, grid_size, h, x, y, u.copy())
    k = 1

    while np.max(np.abs(u - u1)) > eps:
        u = u1
        u1 = _solve(seidel, grid_size, h, x, y, u.copy())
        k += 1

    return u1, k


###################################################################################################


precision = 2
eps = 10 ** (-precision)  # погрешность решения
a = 0  # def_x_start
b = 10  # def_x_end
c = 0  # def_t_start = t_start
d = 10  # def_t_end = t_end
plots_color_theme = "plasma"


def main():
    # Для красивого вывода чисел в массивах numpy.
    np.set_printoptions(linewidth=100, precision=precision, suppress=True, floatmode="fixed")

    # Для интерактивных графиков в matplotlib.
    # Ещё в настройках PyCharm (Settings -> Tools -> Python Scientific) нужно отключить Show plots in tool window.
    matplotlib.use("TkAgg")

    print("...")

    # Решение с выбранными данными.
    grid_size = 10  # размер равномерной сетки (для x и для y)
    seidel = True  # True - использовать метод Зейделя. False - использовать метод простой итерации.

    x = np.linspace(a, b, grid_size)
    y = np.linspace(c, d, grid_size)

    # Вывод ответа.
    u, k = solve(seidel, grid_size, x, y)
    title = "Метод "
    if seidel:
        title += "Зейделя"
    else:
        title += "простой итерации"
    title += f" с размером сетки {grid_size}x{grid_size} (итераций: {k})"
    print(title)
    print(u)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    X, Y = np.meshgrid(x, y)
    ax.plot_surface(X.T, Y.T, u, cmap=plots_color_theme)
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("u(x,y)")
    ax.set_xlim(a, b)
    ax.set_ylim(c, d)
    plt.show()


if __name__ == '__main__':
    main()
