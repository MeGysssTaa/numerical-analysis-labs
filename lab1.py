import math


a = 0.5
b = 1.0
n = 5
precision = 7
eps = 10 ** (-precision)
m = 6.718279


def f(x):
    return math.e ** x - 2 * (x - 2) ** 2


def fmt(x):
    result = str(round(x, precision))
    return result + "0" * max(0, (precision - len(result.split(".")[1])))


def func_vals():
    print("func_vals")
    for i in range(n + 1):
        x = a + (b - a) / n * i
        print(f"    {i + 1}. {fmt(x)} : {fmt(f(x))}")


def newton():
    print("newton")
    c_last = math.inf
    c = a  # c[k-1]
    k = 0  # c[k]
    while abs(c - c_last) >= eps:
        k += 1
        c_last = c
        c = c - f(c) / (math.e ** c - 4 * c + 8)
        print(f"    {k}. {fmt(c)}")


def chord():
    print("chord")
    c_last = math.inf  # c[k-1]
    c = a - f(a) * (b - a) / (f(b) - f(a))  # c[k]
    k = 0
    while abs(c - c_last) >= eps:
        k += 1
        c_last = c
        c = c - f(c) * (c - a) / (f(c) - f(a))
        print(f"    {k}. {fmt(c)}")


def secant():
    print("secant")
    xi = a  # x[k-1]
    xj = b  # x[k]
    k = 0
    while True:
        k += 1
        x = xj - f(xj) * (xj - xi) / (f(xj) - f(xi))
        print(f"    {k}. {fmt(x)}")
        if abs(x - xj) < eps:
            break
        xi = xj
        xj = x


def newton_fin_diff():
    print("newton_fin_diff")
    h = 0.163725
    c_last = math.inf
    c = (a + b) / 2
    k = 0
    while abs(c - c_last) >= eps:
        k += 1
        c_last = c
        c = c - h * f(c) / (f(c + h) - f(c))
        print(f"    {k}. {fmt(c)}")


def steffensen():
    print("steffensen")
    c_last = math.inf
    c = b
    k = 0
    while abs(c - c_last) >= eps:
        k += 1
        c_last = c
        c = c - f(c) ** 2 / (f(c + f(c)) - f(c))
        print(f"    {k}. {fmt(c)}")


def simple_iter():
    print("simple_iter")
    tau = 0.15
    c_last = math.inf
    c = b
    k = 0
    while abs(c - c_last) >= eps:
        k += 1
        c_last = c
        c = c - tau * f(c)
        print(f"    {k}. {fmt(c)}")


###################################################################################################


func_vals()
newton()
chord()
secant()
newton_fin_diff()
steffensen()
simple_iter()
