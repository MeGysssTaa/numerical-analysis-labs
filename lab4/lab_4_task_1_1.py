def f(x):
    return x + 10 / x


def df(x):
    return 1 - 10 / (x ** 2)


def L(x):
    return \
        6.500 * ((x - 4.300) * (x - 4.600) * (x - 4.900)) / ((4.000 - 4.300) * (4.000 - 4.600) * (4.000 - 4.900)) + \
        6.626 * ((x - 4.000) * (x - 4.600) * (x - 4.900)) / ((4.300 - 4.000) * (4.300 - 4.600) * (4.300 - 4.900)) + \
        6.774 * ((x - 4.000) * (x - 4.300) * (x - 4.900)) / ((4.600 - 4.000) * (4.600 - 4.300) * (4.600 - 4.900)) + \
        6.941 * ((x - 4.000) * (x - 4.300) * (x - 4.600)) / ((4.900 - 4.000) * (4.900 - 4.300) * (4.900 - 4.600))


def dfl(x, h):
    return (f(x) - f(x - h)) / h


def dfr(x, h):
    return (f(x + h) - f(x)) / h


def dfc(x, hl, hr):
    return (dfl(x, hl) + dfr(x, hr)) / 2


###################################################################################################


n = 4
tx = [4.000, 4.300, 4.600, 4.900]
ty = [6.500, 6.626, 6.774, 6.941]
X = 4.390

print(f"Значение функции:\n\t"
      f"f({X}) = {f(X)}\n\t"
      f"L({X}) = {L(X)}\n\t"
      f"Погрешность = {abs(f(X) - L(X))}")
print(f"Значение производной f'({X}):\n\t"
      f"Точное = {df(X)}\n\t"
      f"Прибл. лев. = {dfl(X, X - 4.300)} (погрешность = {abs(dfl(X, X - 4.300) - df(X))})\n\t"
      f"Прибл. прав. = {dfr(X, 4.600 - X)} (погрешность = {abs(dfr(X, 4.600 - X) - df(X))})\n\t"
      f"Прибл. центр. = {dfc(X, X - 4.300, 4.600 - X)} (погрешность = {abs(dfc(X, X - 4.300, 4.600 - X) - df(X))}))")
