from math import log


def f(x):
    return log(x / 2)


def df(x):
    return 1 / x


def df2(x):
    return -1 / (x ** 2)


def df2num(x, h):
    return (f(x - h) - 2 * f(x) + f(x + h)) / (h ** 2)


def L(x):
    return \
        0.4055 * ((x - 3.1250) * (x - 3.2500) * (x - 3.3750) * (x - 3.5000)) / ((3.0000 - 3.1250) * (3.0000 - 3.2500) * (3.0000 - 3.3750) * (3.0000 - 3.5000)) + \
        0.4463 * ((x - 3.0000) * (x - 3.2500) * (x - 3.3750) * (x - 3.5000)) / ((3.1250 - 3.0000) * (3.1250 - 3.2500) * (3.1250 - 3.3750) * (3.1250 - 3.5000)) + \
        0.4855 * ((x - 3.0000) * (x - 3.1250) * (x - 3.3750) * (x - 3.5000)) / ((3.2500 - 3.0000) * (3.2500 - 3.1250) * (3.2500 - 3.3750) * (3.2500 - 3.5000)) + \
        0.5232 * ((x - 3.0000) * (x - 3.1250) * (x - 3.2500) * (x - 3.5000)) / ((3.3750 - 3.0000) * (3.3750 - 3.1250) * (3.3750 - 3.2500) * (3.3750 - 3.5000)) + \
        0.5596 * ((x - 3.0000) * (x - 3.1250) * (x - 3.2500) * (x - 3.3750)) / ((3.5000 - 3.0000) * (3.5000 - 3.1250) * (3.5000 - 3.2500) * (3.5000 - 3.3750))


def dfl(x, h):
    return (f(x) - f(x - h)) / h


def dfr(x, h):
    return (f(x + h) - f(x)) / h


def dfc(x, hl, hr):
    return (dfl(x, hl) + dfr(x, hr)) / 2


def aitken_interp(x, _i, _j):
    if _j - _i <= 1:
        v = (ty[_i] * (tx[_j] - x) - ty[_j] * (tx[_i] - x)) / (tx[_j] - tx[_i])
    else:
        v = (aitken_interp(x, _i, _j - 1) * (tx[_j] - x) - aitken_interp(x, _i + 1, _j) * (tx[_i] - x)) / (tx[_j] - tx[_i])
    print(f"\ti = {_i}, j = {_j}  =>  L({x}) = {round(v, precision)}")
    return v


###################################################################################################


precision = 4
a = 3.0000
b = 3.5000
m = 3.0300
n = 5

print(f"Таблица значений функции f(x):")
H = (b - a) / (n - 1)
for i in range(n):
    print(f"\t{a + i * H}\t\t{round(f(a + i * H), precision)}")
print()

tx = [3.0000, 3.1250, 3.2500, 3.3750, 3.5000]
ty = [0.4055, 0.4463, 0.4855, 0.5232, 0.5596]


print("\n\n")
print(f"Значение f(m) = f({m}) (интерполирование по схеме Эйткена): "
      f"L({m}) = {round(aitken_interp(m, 0, n - 1), precision)}")
print("\n\n")


print(f"Производные в узлах:")

for i in range(n):
    print(f"\tx[{i}] = {tx[i]}")

    df_i = df(tx[i])
    df2_i = df2(tx[i])

    dfl_i = dfl(tx[i], tx[i] - tx[i - 1]) if i > 0 else None
    dfr_i = dfr(tx[i], tx[i + 1] - tx[i]) if i < n - 1 else None
    dfc_i = (dfl_i + dfr_i) / 2 if 0 < i < n - 1 else None
    df2num_i = df2num(tx[i], H) if 0 < i < n - 1 else None
    
    print(f"\t\tточная первая = {round(df_i, precision)}")
    print(f"\t\tточная вторая = {round(df2_i, precision)}")
    
    if dfl_i is None:
        print(f"\t\tслева НЕТ")
    else:
        print(f"\t\tслева = {round(dfl_i, precision)}, погрешность {round(abs(dfl_i - df_i), precision)}")
        
    if dfr_i is None:
        print(f"\t\tсправа НЕТ")
    else:
        print(f"\t\tсправа = {round(dfr_i, precision)}, погрешность {round(abs(dfr_i - df_i), precision)}")
        
    if dfc_i is None:
        print(f"\t\tцентр. НЕТ")
    else:
        print(f"\t\tцентр. = {round(dfc_i, precision)}, погрешность {round(abs(dfc_i - df_i), precision)}")

    if df2num_i is None:
        print(f"\t\tвторая НЕТ")
    else:
        print(f"\t\tвторая = {round(df2num_i, precision)}, погрешность {round(abs(df2num_i - df2_i), precision)}")


print(f"Производные в точке m = {m}")
df_m = df(m)
df2_m = df2(m)
dfl_m = dfl(m, m - 3.0000)
dfr_m = dfr(m, 3.1250 - m)
dfc_m = dfc(m, m - 3.0000, 3.1250 - m)
df2num_m = 0
print(f"\tточная первая = {round(df_m, precision)}")
print(f"\tточная вторая = {round(df_m, precision)}")
print(f"\tслева = {round(dfl_m, precision)}, погрешность {round(abs(dfl_m - df_m), precision)}")
print(f"\tсправа = {round(dfr_m, precision)}, погрешность {round(abs(dfr_m - df_m), precision)}")
print(f"\tцентр. = {round(dfc_m, precision)}, погрешность {round(abs(dfc_m - df_m), precision)}")
