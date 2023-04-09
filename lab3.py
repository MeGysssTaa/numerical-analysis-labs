from math import sin, log, tan


def fmt(x):
    result = str(round(x, precision))
    return result if "e" in result else result + "0" * max(0, (precision - len(result.split(".")[1])))


def f(x):
    return x / (sin(3 * x)) ** 2


def F(x):
    return log(sin(3 * x)) / 9 - (x / 3) * (1 / tan(3 * x))


def riemann_left_rule(n, h):
    return h * sum(f(a + i * h) for i in range(n))  # i = 0 .. n-1


def trapezoidal_rule(n, h):
    return h * ((f(a) + f(b)) / 2 + sum(f(a + i * h) for i in range(1, n)))  # i = 1 .. n-1


def simpson_rule(n, h):
    s = f(a) + f(b)
    for i in range(1, n):  # i = 1 .. n-1
        s += 2 * f(a + i * h) if i % 2 == 0 else 4 * f(a + i * h)
    return (h / 3) * s


def integrate(rule, n=4, steps=1, s_last=None):
    h = (b - a) / n
    s = rule(n, h)
    print(f"\tStep {steps}. h = {h}, s = {s}, s_last = {s_last}, "
          f"diff = {abs(s - s_last) if s_last is not None else 'None'}")
    if s_last is not None and abs(s - s_last) < eps:
        return s, steps, n
    return integrate(rule, 2 * n, steps + 1, s)


###################################################################################################


a = 0.2
b = 1
precision = 5
eps = 10 ** (-precision)

s_accurate = F(b) - F(a)
print(f"Accurate integral value: {fmt(s_accurate)}")

print("riemann_left_rule")
res_s, res_steps, res_n = integrate(riemann_left_rule)
print(f"\tAnswer: {fmt(res_s)}. Steps: {res_steps}. Points: {res_n}. "
      f"Relative Error: {fmt(abs((s_accurate - res_s) / s_accurate * 100))}%")

print("trapezoidal_rule")
res_s, res_steps, res_n = integrate(trapezoidal_rule)
print(f"\tAnswer: {fmt(res_s)}. Steps: {res_steps}. Points: {res_n}. "
      f"Relative Error: {fmt(abs((s_accurate - res_s) / s_accurate * 100))}%")

print("simpson_rule")
res_s, res_steps, res_n = integrate(simpson_rule)
print(f"\tAnswer: {fmt(res_s)}. Steps: {res_steps}. Points: {res_n}. "
      f"Relative Error: {fmt(abs((s_accurate - res_s) / s_accurate * 100))}%")
