def aitken_interp(i, j):
    if j - i <= 1:
        v = (ty[i] * (tx[j] - X) - ty[j] * (tx[i] - X)) / (tx[j] - tx[i])
    else:
        v = (aitken_interp(i, j - 1) * (tx[j] - X) - aitken_interp(i + 1, j) * (tx[i] - X)) / (tx[j] - tx[i])
    print(f"\ti = {i}, j = {j}  =>  L({X}) = {round(v, precision)}")
    return v


###################################################################################################


precision = 5
n = 6
tx = [1.00   , 1.08   , 1.20   , 1.27   , 1.31   , 1.38   ]
ty = [1.17520, 1.30254, 1.50946, 1.21730, 1.22361, 1.23470]
X = 1.022
print(f"По интерполяционной схема Эйткена L({X}) = {round(aitken_interp(0, n - 1), precision)}")
