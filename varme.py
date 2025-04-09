import math as m
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Parametere
h = 0.01
k = 0.00005  # da blir k/h^2 = 0.5
L = 1
T = 0.1
nx = int(L / h)
nt = int(T / k)
r = k / h**2

# Initialbetingelse
def f(x):
    return m.sin(m.pi * x)

# Den analytiske løsningen
def u_exact(x, t):
    return m.exp(-m.pi**2 * t) * m.sin(m.pi * x)

u = [f(i*h) for i in range(nx + 1)]
tider_euler_imp = [u.copy()]
tider_crank_nicolson = [u.copy()]
tider_euler_exp = [u.copy()]
tider_analytisk = []

# Thomas-algoritmen for tridiagonalt system
def thomas_algorithm(a, b, c, d):
    n = len(b)
    cp = c.copy()
    bp = b.copy()
    dp = d.copy()

    for i in range(1, n):
        m = a[i-1] / bp[i-1]
        bp[i] = bp[i] - m * cp[i-1]
        dp[i] = dp[i] - m * dp[i-1]

    x = [0] * n
    x[-1] = dp[-1] / bp[-1]
    for i in range(n-2, -1, -1):
        x[i] = (dp[i] - cp[i] * x[i+1]) / bp[i]

    return x

# Løser tidsskritt for tid = 0 til T for de tre metodene
n_indre = nx - 1  # antall indre punkter (uten rand)
for t in range(1, nt):
    # Tridiagonal koeffisienter for Euler Implicit og Crank-Nicolson
    a = [-r] * (n_indre - 1)
    b = [1 + 2*r] * n_indre
    c = [-r] * (n_indre - 1)

    # Euler Implicit
    d_euler_imp = tider_euler_imp[t-1][1:-1]  # høyresiden for implisitt Euler
    u_inner_euler_imp = thomas_algorithm(a, b, c, d_euler_imp)
    u_euler_imp = [0] + u_inner_euler_imp + [0]
    tider_euler_imp.append(u_euler_imp)

    # Crank-Nicolson
    a_cn = [-r/2] * (n_indre - 1)
    b_cn = [1 + r] * n_indre
    c_cn = [-r/2] * (n_indre - 1)
    
    # Beregn d for Crank-Nicolson - håndtere for både gammel og ny tid
    d_cn = []
    for i in range(n_indre):
        if i == 0:
            d_cn.append(0)  # Randbetingelse for venstre punkt
        elif i == n_indre - 1:
            d_cn.append(0)  # Randbetingelse for høyre punkt
        else:
            d_val = (1 - r) * tider_crank_nicolson[t-1][i] + r / 2 * (tider_crank_nicolson[t-1][i+1] + tider_crank_nicolson[t-1][i-1])
            d_cn.append(d_val)

    u_inner_crank_nicolson = thomas_algorithm(a_cn, b_cn, c_cn, d_cn)
    u_crank_nicolson = [0] + u_inner_crank_nicolson + [0]
    tider_crank_nicolson.append(u_crank_nicolson)

    # Euler Explicit
    u_inner_euler_exp = []
    for i in range(1, nx):
        u_inner_euler_exp.append(tider_euler_exp[t-1][i] + r * (tider_euler_exp[t-1][i+1] - 2*tider_euler_exp[t-1][i] + tider_euler_exp[t-1][i-1]))
    u_euler_exp = [0] + u_inner_euler_exp + [0]
    tider_euler_exp.append(u_euler_exp)

    # Beregn analytisk løsning på samme tid
    u_analytisk = [u_exact(i*h, t*k) for i in range(nx + 1)]
    tider_analytisk.append(u_analytisk)

# Animasjon
fig, ax = plt.subplots()
x = [i*h for i in range(nx + 1)]
line_euler_imp, = ax.plot(x, tider_euler_imp[0], label="Euler Implicit")
line_crank_nicolson, = ax.plot(x, tider_crank_nicolson[0], label="Crank-Nicolson")
line_euler_exp, = ax.plot(x, tider_euler_exp[0], label="Euler Explicit")
line_analytisk, = ax.plot(x, tider_analytisk[0], label="Analytisk Løsning", linestyle="--")
ax.set_ylim(-1, 1)
ax.legend()

def update(frame):
    line_euler_imp.set_ydata(tider_euler_imp[frame])
    line_crank_nicolson.set_ydata(tider_crank_nicolson[frame])
    line_euler_exp.set_ydata(tider_euler_exp[frame])
    line_analytisk.set_ydata(tider_analytisk[frame])
    ax.set_title(f"Tid = {frame * k:.4f}")
    return line_euler_imp, line_crank_nicolson, line_euler_exp, line_analytisk

ani = FuncAnimation(fig, update, frames=nt, interval=30, blit=True)
plt.show()
