import math as m
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Parametere
h = 0.001
k = 0.1
L = 1
T = 3
nx = int(L / h)
nt = int(T / k)
r = k / h**2


def f(x):
    return m.sin(m.pi * x)

u = [f(i*h) for i in range(nx + 1)]
tider = [u.copy()]


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


n_indre = nx - 1  
for t in range(1, nt):

    a = [-r] * (n_indre - 1)      
    b = [1 + 2*r] * n_indre       
    c = [-r] * (n_indre - 1)      
    d = u[1:-1]                   

    u_inner = thomas_algorithm(a, b, c, d)
    u = [0] + u_inner + [0]  
    tider.append(u.copy())


fig, ax = plt.subplots()
x = [i*h for i in range(nx + 1)]
line, = ax.plot(x, tider[0])
ax.set_ylim(-1, 1)

def update(frame):
    line.set_ydata(tider[frame])
    ax.set_title(f"Tid = {frame * k:.4f}")
    return line,

ani = FuncAnimation(fig, update, frames=nt, interval=30, blit=True)
plt.show()
