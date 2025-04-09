import math as m
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
tider = []

h = 0.01
k = 0.00005
T = 0.1

def f(x):
    return m.sin(m.pi*x)

t0 = [f(i*h) for i in range(0,int(1/h)+1)]
tider.append(t0)

for j in range(1,int(1/k) -1):
    punkter = []
    for i in range(0,int(1/h)+1):
        if i == 0 or i == int(1/h):
            punkter.append(0)
        else:
            punkter.append(tider[j-1][i]+((tider[j-1][i+1] -2*tider[j-1][i] + tider[j-1][i-1])/h**2)*k )
    tider.append(punkter)
print(tider)

fig, ax = plt.subplots()
x = [i*h for i in range(0, int(1/h)+1)]
line, = ax.plot(x,tider[0])
ax.set_ylim(-1,1)
def update(frame):
    y = tider[frame]
    line.set_ydata(y)
    return line,

ani = FuncAnimation(fig, update, frames=int(1/k)-1, interval=10, blit=True)
plt.show()


