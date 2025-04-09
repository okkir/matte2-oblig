import math as m
def estdf(x,h):
    return (m.exp(x-2*h)-8*m.exp(x-h)+8*m.exp(x+h)-m.exp(x+2*h))/(12*h)

def df(x):
    return m.exp(x)

h = [10**-i for i in range(0,10)]

print(df(1.5))
for i in h:
    print(i, "estimert derivert: ", estdf(1.5, i), "forskjell", estdf(1.5, i) -df(1.5))


