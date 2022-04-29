import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

n = 2
k = 3
u0 = 10

left_border = 0
right_boredr = 1


def func(x, t):
    return 0

def u_an(x, t):
    if (x <= (2*(n+2)/n)**(1/2) * (u0**n * t**(n*k+1))**(1/2)):
        return u0 * t**k * ((n/2*(n+2)) * (2*(n*2)/n - x**2/(u0**n * t**(n*k+1))) )**(1/n)
    else:
        return 0


def gamma_r(t):
    return 0

def gamma_l(t):
    return u0 * t**k

def fi():
    return 0

def coefficients():
    alpha = [1, 1]
    beta = [0, 0]
    return alpha, beta
