import numpy as np
import matplotlib.pyplot as plt


def func(x, t):
    return -2 * (x * np.tanh(x * t) + t**2 * (2 * np.tanh(x * t)**2 - 1)) / np.cosh(x * t)

def u0(x, t):
    return 2 / np.cosh(x * t) - np.exp(x + t - 1) / 4

def gamma_l(t):
    return np.exp(t) / 4 + 2 * (1 + 2 * t * np.tanh(t)) / np.cosh(t)
def gamma_r(t):
    return 2 - np.exp(t - 1) / 4

def fi(x):
    return 2 - np.exp(x - 1) / 4


def coefficients():
    alpha = [0, -2]
    beta = [1, 1]
    return alpha, beta


def ord1(h, NJ, tau, NT, a, alpha, left_border, right_border, k, t):
    alpha, beta = coefficients()
    n = 20
    a_coef = 1
    x, h = np.linspace(left_border, right_border, n, retstep=True)

    prev = fi(x)
    next = np.zeros(len(x))

    a = []
    b = []
    c = []
    d = []
    for j in range(1, NJ-1):
        a.append(-(tau * a_coef / h**2) * k)
        b.append(1 + 2 * (tau * a_coef / h**2) * k)
        c.append(-(tau * a_coef / h**2) * alpha)
        d.append((tau * a_coef / h**2) * (1 - alpha) * prev[j + 1] + (1 - 2 * (tau * a_coef / h**2) * (1 - alpha)) * prev[j]
                 + (tau * a_coef / h**2) * (1 - alpha) * prev[j - 1] + tau * func(j * h, (t - 0.5) * tau))

    d.insert(0, gamma_l((t) * tau))
    d.append(h * gamma_r((t) * tau))
    b.insert(0, -(alpha[0]/h) +beta[0])
    b.append((alpha[0]/h) + beta[0])
    c.append(2)
    C.insert(0, alpha[0]/h)
    A = [0, 0] + A
    A = np.array(A)
    B = np.array(B)
    C = np.array(C)
    F = np.array(F)
    temp = np.vstack([A, B, C])
    ans = scipy.linalg.solve_banded((1,1), temp, F)
    result[i + 1, ::] = ans[::]
    return result
