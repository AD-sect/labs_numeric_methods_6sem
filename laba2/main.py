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


def ord1(tau, left_border, right_border, k, t):

    alpha, beta = coefficients()
    n = 20
    a_coef = 1
    x, h = np.linspace(left_border, right_border, n, retstep=True)

    prev = fi(x)
    # next = np.zeros(len(x))

    a = []
    b = []
    c = []
    d = []
    for j in range(1, n-1):
        a.append(-(tau * a_coef / h**2) * k)
        b.append(1 + 2 * (tau * a_coef / h**2) * k)
        c.append(-(tau * a_coef / h**2) * alpha)
        d.append((tau * a_coef / h**2) * (1 - alpha) * prev[j + 1] + (1 - 2 * (tau * a_coef / h**2) * (1 - alpha)) * prev[j]
                 + (tau * a_coef / h**2) * (1 - alpha) * prev[j - 1] + tau * func(j * h, (t - 0.5) * tau))

    d.insert(0, gamma_l((t) * tau))
    d.append(h * gamma_r((t) * tau))

    b.insert(0, -(alpha[0]/h) +beta[0])
    b.append((alpha[1]/h) + beta[1])

    c.append(0)
    c.insert(0, alpha[0]/h)

    a.insert(0, 0)
    a.append(-alpha[1]/h)

    a = np.array(a)
    b = np.array(b)
    c = np.array(c)
    d = np.array(d)

    A, B, x = run_through_method_right(a, b, c, d, x, n=1)

    return run_through_method_reverse(A, B, x, n=1)


def run_through_method_right(a, b, c, f_ar, x, n):

    A = np.zeros(len(x))
    B = np.zeros(len(x))

    A[0] = -c[0] / b[0]
    B[0] = f_ar[0] / b[0]
    for i in range(1, len(x) - 1):
        A[i] = -c[i] / (b[i] + a[i]*A[i - 1])
    for i in range(1, len(x) - 1):
        B[i] = (f_ar[i] - a[i]*B[i - 1]) / (b[i] + a[i]*A[i - 1])

    A[-1] = 0
    B[-1] = (f_ar[-1] - a[-1]*B[-2]) / (b[-1] + a[-1] * A[-2])

    return A, B, x

def run_through_method_reverse(A, B, x, n):

    alpha, beta, gamma = coefficients()
    y = np.zeros(len(x))

    if (n == 1):
        y[-1] = B[-1]
    if(n == 2):
        y[-1] = gamma[1]
    for i in range(len(x) - 2, -1, -1):
        y[i] = B[i] + A[i] * y[i + 1]
    return y
