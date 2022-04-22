import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def func(x, t):
    return -2 * (x * np.tanh(x * t) + t**2 * (2 * np.tanh(x * t)**2 - 1)) / np.cosh(x * t)

def u0(x, t):
    return 2 / np.cosh(x * t) - np.exp(x + t - 1) / 4

def gamma_r(t):
    return np.exp(t) / 4 + 2 * (1 + 2 * t * np.tanh(t)) / np.cosh(t)

def gamma_l(t):
    return 2 - np.exp(t - 1) / 4

def fi(x):
    return 2 - np.exp(x - 1) / 4

def coefficients():
    alpha = [0, -2]
    beta = [1, 1]
    return alpha, beta


def function(prev, tau, k, t, x, h):
    alpha, beta = coefficients()
    a_coef = 1

    a = np.zeros(len(x))
    b = np.zeros(len(x))
    c = np.zeros(len(x))
    d = np.zeros(len(x))

    d[0] = gamma_l(t*tau)
    d[-1] = gamma_r(t*tau)
    #
    b[0] = -(alpha[0]/h) + beta[0]
    b[-1] = (alpha[1]/h) + beta[1]

    c[0] = alpha[0]/h

    a[-1] = -alpha[1]/h

    for i in range(1, len(x)-1):
        a[i] = (tau * a_coef**2 / h**2) * k
    for i in range(1, len(x)-1):
        b[i] = -1 - 2 * (tau * a_coef**2 / h**2) * k
    for i in range(1, len(x)-1):
        c[i] = (tau * a_coef**2 / h**2) * k
    for i in range(1, len(x)-1):
        d[i] = (tau * a_coef**2 / h ** 2)*(k-1)*(prev[i+1] - 2*prev[i] + prev[i-1])\
                                    -tau * func(x[i], (t - 0.5)*tau) - prev[i]

    A, B = run_through_method_right(a, b, c, d, x)

    return run_through_method_reverse(A, B, x, t)

def run_through_method_right(a, b, c, d, x):

    A = np.zeros(len(x))
    B = np.zeros(len(x))

    A[0] = -c[0] / b[0]
    B[0] = d[0] / b[0]
    for i in range(1, len(x)):
        A[i] = -c[i] / (b[i] + a[i]*A[i - 1])
    for i in range(1, len(x)):
        B[i] = (d[i] - a[i]*B[i - 1]) / (b[i] + a[i]*A[i - 1])

    return A, B

def run_through_method_reverse(A, B, x, t):

    y = np.zeros(len(x))
    gamma1, gamma2 = gamma_l(t), gamma_r(t)

    y[-1] = B[-1]
    for i in range(len(x) - 2, -1, -1):
        y[i] = B[i] + A[i] * y[i + 1]
    return y

def graphic():
    left_border = 0
    right_border = 1
    n = 20
    a = 1
    x = np.linspace(left_border, right_border, n + 1)
    h = (right_border - left_border) / n
    tau = h * 0.5
    prev = fi(x)

    t_min = 0
    t_max = 2

    next = np.zeros(len(x))
    time = np.linspace(t_min, t_max, int((t_max - t_min) // tau))

    for i in range(1, len(time) + 1):
        next = function(prev, tau, 0.5, i, x, h)
        prev = next

    plt.subplot(1, 2, 1)
    plt.title("Метод прогонки для уравнения теплопроводности при t = 2")
    plt.xlabel("x")
    plt.ylabel("u")
    plt.grid()
    plt.plot(x, next, color='r', label="Численное решение")
    plt.plot(x, u0(x, time[-1]), color='g', label="Точное")
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.title("Ошибка на последнем слоe")
    plt.xlabel("x")
    plt.ylabel("|Δu|")
    plt.grid()
    plt.plot(x, abs(u0(x, time[-1]) - next), color='r', label="Ошибка при t = 2")
    plt.legend()

    plt.show()


def animation():

    left_border = 0
    right_border = 1
    t_min = 0
    t_max = 2

    n = 20
    x_arr = np.linspace(left_border, right_border, n + 1)
    h = (right_border - left_border) / n

    a = 1
    tau = 0.5 * h
    t_curr = tau * 2
    u = np.zeros((2, n + 1))

    u[0] = fi(x_arr)
    u[1] = function(u[0], tau, 0.5, 1, x_arr, h)

    fig = plt.figure()
    ax = plt.axes(xlim=(left_border, right_border), ylim=(-3, 3))
    line1, = ax.plot([], [], lw=3)
    line2, = ax.plot([], [], lw=3)

    def init():
        line1.set_data([], [])
        line2.set_data([], [])
        return line1, line2,

    def animate(i):
        x = x_arr
        if i in [0, 1]:
            y1 = u[i]
            y2 = u0(x, i * tau)
        else:

            u[0] = u[1]
            u[1] = function(u[0], tau, 0.5, i, x_arr, h)
            y1 = u[1]
            y2 = u0(x, i*tau)

        line1.set_data(x, y1)
        line2.set_data(x, y2)
        line1.set_color("green")
        line2.set_color("red")
        return line1,

    anim = FuncAnimation(fig, animate, init_func=init,
                         frames= int((t_max - t_min) // tau), interval=30, blit=True)

    anim.save('animation.gif', writer='imagemagick')

def main():
    graphic()
    # animation()

if __name__ == '__main__':
    main()

