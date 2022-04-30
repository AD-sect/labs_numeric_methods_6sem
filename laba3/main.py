import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


n = 1
k = 1
u0 = 1.25

left_border = 0
right_border = 3
t_min = 0
t_max = 2
hi = 1


def alpha():
    return 1/(n**2)

def f(x, t):
    if x == 0 and t != 0:
        return 1
    elif t == 0:
        return 0
    else:
        z = x/(hi * u0**n * t**(k*n+1))**0.5
        if z <= alpha():
            return alpha()*n*(alpha() - z)

def u_an(x, t):
    if (x <= (alpha()*(hi * u0**n * t**(k*n+1))**0.5) ):
        return u0 * t**k * (f(x, t)**(1/n))
    else:
        return 0

def l(u):
    return u**n

def gamma_r(t):
    return 0

def gamma_l(t):
    return u0 * t**k

def fi(x):
    return 0

def coefficients():
    alpha = [0, 0]
    beta = [1, 1]
    return alpha, beta

def function(prev, prev_iter,  tau, t, x, h):
    alpha, beta = coefficients()

    a = np.zeros(len(x))
    b = np.zeros(len(x))
    c = np.zeros(len(x))
    d = np.zeros(len(x))

    d[0] = gamma_l(t*tau)
    d[-1] = gamma_r(t*tau)

    b[0] = -(alpha[0]/h) + beta[0]
    b[-1] = (alpha[1]/h) + beta[1]

    c[0] = alpha[0]/h

    a[-1] = -alpha[1]/h


    for i in range(1, len(x)-1):
        a[i] = tau/(2*h**2) * ( l(prev_iter[i]) + l(prev_iter[i - 1]))
        b[i] = -tau/(2*h**2) * (2*l(prev_iter[i]) + l(prev_iter[i - 1]) + l(prev_iter[i + 1])) - 1
        c[i] = tau/(2*h**2) * (l(prev_iter[i]) + l(prev_iter[i + 1]))
        d[i] = -prev[i]

    A, B = run_through_method_right(a, b, c, d, x)

    return run_through_method_reverse(A, B, x)

def run_through_method_right(a, b, c, d, x):

    A = np.zeros(len(x))
    B = np.zeros(len(x))

    A[0] = -c[0] / b[0]
    B[0] = d[0] / b[0]
    for i in range(1, len(x)):
        A[i] = -c[i] / (b[i] + a[i]*A[i - 1])
        B[i] = (d[i] - a[i]*B[i - 1]) / (b[i] + a[i]*A[i - 1])

    return A, B

def run_through_method_reverse(A, B, x):

    y = np.zeros(len(x))

    y[-1] = B[-1]
    for i in range(len(x) - 2, -1, -1):
        y[i] = B[i] + A[i] * y[i+1]
    return y


def graphic():
    m = 60
    tau = 0.01


    x = np.linspace(left_border, right_border, m + 1)
    h = (right_border - left_border) / m
    time = np.linspace(t_min, t_max, int((t_max - t_min) // tau))

    u = np.zeros(len(x))
    for i in range(len(x)):
        u[i] = u_an(x[i], time[-1])

    prev = np.zeros(len(x))
    prev_iter = np.zeros(len(x))
    next = np.zeros(len(x))


    for i in range(1, len(time)):
        for j in range(5):
            prev_iter = function(prev, prev_iter, tau, i, x, h)
        prev = prev_iter

    next = prev_iter

    plt.subplot(1, 2, 1)
    plt.title("Задача нелинейной теплопроводности при t = 2")
    plt.xlabel("x")
    plt.ylabel("u")
    plt.grid()
    plt.plot(x, next, color='r', label="Численное решение")
    plt.plot(x, u, color='g', label="Точное")
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.title("Ошибка на последнем слоe")
    plt.xlabel("x")
    plt.ylabel("|Δu|")
    plt.grid()
    plt.plot(x, abs(u - next), color='r', label="Ошибка при t = 2")
    plt.legend()

    plt.show()

def anim():


    m = 60
    x_arr = np.linspace(left_border, right_border, m + 1)
    h = (right_border - left_border) / m

    tau = 0.01

    u = np.zeros((2, m + 1))

    u[0] = np.zeros(len(x_arr))
    u[1] = np.zeros(len(x_arr))

    for j in range(5):
        u[1] = function(u[0], u[1], tau, 1, x_arr, h)

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
        y2 = np.zeros(len(x))
        if i in [0, 1]:
            y1 = u[i]
            for k in range(len(x)):
                y2[k] = u_an(x[k], i*tau)
        else:
            u[0] = u[1]
            for j in range(5):
                u[1] = function(u[0], u[1], tau, i, x, h)
            y1 = u[1]
            for k in range(len(x)):
                y2[k] = u_an(x[k], i * tau)


        line1.set_data(x, y1)
        line2.set_data(x, y2)
        line1.set_color("green")
        line2.set_color("red")
        return line1, line2

    anim = FuncAnimation(fig, animate, init_func=init,
                         frames= int((t_max - t_min) // tau), interval=50, blit=True)

    anim.save('sol.gif', writer='imagemagick')

def main():
    # graphic()
    anim()
#
if __name__ == '__main__':
    main()
