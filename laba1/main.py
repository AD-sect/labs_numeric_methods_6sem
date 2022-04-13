import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def fi1(x):
    return (1/4)*np.exp(x)

def fi2(x):
    return (1/4)*np.exp(x)

def func(x, t):
    return (1/8)*np.exp(x+t)

def coefficients():
    alpha = [-1, -1]
    beta = [2, 3]
    return alpha, beta

def gamma(t):
    gamma1 = (1/4)*np.exp(t)
    gamma2 = (1/2)*np.exp(1+t)
    return gamma1, gamma2

def u0(x, t):
    return (1/4)*np.exp(x+t)

def method_1_ord(left_border, right_border,prev, curr, x,  h, t, t_curr):
    alpha, beta = coefficients()
    a = np.sqrt(0.5)

    next = np.zeros(len(x))

    for i in range(1, len(x) - 1):
        next[i] = ((a * t**2) / h**2) * (curr[i+1] - 2*curr[i] + curr[i - 1]) \
                  + t*t*func(x[i], t_curr - t) + 2*curr[i] - prev[i]

    gamma1, gamma2 = gamma(t_curr)
    # next[0] = (gamma1 - next[1] * (alpha[0]/h)) * (h / (beta[0] * h - alpha[0]))
    # next[-1] = (gamma2 + next[-2] * (alpha[1]/h)) * (h / (beta[1] * h + alpha[1]))

    next[0] = (gamma1 - beta[0] * next[1] / h) / (alpha[0] - beta[0] / h)
    next[-1] = (gamma2 + beta[1] * next[-2] / h) / (alpha[1] + beta[1] / h)
    return next

def method_2_ord(left_border, right_border,prev, curr, x, h, t, t_curr):
    alpha, beta = coefficients()
    a = np.sqrt(0.5)

    next = np.zeros(len(x))

    for i in range(1, len(x) - 1):
        next[i] = ((a*a * t**2) / h**2) * (curr[i + 1] - 2 * curr[i] + curr[i - 1]) \
                  + t * t * func(x[i], t_curr - t) + 2 * curr[i] - prev[i]

    gamma1, gamma2 = gamma(t_curr)
    next[0] = (gamma1 + (alpha[0]/(2*h))*(-4*next[1] + next[2])) * ((2*h) / (-3*alpha[0] + 2*h*beta[0]))
    next[-1] = (gamma2 + (alpha[1]/(2*h))*(-next[-3] + 4*next[-2] )) * ((2*h) / (3*alpha[1] + 2*h*beta[1]))


    return next

def graphic1():
    left_border = 0
    right_border = 1

    n = 20
    a = np.sqrt(0.5)
    x, h = np.linspace(left_border, right_border, n, retstep=True)
    # print(h)

    t = h * 0.5 / a

    curr = np.zeros(len(x))
    prev = np.zeros(len(x))

    prev = fi1(x)
    curr = fi2(x) * t + fi1(x)

    plt.subplot(2, 2, 1)
    plt.title("Метод первого порядка точности")
    plt.xlabel("x")
    plt.ylabel("u")
    plt.grid()
    plt.plot(x, method_1_ord(left_border, right_border, prev, curr, x, h, t, t * 2), color='r', label="Численное решение")
    plt.plot(x, u0(x, 2*t), color='g', label="Точное")
    plt.legend()

    plt.subplot(2, 2, 3)
    plt.title("Метод первого порядка точности")
    plt.xlabel("x")
    plt.ylabel("|Δu|")
    plt.grid()
    plt.plot(x, abs(u0(x, t*2) - method_1_ord(left_border, right_border, prev, curr, x, h, t, t * 2)), color='r', label="Ошибка для одного слоя")
    plt.legend()

    prev = fi1(x)
    curr = fi2(x) * t + fi1(x) + ((t ** 2) / 2) * (a * a * (1/4)*np.exp(x) + func(x, 0))

    plt.subplot(2, 2, 2)
    plt.title("Метод второго порядка точности")
    plt.xlabel("x")
    plt.ylabel("u")
    plt.grid()
    plt.plot(x, method_2_ord(left_border, right_border, prev, curr, x, h, t, t * 2), color='r', label="Численное решение")
    plt.plot(x, u0(x, t * 2), color='g', label="Точное решение")
    plt.legend()

    plt.subplot(2, 2, 4)
    plt.title("Метод первого порядка точности")
    plt.xlabel("x")
    plt.ylabel("|Δu|")
    plt.grid()
    plt.plot(x, abs(u0(x, t*2) - method_2_ord(left_border, right_border, prev, curr, x, h, t, t * 2)), color='r', label="Ошибка для одного слоя")
    plt.legend()


def graphic_error():
    left_border = 0
    right_border = 1

    n = 20
    a = np.sqrt(0.5)
    x, h = np.linspace(left_border, right_border, n, retstep=True)
    # print(h)
    t_min = 0
    t_max = 2

    t = h * 0.5 / a

    curr = np.zeros(len(x))
    prev = np.zeros(len(x))

    prev = fi1(x)
    curr = fi2(x) * t + fi1(x)
    next = np.zeros(len(x))

    time = np.linspace(t_min, t_max,int((t_max - t_min)//t))

    for i in range(2, len(time)+1):
        next = method_1_ord(left_border, right_border, prev, curr, x, h, t, i*t)
        prev = curr
        curr = next
        print(i)
    print(53*t)

    print("len ", len(time))

    print(len(time)*t)

    plt.subplot(1, 2, 1)
    plt.title("Метод первого порядка точности")
    plt.xlabel("x")
    plt.ylabel("|Δu|")
    plt.grid()
    plt.plot(x, abs(u0(x, len(time)*t) - next), color='r',
             label="Ошибка при t ~ 2")
    plt.legend()

    prev = fi1(x)
    curr = fi2(x) * t + fi1(x) + ((t ** 2) / 2) * (a * a * (1 / 4) * np.exp(x) + func(x, 0))
    for i in range(2, len(time)+1):
        next = method_2_ord(left_border, right_border, prev, curr, x, h, t, i*t)
        prev = curr
        curr = next

    plt.subplot(1, 2, 2)
    plt.title("Метод второго порядка точности")
    plt.xlabel("x")
    plt.ylabel("|Δu|")
    plt.grid()
    plt.plot(x, abs(u0(x, len(time)*t) - next), color='r',
             label="Ошибка при t ~ 2")
    plt.legend()
    plt.show()

def animation(ord):

    left_border = 0
    right_border = 1
    t_min = 0
    t_max = 2

    n = 20
    x_mass, h = np.linspace(left_border, right_border, n, retstep=True)

    C = 0.5
    a = np.sqrt(0.5)
    t = C * h / a
    t_curr = t * 2
    u = np.zeros((3, n))

    u[0] = fi1(x_mass)
    if(ord == 1):
        u[1] = fi2(x_mass) * t + fi1(x_mass)
        u[2] = method_1_ord(left_border, right_border, u[0], u[1], x_mass,  h, t, t_curr)
    if(ord == 2):
        u[1] = fi2(x_mass) * t + fi1(x_mass) + ((t**2)/2)*(a*a*(1/4)*np.exp(x_mass) + func(x_mass, 0))
        u[2] = method_2_ord(left_border, right_border, u[0], u[1], x_mass, h, t, t_curr)

    fig = plt.figure()
    ax = plt.axes(xlim=(left_border, right_border), ylim=(-3, 3))
    line1, = ax.plot([], [], lw=3)
    line2, = ax.plot([], [], lw=3)

    def init():
        line1.set_data([], [])
        line2.set_data([], [])
        return line1, line2,

    def animate(i):
        x = x_mass

        if i in [0, 1, 2]:
            y1 = u[i]
            y2 = u0(x, i * t)
        else:
            time = i * t
            u[0] = u[1]
            u[1] = u[2]
            if(ord == 1):
                u[2] = method_1_ord(left_border, right_border, u[0], u[1], x, h, t, time)
            if (ord == 2):
                u[2] = method_2_ord(left_border, right_border, u[0], u[1], x, h, t, time)
            y1 = u[2]
            y2 = u0(x, time)

        line1.set_data(x, y1)
        line2.set_data(x, y2)
        line1.set_color("green")
        line2.set_color("red")
        return line1, line2,

    anim = FuncAnimation(fig, animate, init_func=init,
                         frames=int((t_max - t_min) // t), interval=20, blit=True)

    anim.save('ord2.gif', writer='imagemagick')

if __name__ == '__main__':
     # graphic1()
     # animation(2)
     graphic_error()


#
