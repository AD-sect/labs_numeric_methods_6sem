import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


left_border = 0
right_border = 10
t_min = 0
t_max = 3
a = 5
C = 0.5


def a(x):
    return x
def f(x):
    return (x**2)/2

def gamma_r(x, t):
    return 0

def gamma_l(x, t):
    return 0

def fi(x):
    # return np.heaviside(x-1, 1)* np.heaviside(2-x, 1)
    return np.heaviside(x, 1)* np.heaviside(2-x, 1)*np.sin(np.pi*x/2)**2

def method(prev, x,  h, tau):

    next = np.zeros(len(x))

    for i in range(1, len(x) - 1):
        next[i] = 0.5*(prev[i+1] + prev[i-1]) - tau/(2*h)*(f(prev[i+1]) - f(prev[i-1]))

    next[0] = gamma_l(x, tau)
    next[-1] = gamma_r(x, tau)
    return next

# def graphic():
#     m = 2000
#     x = np.linspace(left_border, right_border, m)
#     h = (right_border - left_border) / m
#     tau = C * h/a
#
#     time = np.linspace(t_min, t_max, int((t_max - t_min) // tau))
#
#     prev = fi(x)
#
#
#     for i in range(1, len(time)):
#         print(i)
#         next = method(prev, x, h, tau)
#         prev = next
#
#
#     plt.subplot(1, 2, 1)
#     plt.title("Задача нелинейной теплопроводности при t = 2")
#     plt.xlabel("x")
#     plt.ylabel("u")
#     plt.grid()
#     plt.plot(x, next, color='r', label="Численное решение")
#     plt.legend()
#
#     plt.subplot(1, 2, 2)
#     plt.title("Задача нелинейной теплопроводности при t = 2")
#     plt.xlabel("x")
#     plt.ylabel("u")
#     plt.grid()
#     plt.plot(x, fi(x), color='r', label="Численное решение")
#     plt.legend()
#
#     plt.show()

def anim():

    m = 1000
    x_arr = np.linspace(left_border, right_border, m)
    h = (right_border - left_border) / m

    # tau = C * h / a

    u = np.zeros((2, m))

    u[0] = fi(x_arr)
    tau = (C * h) / max(a(u[0]))
    u[1] = method(u[0], x_arr, h, tau)

    fig = plt.figure()
    ax = plt.axes(xlim=(left_border, right_border), ylim=(-3, 3))
    line1, = ax.plot([], [], lw=3)

    def init():
        line1.set_data([], [])
        return line1,

    def animate(i):
        x = x_arr
        if i in [0, 1]:
            y = u[i]
        else:
            u[0] = u[1]
            tau = (C * h) / max(a(u[0]))
            u[1] = method(u[0], x_arr, h, tau)
            y = u[1]

        line1.set_data(x, y)
        line1.set_color("red")
        return line1,

    anim = FuncAnimation(fig, animate, init_func=init,
                         frames=int((t_max - t_min) // tau), interval=30, blit=True)

    # anim.save('sol_transferring_rectangular_profile.gif', writer='imagemagick')
    anim.save('sol_transferring_with_discontinuous_solution.gif', writer='imagemagick')

# graphic()

anim()
