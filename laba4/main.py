import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


left_border = 0
right_border = 10
t_min = 0
t_max = 2
a = 5
C = 0.5


def a(x):
    # return x
    return 5*1/(x+1) - 0.7*x
def f(x):
    # return (x**2)/2
    return 5*np.log(x+1) - 0.7*x**2/2

def gamma_r(x, t):
    return 0

def gamma_l(x, t):
    return 0

def fi(x):
    # return np.heaviside(x-1, 1)* np.heaviside(2-x, 1)
    # return np.heaviside(x, 1)* np.heaviside(2-x, 1)*np.sin(np.pi*x/2)**2
    return np.heaviside(4-x, 1) * np.cos(np.pi*x)**2

def method(prev, x,  h, tau):

    next = np.zeros(len(x))

    for i in range(1, len(x) - 1):
        next[i] = 0.5*(prev[i+1] + prev[i-1]) - tau/(2*h)*(f(prev[i+1]) - f(prev[i-1]))

    next[0] = gamma_l(x, tau)
    next[-1] = gamma_r(x, tau)
    return next

def anim():

    m = 200
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
    anim.save('sol_3_1.gif', writer='imagemagick')

anim()
