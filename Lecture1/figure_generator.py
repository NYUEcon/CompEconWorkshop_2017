import interpolation.smolyak as smolyak
import matplotlib.pyplot as plt
import numpy as np

from warnings import warn

#
# Curse of dimensionality pictures
#

# Generate ranges...
d = np.arange(1, 7)  # 1, 2, 3, 4, 5, 6
n = 5

# Compute number of points needed for tensor grid
npoints = n**d

fig, ax = plt.subplots()
ax.semilogy(d, npoints)

ax.set_xlabel("Number of dimensions")
ax.set_ylabel("Total Points")
ax.set_title("Curse of dimensionality ({} points in each dimension)".format(n))

fig.savefig("./images/curseofdimensionality.pgf")


sg = smolyak.SmolyakGrid(2, 3)

npoints_smolyak = []
for _d in d[1:]:
    _sg = smolyak.SmolyakGrid(_d, 3)
    npoints_smolyak.append(_sg.grid.shape[0])

fig = plt.figure()

ax1 = fig.add_subplot(2, 2, 1)
ax2 = fig.add_subplot(2, 2, 2)
ax3 = fig.add_subplot(2, 1, 2)

ax1.scatter(sg.grid[:, 0], sg.grid[:, 1])
ax1.set_title("Smolyak Grid in 2d")
ax1.set_xlabel("x")
ax1.set_ylabel("y")

n_rand_pts = 25
rand_grid = np.random.rand(n_rand_pts, 2)
ax2.scatter(rand_grid[:, 0], rand_grid[:, 1])
ax2.set_title("Random Grid in 2d")
ax2.set_xlabel("x")
ax2.set_ylabel("y")

ax3.semilogy(d, npoints, label="Tensor")
ax3.semilogy(d[1:], npoints_smolyak, label="Smolyak")
ax3.semilogy(d, np.ones(d.size)*n_rand_pts, label="Random Grid")
ax3.set_title("Points Comparison")
ax3.set_xlabel("Number of dimensions")
ax3.set_ylabel("Total Points")
ax3.legend()

fig.tight_layout()

fig.savefig("./images/smolyakrandomtensor.pgf")

#
# Root Finding and Optimization Pictures
#
def bisection(f, a, b, ftol=1e-4, xtol=1e-5):
    """
    Uses bisection method to find zeros of f

    Parameters
    ----------
    f : Function
        Function that we want to find the zero of
    a, b : scalar(Float64)
        The endpoints that we want to bisect over... Must be the case
        that f(a)*f(b) < 0
    """
    fa, fb = f(a), f(b)
    assert fa*fb < 0

    # Iterate until done
    dist = 10.0
    xtrace = []  # Keep track of each c
    while dist > xtol:
        # Compute midpoint
        c = (a + b)/2
        fc = f(c)

        if fa*fc < 0:
            b = c
            fb = fc
        else:
            a = c
            fa = fc

        dist = abs(b - a)
        xtrace.append(c)

    # Check whether fc is close to zero
    warn("f(c) is not close to zero") if abs(fc) > ftol else None

    return c, fc, xtrace


def newton(f, df, x0, ftol=1e-4, xtol=1e-5):
    """
    Uses Newton's method to find zeros of f

    Parameters
    ----------
    f : Function
        Function that we want to find the zero of
    df : Function
        Derivative of function that we want to find the zero of
    x0 : scalar(Float64)
        Initial guess of the zero
    """
    # Initialize variables
    xk = x0
    fk = f(x0)
    dfk = df(x0)

    # Iterate till done
    dist = 10.0
    xtrace = []  # Keep track of each xk
    while dist > xtol:
        # Compute x_{k+1}
        xkp1 = xk - fk/dfk
        dist = abs(xkp1 - xk)
        xtrace.append(xkp1)

        # Update all values
        xk = xkp1
        fk = f(xkp1)
        dfk = df(xkp1)

    # Check whether fc is close to zero
    warn("f(x_k) is not close to zero") if abs(fk) > ftol else None

    return xk, fk, xtrace


def secant(f, x0, xm1, ftol=1e-4, xtol=1e-5):
    """
    Uses secant method to find zeros of f

    Parameters
    ----------
    f : Function
        Function that we want to find the zero of
    x0, xm1 : scalar(Float64)
        Initial guesses of the zero
    """
    xk, xkm1 = x0, xm1
    fk, fkm1 = f(xk), f(xkm1)

    # Iterate till done
    dist = 10.0
    xtrace = []  # Keep track of each xk
    while dist > xtol:
        # Compute x_{k+1}
        xkp1 = xk - fk*(xk - xkm1)/(fk - fkm1)
        dist = abs(xkp1 - xk)
        xtrace.append(xkp1)

        # Update all values
        xkm1, fkm1 = xk, fk
        xk = xkp1
        fk = f(xkp1)

    # Check whether fc is close to zero
    warn("f(x_k) is not close to zero") if abs(fk) > ftol else None

    return xk, fk, xtrace


f = lambda x: np.sin(x)
df = lambda x: np.cos(x)

xb, fb, xtraceb = bisection(f, -0.1, 0.5)
xn, fn, xtracen = newton(f, df, -0.1)
xs, fs, xtraces = secant(f, -0.1, -0.09)

