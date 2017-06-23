import numpy as np
import scipy.linalg as la

from collections import namedtuple
from math import sqrt
from numba import jit

#
# Getting to know your computer
#

# 1. Compute the fp_num
bin_fp_num = "1101011101010000"

sign = int(bin_fp_num[0])
exponent = bin_fp_num[1:6]
mantissa = bin_fp_num[6:]

eval_mantissa = 0.0
for (i_m, m) in enumerate(mantissa):
    eval_mantissa += 2**(-i_m - 1) * int(m)

fp_num = (-1)**sign * (1 + eval_mantissa) * 2**(int(exponent, 2) - 15)
print("The 16 bit binary number is {}".format(fp_num))

# 2. Compute eigen-values and condition number
A = np.array([[5, 3], [0, 5]])
eig_vals, eig_vecs = la.eig(A)
eig_vals.sort()

print("Eigen values are {} and {}".format(eig_vals[0], eig_vals[1]))
print("Condition number is {}".format(eig_vals[1] / eig_vals[0]))

# 3. Summation examples

@jit(nopython=True)
def sum_cols_first(X):
    """
    Sum a matrix with cols as outer loop and rows as
    inner loop
    """
    total = 0.0
    ncol, nrow = X.shape
    for i in range(ncol):
        for j in range(nrow):
            total += X[i, j]

    return total


@jit(nopython=True)
def sum_rows_first(X):
    """
    Sum a matrix with rows as outer loop and cols as
    inner loop
    """
    total = 0.0
    ncol, nrow = X.shape
    for j in range(nrow):
        for i in range(ncol):
            total += X[i, j]

    return total


X = np.random.randn(1000, 1000)

col_first = sum_cols_first(X)
row_first = sum_rows_first(X)

#
# Root finding and optimization
#

def secant(f, x0, xm1, ftol=1e-6, xtol=1e-8):
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
    while dist > xtol:
        # Compute x_{k+1}
        xkp1 = xk - fk*(xk - xkm1)/(fk - fkm1)
        dist = abs(xkp1 - xk)

        # Update all values
        xkm1, fkm1 = xk, fk
        xk = xkp1
        fk = f(xkp1)

    # Check whether fc is close to zero
    warn("f(x_k) is not close to zero") if abs(fk) > ftol else None

    return xk


def golden_search(f, a, b, xtol=1e-8):
    """
    Use golden search to find the maximum of a function f

    Parameters
    ----------
    f : Function
        Function that we want to find the zero of
    a, b : scalar(Float64)
        Bracket that we want to find maximum in
    """
    # Golden ratio
    phi = (1 + sqrt(5))/2

    # Evaluate function at a and b
    fa, fb = f(a), f(b)

    dist = 10.0
    while dist > xtol:
        # Compute new values
        c = b - (b-a)/phi
        d = a + (b-a)/phi

        # Evaluate f at c and d
        fc = f(c)
        fd = f(d)

        if fc > fd:
            b, fb = d, fd
        else:
            a, fa = c, fc

        dist = abs(a - b)

    return (a+b)/2


# Solve the problem
model = namedtuple("model", ["A", "B", "x", "w"])


def BC_A(m, a, p):
    "Uses budget constraint to say how many bananas are bought by A"
    b = (m.A - a)/p

    return b

def BC_B(m, a, p):
    "Uses budget constraint to say how many bananas are bought by A"
    b = (p*m.B - a)/p

    return b


def utility_A(m, a, p):
    """
    Country A's utility function
    """
    # Unpack parameters
    A, B, x, w = m.A, m.B, m.x, m.w

    # Use budget constraint to get b in terms of a
    b = BC_A(m, a, p)

    return (m.w * a**x + (1 - m.w) * b**x)**(1.0/x)


def utility_B(m, a, p):
    """
    Country B's utility function
    """
    # Unpack parameters
    A, B, x, w = m.A, m.B, m.x, m.w

    # Use budget constraint to get b in terms of a
    b = BC_B(m, a, p)

    return (m.w * b**x + (1 - m.w) * a**x)**(1.0/x)


def compute_demand(m, p):
    "Gets demand for apples in both countries given prices"
    # Make sure p is positive
    p = max(p, 1e-2)

    # Compute demand in country A
    opt_A = lambda a: utility_A(m, a, p)
    a_A_star = golden_search(opt_A, 1e-3, m.A-1e-3)

    # Compute demand in country B
    opt_B = lambda a: utility_B(m, a, p)
    a_B_star = golden_search(opt_B, 1e-3, p*m.B-1e-3)

    return a_A_star, a_B_star


def excess_demand(m, p):
    "Computes excess demand"
    a_A_star, a_B_star = compute_demand(m, p)

    return a_A_star + a_B_star - m.A


def solve_model(m):
    """
    Solve model described in practice question document

    Substitute out to get optimization in terms of only apples in both
    maximization problems. Root finding searches over relative price in
    bananas to apples
    """
    # Solve for market clearing price
    solve_me = lambda p: excess_demand(m, p)
    pstar = secant(solve_me, 1.0, 1.5)

    # Pull out corresponding demand for a and b in each country
    a_A, a_B = compute_demand(m, pstar)
    b_A, b_B = BC_A(m, a_A, pstar), BC_B(m, a_B, pstar)

    return pstar, a_A, b_A, a_B, b_B


m1 = model(5.0, 10.0, -0.5, 0.9)
m2 = model(5.0, 5.0, -0.5, 0.9)

p1, a_A1, b_A1, a_B1, b_B1 = solve_model(m1)
p2, a_A2, b_A2, a_B2, b_B2 = solve_model(m2)

