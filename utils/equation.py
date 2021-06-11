import numpy as np


def f(x, A):
    return A * x


def p(x):
    return -1


def q(x, C):
    return C * (3*x + np.cos(x**2))