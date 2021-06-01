import numpy as np


# g(x)
# -u'' + A(3u + cos(u)^2) = C   =>   u'' = A(3u + cos(u)^2) - C
def g(A, C, h, i, x):
    return A * (3*x + np.cos(x)) - C


# производная g(x)
def dg(A, x):
    return A * (3 - np.sin(x))