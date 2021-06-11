import numpy as np
import tridiagonalMatrixAlgorithm as tma
import utils.function as func


'''I use the piecewise linear Rayleigh-Ritz method to solve a second order  ordinary differential equation with a nonzero RHS
The d.e is of the form :  - d/dx(  p(x)dy/dx) + q(x)y = f(x).
The integration is performed on the interval 0<=x<=1'''


def rhs(x, A):
    '''This computes the rhs of the d.e.
    INPUT: the variable x'''
    return A


# basis_function(x)

# coefficient functions
def p(x):
    '''Defining the coefficient p in the d.e.'''
    return -1


def q(x, C):
    '''Defining the coefficient q in the d.e.'''
    return C * (3*x + np.cos(x**2))

def solve(A, C, N, h, x):
    # initialize vectors for computing integrals
    Q1 = []
    Q2 = []
    Q3 = []
    Q4 = []
    Q5 = []
    Q6 = []
    # approximating the 6 integrals
    for i in range(N-1):
        q1 = (h[i] / 12) * (q(x[i],C) + q(x[i+1],C))
        Q1.append(q1)
        q2 = (h[i-1] / 12) * (3 * q(x[i],C) + q(x[i-1],C))
        Q2.append(q2)
        q3 = (h[i] / 12) * (3 * q(x[i],C) + q(x[i+1],C))
        Q3.append(q3)
        q4 = (1 / (2 * h[i-1])) * (p(x[i]) + p(x[i-1]))
        Q4.append(q4)
        q5 = (h[i-1] / 6) * (2 * rhs(x[i],A) + rhs(x[i-1],A))
        Q5.append(q5)
        q6 = (h[i] / 6) * (2 * rhs(x[i],A) + rhs(x[i+1],A))
        Q6.append(q6)

    # now compute  Q1,n , Q2,n , Q3,n , Q4,n Q4,n+1, Q5n, Q6,n

    q1n = (h[N]/12)*(q(x[N-1],C) + q(x[N],C))
    Q1.append(q1n)
    q2n = (h[N-2] / 12) * (3 * q(x[N-1],C) + q(x[N-2],C))
    Q2.append(q2n)
    q3n = (h[N] / 12) * (3 * q(x[N-1],C) + q(x[N],C))
    Q3.append(q3n)
    q4n_last = (1 / (2 * h[N-1]) ) * (p(x[N]) + p(x[N-1]))
    q4n_second_last = (1 / (2 * h[N-2]) ) * (p(x[N-1]) + p(x[N-2]))
    Q4.append(q4n_second_last)
    Q4.append(q4n_last)
    q5n = (h[N-1] / 6) * (2 * rhs(x[N-1],A) + rhs(x[N-2],A))
    Q5.append(q5n)
    q6n = (h[N] / 6) * (2 * rhs(x[N-1],A) + rhs(x[N],A))
    Q6.append(q6n)

    alpha = np.zeros(N+1)
    beta = np.zeros(N+1)
    b = np.zeros(N+1)
    a = np.zeros(N+1)
    for i in range(N-1):
        alpha[i] = Q4[i] + Q4[i+1] + Q2[i] + Q3[i]
        beta[i] = - Q4[i+1] + Q1[i]
        b[i] = Q5[i] + Q6[i]
    alpha[N-1] = Q4[N-1] + Q4[N] + Q2[N-1] + Q3[N-1]
    b[N-1] = Q5[N-1] + Q6[N-1]
    a[0] = alpha[0]

    c = tma.solve(N, a, b, alpha, beta)

    u = func.basis(N, h, x)

    # compute the error
    y = np.dot(u, c)

    return y
