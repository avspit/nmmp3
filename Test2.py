import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd


'''I use the piecewise linear Rayleigh-Ritz method to solve a second order  ordinary differential equation with a nonzero RHS
The d.e is of the form :  - d/dx(  p(x)dy/dx) + q(x)y = f(x).
The integration is performed on the interval 0<=x<=1'''


def rhs(x):
    '''This computes the rhs of the d.e.
    INPUT: the variable x'''
    return 3*x + np.cos(x**2)


N = 10
a = 0
b = 1
x = np.linspace(a, b, N+1)
h = [(b-a)/N]*(N+1)

def basis_function(x):
    '''Basis functions of the Rayleigh-Ritz method'''

    phi = np.zeros([N+1, N+1])  # initialize basis function
    for j in range(N):
        for idx in range(N):
            if (0 <= x[idx] <= x[j-1]):
                phi[j, idx] = 0
            elif (x[j-1] < x[idx] <= x[j]):
                phi[j, idx] = (x[idx] - x[j-1]) / h[j-1]
            elif (x[j] < x[idx] <= x[j+1]):
                phi[j, idx] = ((x[j+1] - x[idx]) / h[j])
            elif (x[j+1] < x[idx] <= 1):
                phi[j, idx] = 0
    return phi


# basis_function(x)

# coefficient functions
def p(x):
    '''Defining the coefficient p in the d.e.'''
    return -1


def q(x):
    '''Defining the coefficient q in the d.e.'''
    return 0


# initialize vectors for computing integrals
Q1 = []
Q2 = []
Q3 = []
Q4 = []
Q5 = []
Q6 = []
# approximating the 6 integrals
for i in range(N-1):
    q1 = (h[i] / 12) * (q(x[i]) + q(x[i+1]))
    Q1.append(q1)
    q2 = (h[i-1] / 12) * (3 * q(x[i]) + q(x[i-1]))
    Q2.append(q2)
    q3 = (h[i] / 12) * (3 * q(x[i]) + q(x[i+1]))
    Q3.append(q3)
    q4 = (1 / (2 * h[i-1])) * (p(x[i]) + p(x[i-1]))
    Q4.append(q4)
    q5 = (h[i-1] / 6) * (2 * rhs(x[i]) + rhs(x[i-1]))
    Q5.append(q5)
    q6 = (h[i] / 6) * (2 * rhs(x[i]) + rhs(x[i+1]))
    Q6.append(q6)

# now compute  Q1,n , Q2,n , Q3,n , Q4,n Q4,n+1, Q5n, Q6,n

q1n = (h[N]/12)*(q(x[N-1]) + q(x[N]))
Q1.append(q1n)
q2n = (h[N-2] / 12) * (3 * q(x[N-1]) + q(x[N-2]))
Q2.append(q2n)
q3n = (h[N] / 12) * (3 * q(x[N-1]) + q(x[N]))
Q3.append(q3n)
q4n_last = (1 / (2 * h[N-1]) ) * (p(x[N]) + p(x[N-1]))
q4n_second_last = (1 / (2 * h[N-2]) ) * (p(x[N-1]) + p(x[N-2]))
Q4.append(q4n_second_last)
Q4.append(q4n_last)
q5n = (h[N-1] / 6) * (2 * rhs(x[N-1]) + rhs(x[N-2]))
Q5.append(q5n)
q6n = (h[N] / 6) * (2 * rhs(x[N-1]) + rhs(x[N]))
Q6.append(q6n)

alpha = np.zeros(N+1)
beta = np.zeros(N+1)
b = np.zeros(N+1)
a = np.zeros(N+1)
for i in range(N-1):
    alpha[i] = Q4[i] + Q4[i+1] + Q2[i] + Q3[i]
    beta[i] = Q1[i] - Q4[i+1]
    b[i] = Q5[i] + Q6[i]
alpha[N-1] = Q4[N-1] + Q4[N] + Q2[N-1] + Q3[N-1]
b[N-1] = Q5[N-1] + Q6[N-1]
a[0] = alpha[0]

# solve a symmetric tridiagonal linear system
zeta = np.zeros(N+1)
z = np.zeros(N+1)
c = np.zeros(N+1)  # the weights
zeta[0] = beta[0] / alpha[0]
z[0] = b[0] / a[0]
for i in range(1, N-1):
    a[i] = alpha[i] - beta[i-1] * zeta[i-1]
    zeta[i] = beta[i] / a[i]
    z[i] = (b[i] - beta[i-1] * z[i-1]) / a[i]
a[N-1] = alpha[N-1] - beta[N-2] * zeta[N-2]
z[N-1] = (b[N-1] - beta[N-2] * z[N-2]) / a[N-1]

# indices for the weights
# indices must be reversed
indices = []
for i in range(N):
    a = N - 1 - i
    indices.append(a)
# print(indices)
c[N] = z[N]

for i in indices:
    c[i] = z[i] - zeta[i] * c[i+1]
u = basis_function(x)
# plt.plot(u)
# plt.show()

# compute the error
phi_new = np.dot(u, c)

df = pd.DataFrame(c, phi_new)
print(df)

# plotting
plt.grid()
plt.plot(x, phi_new, 'r')
plt.legend(['Rayleigh-Ritz'])
plt.xlabel(r'x')
plt.ylabel(r'$\phi(x)$')
plt.savefig('Results')
plt.show()

