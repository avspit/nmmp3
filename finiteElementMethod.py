import numpy as np
import tridiagonalMatrixAlgorithm as tma
import utils.function as func
import utils.equation as eq


'''
Метод конечных элементов (метод Ритца) с помощью функций-крышек
Вид уравнения :  - d/dx(  p(x)dy/dx) + C * q(x)y = A * f(x)
'''
def solve(A, C, N, h, x):
    # Инициализируем векторы для вычисления интегралов
    Q1 = []
    Q2 = []
    Q3 = []
    Q4 = []
    Q5 = []
    Q6 = []
    # аппроксимация 6 интегралов
    for i in range(N-1):
        q1 = (h[i] / 12) * (eq.q(x[i],C) + eq.q(x[i+1],C))
        Q1.append(q1)
        q2 = (h[i-1] / 12) * (3 * eq.q(x[i],C) + eq.q(x[i-1],C))
        Q2.append(q2)
        q3 = (h[i] / 12) * (3 * eq.q(x[i],C) + eq.q(x[i+1],C))
        Q3.append(q3)
        q4 = (1 / (2 * h[i-1])) * (eq.p(x[i]) + eq.p(x[i-1]))
        Q4.append(q4)
        q5 = (h[i-1] / 6) * (2 * eq.f(x[i],A) + eq.f(x[i-1],A))
        Q5.append(q5)
        q6 = (h[i] / 6) * (2 * eq.f(x[i],A) + eq.f(x[i+1],A))
        Q6.append(q6)

    q1n = (h[N]/12)*(eq.q(x[N-1],C) + eq.q(x[N],C))
    Q1.append(q1n)
    q2n = (h[N-2] / 12) * (3 * eq.q(x[N-1],C) + eq.q(x[N-2],C))
    Q2.append(q2n)
    q3n = (h[N] / 12) * (3 * eq.q(x[N-1],C) + eq.q(x[N],C))
    Q3.append(q3n)
    q4n_last = (1 / (2 * h[N-1]) ) * (eq.p(x[N]) + eq.p(x[N-1]))
    q4n_second_last = (1 / (2 * h[N-2]) ) * (eq.p(x[N-1]) + eq.p(x[N-2]))
    Q4.append(q4n_second_last)
    Q4.append(q4n_last)
    q5n = (h[N-1] / 6) * (2 * eq.f(x[N-1],A) + eq.f(x[N-2],A))
    Q5.append(q5n)
    q6n = (h[N] / 6) * (2 * eq.f(x[N-1],A) + eq.f(x[N],A))
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

    y = np.dot(u, c)

    return y
