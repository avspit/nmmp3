import numpy as np
import tridiagonalMatrixAlgorithm as tma
import utils.function as func
import utils.equation as eq


'''
Метод конечных элементов с помощью функций-крышек
Вид уравнения :  - d/dx(  p(x)dy/dx) + C * q(x)y = A * f(x)
'''
def solve(A, C, n, h, x):
    # Инициализируем векторы для вычисления интегралов
    Q1 = []
    Q2 = []
    Q3 = []
    Q4 = []
    Q5 = []
    Q6 = []
    # аппроксимация 6 интегралов
    for i in range(n-1):
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

    q1n = (h[n]/12)*(eq.q(x[n-1],C) + eq.q(x[n],C))
    Q1.append(q1n)
    q2n = (h[n-2] / 12) * (3 * eq.q(x[n-1],C) + eq.q(x[n-2],C))
    Q2.append(q2n)
    q3n = (h[n] / 12) * (3 * eq.q(x[n-1],C) + eq.q(x[n],C))
    Q3.append(q3n)
    q4n_last = (1 / (2 * h[n-1])) * (eq.p(x[n]) + eq.p(x[n-1]))
    q4n_second_last = (1 / (2 * h[n-2])) * (eq.p(x[n-1]) + eq.p(x[n-2]))
    Q4.append(q4n_second_last)
    Q4.append(q4n_last)
    q5n = (h[n-1] / 6) * (2 * eq.f(x[n-1],A) + eq.f(x[n-2],A))
    Q5.append(q5n)
    q6n = (h[n] / 6) * (2 * eq.f(x[n-1],A) + eq.f(x[n],A))
    Q6.append(q6n)

    alpha = np.zeros(n+1)
    beta = np.zeros(n+1)
    b = np.zeros(n+1)
    a = np.zeros(n+1)
    for i in range(n-1):
        alpha[i] = Q4[i] + Q4[i+1] + Q2[i] + Q3[i]
        beta[i] = - Q4[i+1] + Q1[i]
        b[i] = Q5[i] + Q6[i]
    alpha[n-1] = Q4[n-1] + Q4[n] + Q2[n-1] + Q3[n-1]
    b[n-1] = Q5[n-1] + Q6[n-1]
    a[0] = alpha[0]

    c = tma.solve(n, a, b, alpha, beta)
    u = func.basis(n, h, x)

    y = np.dot(u, c)

    return y
