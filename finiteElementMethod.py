import numpy as np
import constants.variables as const
import tridiagonalMatrixAlgorithm as tma
import utils.logger as logger
import utils.util as util
import utils.equation as equation
import utils.function as func

def simpson(f, a, b, n):
    h=(b-a)/n
    k=0.0
    x=a + h
    for i in range(1,int(n/2 + 1)):
        k += 4*equation.g(1,0,0,0,x)
        x += 2*h

    x = a + 2*h
    for i in range(1,int(n/2)):
        k += 2*equation.g(1,0,0,0,x)
        x += 2*h
    return (h/3)*(equation.g(1,0,0,0,a)+equation.g(1,0,0,0,b)+k)

# метод Ритца
def solve(inputA, inputC, n, x):
    y = np.zeros((n + 1, 1))  # результат решения
    a = np.empty(shape=(n+1,n+1)) #np.zeros((n + 1, 1));
    b = np.zeros((n + 1, 1));

    for i in range(0, n+1, 1):
        for j in range(i, min(i+3,n+1), 1):
            if (j-2 < 0):
                L=0
            else:
                L = max(x[j-2],0);
            if (j+2 > len(x)-1):
                len(x)-1
            else:
                U = min(x[j+2],1);

            # ... вычисляем a[i][j]
            a[i][j] = simpson(func.dfi(n, 0.1, x[i], i)*func.dfi(n, 0.1, x[j], j) + func.fi(n, 0.1, x[i], i)*func.fi(n, 0.1, x[j], j), L, U, n)

            if (i != j):
                a[j][i] = a[i][j];

        if (i <= n-3):
            for j in range(i+4, n+1, 1):
                a[i][j] = 0;

        if (i-2 < 0):
            L=0
        else:
            L = max(x[i-2], 0);
        if (i+2 > len(x)-1):
            len(x)-1
        else:
            U = min(x[i+2], 1);

        # ... вычисляем b[i]
        b[i] = simpson(
            equation.g(inputA, inputC, 0.1, i, x[i]) * func.fi(n, 0.1, x[i], i),
            L, U, n)

    # вычисление вектора-решения методом прогонки
    c = tma.solve(a, b);

    for i in range(0, n+1, 1):
        y[i] = c[i]

    return y