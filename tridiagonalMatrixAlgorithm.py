import numpy as np
import utils.logger as logger

'''
Метод прогонки

A - матрица уравнений
B - решения
'''
def solve(A, B):
    logger.log(text="Начинаем расчет СЛАУ методом прогонки", force=True)

    n = B.size

    # 1. прямой ход
    # вычисление альфы и беты для первой строки матрицы
    b0 = A[0][0]
    c0 = A[0][1]
    d0 = B[0]
    y0 = b0
    alfa0 = -c0 / y0
    beta0 = d0 / y0
    alfas = [alfa0]
    betas = [beta0]
    # вычисление альф и бет для оставшихся строк матрицы
    for i in range(1, n, 1):
        a = A[i][i-1]
        b = A[i][i]
        c = A[i][i+1] if i < n-1 else 0
        d = B[i]
        y = b + a * alfas[i-1]
        alfa = -c / y
        beta = (d - a * betas[i-1]) / y
        alfas.append(alfa)
        betas.append(beta)

    # 2. обратный ход
    # вычисление первого х
    x = [betas[betas.__len__()-1]]
    # вычисление оставшихся х
    for i in range(n-2, -1, -1):
        x.append(alfas[i] * x[n-i-2] + betas[i])
    # переворачиваем список х-ов
    x = np.array(x)[::-1]

    logger.log(text="Расчет СЛАУ методом прогонки завершен", force=True)

    return x