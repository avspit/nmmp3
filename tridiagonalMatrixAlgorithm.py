import numpy as np
import utils.logger as logger

'''
Расчет СЛАУ
'''
def solve(n, a, b, alpha, beta):
    logger.log(text="Начинаем расчет СЛАУ", force=True)

    zeta = np.zeros(n + 1)
    z = np.zeros(n + 1)
    c = np.zeros(n + 1)
    zeta[0] = beta[0] / alpha[0]
    z[0] = b[0] / a[0]
    for i in range(1, n - 1):
        a[i] = alpha[i] - beta[i - 1] * zeta[i - 1]
        zeta[i] = beta[i] / a[i]
        z[i] = (b[i] - beta[i - 1] * z[i - 1]) / a[i]
    a[n - 1] = alpha[n - 1] - beta[n - 2] * zeta[n - 2]
    z[n - 1] = (b[n - 1] - beta[n - 2] * z[n - 2]) / a[n - 1]

    # переворачиваем список
    indices = []
    for i in range(n):
        a = n - 1 - i
        indices.append(a)
    c[n] = z[n]

    for i in indices:
        c[i] = z[i] - zeta[i] * c[i + 1]

    logger.log(text="Расчет СЛАУ завершен", force=True)

    return c