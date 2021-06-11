import numpy as np
import utils.logger as logger

'''
Расчет СЛАУ
'''
def solve(N, a, b, alpha, beta):
    logger.log(text="Начинаем расчет СЛАУ", force=True)

    zeta = np.zeros(N + 1)
    z = np.zeros(N + 1)
    c = np.zeros(N + 1)
    zeta[0] = beta[0] / alpha[0]
    z[0] = b[0] / a[0]
    for i in range(1, N - 1):
        a[i] = alpha[i] - beta[i - 1] * zeta[i - 1]
        zeta[i] = beta[i] / a[i]
        z[i] = (b[i] - beta[i - 1] * z[i - 1]) / a[i]
    a[N - 1] = alpha[N - 1] - beta[N - 2] * zeta[N - 2]
    z[N - 1] = (b[N - 1] - beta[N - 2] * z[N - 2]) / a[N - 1]

    # переворачиваем список
    indices = []
    for i in range(N):
        a = N - 1 - i
        indices.append(a)
    c[N] = z[N]

    for i in indices:
        c[i] = z[i] - zeta[i] * c[i + 1]

    logger.log(text="Расчет СЛАУ завершен", force=True)

    return c