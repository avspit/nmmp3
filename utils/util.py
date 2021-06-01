import numpy as np
from scipy.sparse import spdiags
import constants.variables as const
from utils import logger


# инициализация иксов
def init_x(h, n):
    x = []
    for i in range(0, n+1, 1):
        x.append(const.A + h * i)
    logger.log(text='иксы', value=x)
    return x


# инициализация h
def init_h(n):
    return (const.B - const.A) / n


# инициализация матрицы А
def init_A(n):
    Avalues = np.array([[-1] * (n-1), [2] * (n-1), [-1] * (n-1)])
    Adiags = np.array([-1, 0, 1])
    A = spdiags(Avalues, Adiags, (n-1), (n-1))
    logger.log(text='Матрица А', value=A.todense())
    return A.toarray()


# проверка на пустую строку
def is_not_blank(s):
    return isinstance(s, str) and bool(s and not s.isspace())