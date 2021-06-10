import numpy as np

# Базисные функции (функции-крышки)
def basis(N, h, x):
    phi = np.zeros([N+1, N+1])
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