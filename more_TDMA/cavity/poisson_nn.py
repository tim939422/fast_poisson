from scipy.fftpack import dct
import numpy as np
from linear_solver import tridag


def poisson_v1(Nx, Ny, dx, dy, f):
    a = np.ones(Ny)/dy**2
    c = np.ones(Ny)/dy**2
    b = -(a + c)
    b[0]  += a[0]
    b[-1] += c[-1]

    lambda_x = -4*(np.sin(np.arange(Nx)*np.pi/(2*Nx)))**2
    laplacian_x = lambda_x/dx**2
    norm = 1/(2*Nx)

    f_hat_x = dct(f, type=2, axis=1)
    f_hat = np.copy(f_hat_x)
    buffer = np.zeros(Ny, dtype=float)
    for i in range(Nx):
        bb = b + laplacian_x[i]
        buffer[:] = f_hat[:, i]
        tridag(a, bb, c, buffer, Ny)
        f_hat[:, i] = buffer[:]

    p = norm*dct(f_hat, type=3, axis=1)
    p = p
    return p

def poisson_v2(Nx, Ny, dx, dy, f):
    lambda_xx = -4*(np.sin(np.arange(Nx)*np.pi/(2*Nx)))**2
    lambda_yy = -4*(np.sin(np.arange(Ny)*np.pi/(2*Ny)))**2
    lambda_x, lambda_y = np.meshgrid(lambda_xx, lambda_yy)
    laplacian = lambda_x/dx**2 + lambda_y/dy**2
    laplacian = np.where(laplacian != 0.0, laplacian, 1)
    norm = 1/(2*Nx)*1/(2*Ny)
    f_hat = dct(dct(f, type=2, axis=1), type=2, axis=0)
    f_hat[0, 0] = 0.0 # force mean to zero
    p = norm*dct(dct(f_hat/laplacian, type=3, axis=0), type=3, axis=1)
    return p