import numpy as np
from helpers import *
import sys

if __name__ == '__main__':
    nx = int(sys.argv[1])
    ny = nx
    
    Lx = 4*np.pi
    Ly = 2
    beta = 1.2

    # grid
    xf, xc = staggered_uniform_grid(nx, Lx)
    yf, yc = staggered_twoside_stretched_grid(ny, Ly, beta)

    # metric
    dx = Lx/nx
    dyf, dyc = staggered_metric(ny, yf, yc)

    # laplacian
    a, b, c = create_laplacian(ny, dyf, dyc)
    laplacian_x = 2*(np.cos(2*np.pi*np.arange(nx)/nx) - 1)/dx**2

    # source
    phi = np.zeros((ny + 2, nx + 2))
    for j in range(1, ny + 1):
        for i in range(1, nx + 1):
            phi[j, i] = -(1 + 0.25*np.pi**2)*np.sin(xc[i])*np.cos(0.5*np.pi*yc[j])

    # solve Poisson
    work = np.zeros((ny, nx))
    work[:, :] = phi[1:ny + 1, 1:nx + 1]

    
    # forward transform
    for j in range(ny):
        work[j, :] = fftw_r2hc(work[j, :])
    
    # solve linear system
    for i in range(nx):
        bb = b + laplacian_x[i]
        tridag(a, bb, c, work[:, i], ny)
    
    # transform back
    for j in range(ny):
        work[j, :] = fftw_hc2r(work[j, :])


    phi[1:ny + 1, 1:nx + 1] = work[:, :]/nx

    # BC
    phi[1:-1, 0] = phi[1:-1, -2]
    phi[1:-1, -1] = phi[1:-1, 1]
    phi[0, :] = phi[1, :]
    phi[-1, :] = phi[-2, :]

    ref = np.zeros((ny + 2, nx + 2))
    sol = np.zeros((ny + 2, nx + 2))
    for j in range(1, ny + 1):
        for i in range(nx + 1):
            ref[j, i] = np.cos(xf[i])*np.cos(0.5*np.pi*yc[j])
            sol[j, i] = (phi[j, i + 1] - phi[j, i])/dx

    relative_error = np.linalg.norm(ref[1:ny+1, :nx + 1] - sol[1:ny+1, :nx + 1])/np.linalg.norm(ref[1:ny+1, :nx + 1])
    print(f'Relative error in dphi/dx {relative_error:.15e}')

    for j in range(ny+1):
        for i in range(1, nx + 1):
            ref[j, i] = -0.5*np.pi*np.sin(xc[i])*np.sin(0.5*np.pi*yf[j])
            sol[j, i] = (phi[j + 1, i] - phi[j, i])/dyf[j]

    relative_error = np.linalg.norm(ref[:ny+1, 1:nx + 1] - sol[:ny+1, 1:nx + 1])/np.linalg.norm(ref[:ny+1, 1:nx + 1])
    print(f'Relative error in dphi/dy {relative_error:.15e}')
    