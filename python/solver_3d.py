import numpy as np
from helpers import *
import sys

if __name__ == '__main__':
    nx = int(sys.argv[1])
    ny = nx
    nz = nx
    
    Lx = 4*np.pi
    Ly = 2*np.pi
    Lz = 2
    beta = 1.2

    # grid
    xf, xc = staggered_uniform_grid(nx, Lx)
    yf, yc = staggered_uniform_grid(ny, Ly)
    zf, zc = staggered_twoside_stretched_grid(nz, Lz, beta)

    # metric
    dx = Lx/nx
    dy = Ly/ny
    dzf, dzc = staggered_metric(nz, zf, zc)

    # laplacian
    a, b, c = create_laplacian(ny, dzf, dzc)
    laplacian_x = 2*(np.cos(2*np.pi*np.arange(nx)/nx) - 1)/dx**2
    laplacian_y = 2 * (np.cos(2 * np.pi * np.arange(ny)/ny) - 1)/dy**2

    # source
    phi = np.zeros((nz + 2, ny + 2, nx + 2))
    for k in range(1, nz + 1):
        for j in range(1, ny + 1):
            for i in range(1, nx + 1):
                phi[k, j, i] = -(2 + 0.25*np.pi**2)*np.sin(xc[i])*np.sin(yc[j])*np.cos(0.5*np.pi*zc[k])



    # solver Poisson
    work = np.zeros((nz, ny, nx))
    work[:, :, :] = phi[1:nz + 1, 1:ny + 1, 1:nx + 1]

    # forward transform
    # x transform
    for k in range(nz):
        for j in range(ny):
            work[k, j, :] = fftw_r2hc(work[k, j, :])
    
    # y transform
    for k in range(nz):
        for i in range(nx):
            work[k, :, i] = fftw_r2hc(work[k, :, i])

    # solve linear system
    for j in range(ny):
        for i in range(nx):
            bb = b + laplacian_x[i] + laplacian_y[j]
            tridag(a, bb, c, work[:, j, i], nz)

    # forward transform
    # y transform
    for k in range(nz):
        for i in range(nx):
            work[k, :, i] = fftw_hc2r(work[k, :, i])


    # x transform
    for k in range(nz):
        for j in range(ny):
            work[k, j, :] = fftw_hc2r(work[k, j, :])

    phi[1:nz + 1, 1:ny + 1, 1:nx + 1] = work[:, :, :]/nx/ny

    # BC
    # periodic in x
    phi[1:-1, 1:-1, 0] = phi[1:-1, 1:-1, -2]
    phi[1:-1, 1:-1, -1] = phi[1:-1, 1:-1, 1]
    # periodic in y
    phi[1:-1, 0, :] = phi[1:-1, -2, :]
    phi[1:-1, -1, :] = phi[1:-1, 1, :]
    # Neumann in z
    phi[0, :, :] = phi[1, :, :]
    phi[-1, :, :] = phi[-2, :, :]

    ref = np.zeros((nz + 2, ny + 2, nx + 2))
    sol = np.zeros((nz + 2, ny + 2, nx + 2))

    for k in range(1, nz + 1):
        for j in range(1, ny + 1):
            for i in range(nx + 1):
                ref[k, j, i] = np.cos(xf[i])*np.sin(yc[j])*np.cos(0.5*np.pi*zc[k])
                sol[k, j, i] = (phi[k, j, i + 1] - phi[k, j, i])/dx

    relative_error = np.linalg.norm(ref[1:nz+1, 1:ny+1, :nx + 1] - sol[1:nz+1, 1:ny+1, :nx + 1])/np.linalg.norm(ref[1:nz+1, 1:ny+1, :nx + 1])
    print(f'Relative error in dphi/dx {relative_error:.15e}')

    
    for k in range(1, nz + 1):
        for j in range(ny + 1):
            for i in range(1, nx + 1):
                ref[k, j, i] = np.sin(xc[i])*np.cos(yf[j])*np.cos(0.5*np.pi*zc[k])
                sol[k, j, i] = (phi[k, j + 1, i] - phi[k, j, i])/dy

    relative_error = np.linalg.norm(ref[1:nz+1, :ny+1, 1:nx + 1] - sol[1:nz+1, :ny+1, 1:nx + 1])/np.linalg.norm(ref[1:nz+1, :ny+1, 1:nx + 1])
    print(f'Relative error in dphi/dy {relative_error:.15e}')


    for k in range(nz + 1):
        for j in range(1, ny + 1):
            for i in range(1, nx + 1):
                ref[k, j, i] = -0.5*np.pi*np.sin(xc[i])*np.sin(yc[j])*np.sin(0.5*np.pi*zf[k])
                sol[k, j, i] = (phi[k + 1, j, i] - phi[k, j, i])/dzf[k]

    relative_error = np.linalg.norm(ref[:nz+1, 1:ny+1, 1:nx + 1] - sol[:nz+1, 1:ny+1, 1:nx + 1])/np.linalg.norm(ref[:nz+1, 1:ny+1, 1:nx + 1])
    print(f'Relative error in dphi/dz {relative_error:.15e}')


