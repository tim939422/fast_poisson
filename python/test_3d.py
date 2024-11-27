import numpy as np
from helpers import *
from poisson import Poisson3D
from derivative import gradpx_3d, gradpy_3d, gradpz_3d
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

    # source
    phi = np.zeros((nz + 2, ny + 2, nx + 2))
    for k in range(1, nz + 1):
        for j in range(1, ny + 1):
            for i in range(1, nx + 1):
                phi[k, j, i] = -(2 + 0.25*np.pi**2)*np.sin(xc[i])*np.sin(yc[j])*np.cos(0.5*np.pi*zc[k])


    # solve pressure Poisson equation
    potential_solver = Poisson3D(nx, ny, ny, dx, dy, dzf, dzc)
    phi = potential_solver.solve(phi)
    
    # Verification

    sol = np.zeros((nz + 2, ny + 2, nx + 2))
    sol = gradpx_3d(nx, ny, nz, dx, phi, sol)
    ref = np.zeros((nz + 2, ny + 2, nx + 2))
    for k in range(1, nz + 1):
        for j in range(1, ny + 1):
            for i in range(nx + 1):
                ref[k, j, i] = np.cos(xf[i])*np.sin(yc[j])*np.cos(0.5*np.pi*zc[k])
    relative_error = np.linalg.norm(ref[1:nz+1, 1:ny+1, :nx + 1] - sol[1:nz+1, 1:ny+1, :nx + 1])/np.linalg.norm(ref[1:nz+1, 1:ny+1, :nx + 1])
    print(f'Relative error in dphi/dx {relative_error:.15e}')

    
    sol = np.zeros((nz + 2, ny + 2, nx + 2))
    sol = gradpy_3d(nx, ny, nz, dy, phi, sol)
    ref = np.zeros((nz + 2, ny + 2, nx + 2))
    for k in range(1, nz + 1):
        for j in range(ny + 1):
            for i in range(1, nx + 1):
                ref[k, j, i] = np.sin(xc[i])*np.cos(yf[j])*np.cos(0.5*np.pi*zc[k])
    relative_error = np.linalg.norm(ref[1:nz+1, :ny+1, 1:nx + 1] - sol[1:nz+1, :ny+1, 1:nx + 1])/np.linalg.norm(ref[1:nz+1, :ny+1, 1:nx + 1])
    print(f'Relative error in dphi/dy {relative_error:.15e}')


    sol = np.zeros((nz + 2, ny + 2, nx + 2))
    sol = gradpz_3d(nx, ny, nz, dzf, phi, sol)
    ref = np.zeros((nz + 2, ny + 2, nx + 2))
    for k in range(nz + 1):
        for j in range(1, ny + 1):
            for i in range(1, nx + 1):
                ref[k, j, i] = -0.5*np.pi*np.sin(xc[i])*np.sin(yc[j])*np.sin(0.5*np.pi*zf[k])
    relative_error = np.linalg.norm(ref[:nz+1, 1:ny+1, 1:nx + 1] - sol[:nz+1, 1:ny+1, 1:nx + 1])/np.linalg.norm(ref[:nz+1, 1:ny+1, 1:nx + 1])
    print(f'Relative error in dphi/dz {relative_error:.15e}')


