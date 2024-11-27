import numpy as np
#from helpers import *
from helpers import staggered_uniform_grid, staggered_twoside_stretched_grid, staggered_metric
from poisson import Poisson2D
from derivative import gradpx_2d, gradpy_2d
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


    # setup source
    phi = np.zeros((ny + 2, nx + 2))
    for j in range(1, ny + 1):
        for i in range(1, nx + 1):
            phi[j, i] = -(1 + 0.25*np.pi**2)*np.sin(xc[i])*np.cos(0.5*np.pi*yc[j])

    # solve pressure Poisson equation
    potential_solver = Poisson2D(nx, ny, dx, dyf, dyc)
    phi = potential_solver.solve(phi)

    # Verification
    sol = np.zeros((ny + 2, nx + 2))
    sol = gradpx_2d(nx, ny, dx, phi, sol)
    ref = np.zeros((ny + 2, nx + 2))
    for j in range(1, ny + 1):
        for i in range(nx + 1):
            ref[j, i] = np.cos(xf[i])*np.cos(0.5*np.pi*yc[j])
    relative_error = np.linalg.norm(ref[1:ny+1, :nx + 1] - sol[1:ny+1, :nx + 1])/np.linalg.norm(ref[1:ny+1, :nx + 1])
    print(f'Relative error in dphi/dx {relative_error:.15e}')


    sol = np.zeros((ny + 2, nx + 2))
    sol = gradpy_2d(nx, ny, dyf, phi, sol)
    ref = np.zeros((ny + 2, nx + 2))
    for j in range(ny+1):
        for i in range(1, nx + 1):
            ref[j, i] = -0.5*np.pi*np.sin(xc[i])*np.sin(0.5*np.pi*yf[j])
    relative_error = np.linalg.norm(ref[:ny+1, 1:nx + 1] - sol[:ny+1, 1:nx + 1])/np.linalg.norm(ref[:ny+1, 1:nx + 1])
    print(f'Relative error in dphi/dy {relative_error:.15e}')
    