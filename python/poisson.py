import numpy as np
from helpers import create_laplacian, fftw_r2hc, fftw_hc2r, tridag

class Poisson2D:
    def __init__(self, nx, ny, dx, dyf, dyc):
        self.nx = nx; self.ny = ny
        # DFT laplacian operator in x
        laplacian_x = 2*(np.cos(2*np.pi*np.arange(nx)/nx) - 1)/dx**2
        # construct tridiagonal matrix
        a, b, c = create_laplacian(ny, dyf, dyc)

        self.laplacian_x = laplacian_x
        self.a = a; self.b = b; self.c = c


        self.work = np.zeros((ny, nx))

    def solve(self, phi):
        self.work[:, :] = phi[1:self.ny + 1, 1:self.nx + 1]
        

        # forward transform
        for j in range(self.ny):
            self.work[j, :] = fftw_r2hc(self.work[j, :])

        self.work.tofile('demo_2d_python.bin')

        # solve linear system
        for i in range(self.nx):
            bb = self.b + self.laplacian_x[i]
            tridag(self.a, bb, self.c, self.work[:, i], self.ny)

        
        # transform back
        for j in range(self.ny):
            self.work[j, :] = fftw_hc2r(self.work[j, :])

        phi[1:self.ny + 1, 1:self.nx + 1] = self.work[:, :]/self.nx


        
        # BC (periodic in x and Neumann in y)
        phi[1:-1, 0] = phi[1:-1, -2]
        phi[1:-1, -1] = phi[1:-1, 1]
        phi[0, :] = phi[1, :]
        phi[-1, :] = phi[-2, :]

        return phi
    
class Poisson3D:
    def __init__(self, nx, ny, nz, dx, dy, dzf, dzc):
        self.nx = nx; self.ny = ny; self.nz = nz

        # DFT laplacian operator in x and y
        laplacian_x = 2*(np.cos(2*np.pi*np.arange(nx)/nx) - 1)/dx**2
        laplacian_y = 2 * (np.cos(2 * np.pi * np.arange(ny)/ny) - 1)/dy**2

        # construct tridiagonal matrix
        a, b, c = create_laplacian(ny, dzf, dzc)

        self.laplacian_x = laplacian_x
        self.laplacian_y = laplacian_y
        self.a = a; self.b = b; self.c = c

        self.work = np.zeros((nz, ny, nx))

    def solve(self, phi):
        self.work[:, :, :] = phi[1:self.nz + 1, 1:self.ny + 1, 1:self.nx + 1]

        self.work.tofile('scratch/python/work_before.bin')

        # forward transform
        # x transform
        for k in range(self.nz):
            for j in range(self.ny):
                self.work[k, j, :] = fftw_r2hc(self.work[k, j, :])
        

        # y transform
        for k in range(self.nz):
            for i in range(self.nx):
                self.work[k, :, i] = fftw_r2hc(self.work[k, :, i])

        self.work.tofile('scratch/python/work_after.bin')

        # solve linear system
        for j in range(self.ny):
            for i in range(self.nx):
                bb = self.b + self.laplacian_x[i] + self.laplacian_y[j]
                tridag(self.a, bb, self.c, self.work[:, j, i], self.nz)


        # backward transform
        # y transform
        for k in range(self.nz):
            for i in range(self.nx):
                self.work[k, :, i] = fftw_hc2r(self.work[k, :, i])


        # x transform
        for k in range(self.nz):
            for j in range(self.ny):
                self.work[k, j, :] = fftw_hc2r(self.work[k, j, :])

        phi[1:self.nz + 1, 1:self.ny + 1, 1:self.nx + 1] = self.work[:, :, :]/self.nx/self.ny

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

        return phi
