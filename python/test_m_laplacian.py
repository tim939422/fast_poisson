import numpy as np
import os
from helpers import create_laplacian, staggered_twoside_stretched_grid, staggered_metric

if __name__ == '__main__':

    data_root = '/Users/duosifan/Desktop/code_development/tutorials/fast_poisson/run/test_m_laplacian'
    # DFT operator in x
    Lx = 4*np.pi
    nx = 128
    dx = Lx/nx
    laplacian = 2*(np.cos(2*np.pi*np.arange(nx)/nx) - 1)/dx**2
    laplacian.tofile(os.path.join(data_root, 'ref/laplacian.bin'))

    # a, b, c in y
    Ly= 2
    ny = 128
    beta = 1.2
    yf, yc = staggered_twoside_stretched_grid(ny, Ly, beta)
    dyf, dyc = staggered_metric(ny, yf, yc)
    dyf.tofile(os.path.join(data_root, 'dyf.bin'))
    dyc.tofile(os.path.join(data_root, 'dyc.bin'))
    a, b, c = create_laplacian(ny, dyf, dyc)
    a.tofile(os.path.join(data_root, 'ref/a.bin'))
    b.tofile(os.path.join(data_root, 'ref/b.bin'))
    c.tofile(os.path.join(data_root, 'ref/c.bin'))
    
