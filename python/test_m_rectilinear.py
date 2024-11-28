import numpy as np
from helpers import staggered_uniform_grid, staggered_twoside_stretched_grid, staggered_metric
import os

if __name__ == '__main__':

    nx = 128; ny = 128
    Lx = 4*np.pi
    Ly = 2
    beta = 1.2
    
    # grid
    xf, xc = staggered_uniform_grid(nx, Lx)
    yf, yc = staggered_twoside_stretched_grid(ny, Ly, beta)

    # metric
    dx = Lx/nx
    dyf, dyc = staggered_metric(ny, yf, yc)

    data_root = '/Users/duosifan/Desktop/code_development/tutorials/fast_poisson/run/test_m_rectilinear/ref'
    xf.tofile(os.path.join(data_root, 'xf.bin'))
    xc.tofile(os.path.join(data_root, 'xc.bin'))
    dx = np.array([dx])
    dx.tofile(os.path.join(data_root, 'dx.bin'))

    yf.tofile(os.path.join(data_root, 'yf.bin'))
    yc.tofile(os.path.join(data_root, 'yc.bin'))
    dyf.tofile(os.path.join(data_root, 'dyf.bin'))
    dyc.tofile(os.path.join(data_root, 'dyc.bin'))
