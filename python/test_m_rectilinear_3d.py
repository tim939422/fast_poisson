import numpy as np
from helpers import staggered_uniform_grid, staggered_twoside_stretched_grid, staggered_metric
import os

if __name__ == '__main__':

    nx = 128; ny = 128; nz = 128
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

    data_root = '/Users/duosifan/Desktop/code_development/tutorials/fast_poisson/run/test_m_rectilinear_3d/ref'
    xf.tofile(os.path.join(data_root, 'xf.bin'))
    xc.tofile(os.path.join(data_root, 'xc.bin'))
    dx = np.array([dx])
    dx.tofile(os.path.join(data_root, 'dx.bin'))

    yf.tofile(os.path.join(data_root, 'yf.bin'))
    yc.tofile(os.path.join(data_root, 'yc.bin'))
    dy = np.array([dy])
    dy.tofile(os.path.join(data_root, 'dy.bin'))


    zf.tofile(os.path.join(data_root, 'zf.bin'))
    zc.tofile(os.path.join(data_root, 'zc.bin'))
    dzf.tofile(os.path.join(data_root, 'dzf.bin'))
    dzc.tofile(os.path.join(data_root, 'dzc.bin'))
