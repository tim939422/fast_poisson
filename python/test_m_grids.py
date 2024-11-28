from helpers import staggered_twoside_stretched_grid, staggered_uniform_grid
import numpy as np
import os

if __name__ == "__main__":
    data_root = '/Users/duosifan/Desktop/code_development/tutorials/fast_poisson/run/test_m_grid/ref'
    
    ''' Generate data for uniform grid
    '''
    L = 4*np.pi
    n = 128
    xf, xc = staggered_uniform_grid(n, L)
    xf.tofile(os.path.join(data_root, 'staggered_uniform_grid_xf.bin'))
    xc.tofile(os.path.join(data_root, 'staggered_uniform_grid_xc.bin'))

    ''' Generate data for twoside stretched grid
    '''
    L = 2
    n = 128
    beta = 1.2
    xf, xc = staggered_twoside_stretched_grid(n, L, beta)
    xf.tofile(os.path.join(data_root, 'staggered_twoside_stretched_grid_xf.bin'))
    xc.tofile(os.path.join(data_root, 'staggered_twoside_stretched_grid_xc.bin'))

    
