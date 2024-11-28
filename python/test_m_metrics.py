from helpers import staggered_twoside_stretched_grid, staggered_metric
import numpy as np
import os

if __name__ == "__main__":
    data_root = "/Users/duosifan/Desktop/code_development/tutorials/fast_poisson/run/test_m_metrics"

    L = 2
    n = 128
    beta = 1.2
    xf, xc = staggered_twoside_stretched_grid(n, L, beta)
    xf.tofile(os.path.join(data_root, 'xf.bin'))
    xc.tofile(os.path.join(data_root, 'xc.bin'))

    dxf, dxc = staggered_metric(n, xf, xc)

    dxf.tofile(os.path.join(data_root, 'ref/dxf.bin'))
    dxc.tofile(os.path.join(data_root, 'ref/dxc.bin'))
