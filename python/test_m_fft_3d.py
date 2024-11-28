import numpy as np
from helpers import fftw_hc2r, fftw_r2hc
import os

if __name__ == "__main__":
    data_root = '/Users/duosifan/Desktop/code_development/tutorials/fast_poisson/run/test_m_fft_3d'

    nx = 128
    ny = 64
    nz = 128
    data = np.random.random((nz, ny, nx))
    data.tofile(os.path.join(data_root, 'data.bin'))

    work = np.zeros((nz, ny, nx))

    # forward transform in x
    work[:, :, :] = data[:, :, :]
    for k in range(nz):
        for j in range(ny):
            work[k, j, :] = fftw_r2hc(work[k, j, :])
    work.tofile(os.path.join(data_root, 'ref/fwdx.bin'))

    # backward transform in x
    work[:, :, :] = data[:, :, :]
    for k in range(nz):
        for j in range(ny):
            work[k, j, :] = fftw_hc2r(work[k, j, :])
    work.tofile(os.path.join(data_root, 'ref/bwdx.bin'))

    # forward transform in y
    work[:, :, :] = data[:, :, :]
    for k in range(nz):
        for i in range(nx):
            work[k, :, i] = fftw_r2hc(work[k, :, i])
    work.tofile(os.path.join(data_root, 'ref/fwdy.bin'))


    # backward transform in y
    work[:, :, :] = data[:, :, :]
    for k in range(nz):
        for i in range(nx):
            work[k, :, i] = fftw_hc2r(work[k, :, i])
    work.tofile(os.path.join(data_root, 'ref/bwdy.bin'))
