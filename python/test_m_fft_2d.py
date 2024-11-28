import numpy as np
from helpers import fftw_hc2r, fftw_r2hc
import os

if __name__ == "__main__":
    data_root = '/Users/duosifan/Desktop/code_development/tutorials/fast_poisson/run/test_m_fft_2d'

    nx = 256
    ny = 256
    data = np.random.random((ny, nx))
    data.tofile(os.path.join(data_root, 'data.bin'))

    work = np.zeros((ny, nx))

    # forward transform in x
    work[:, :] = data[:, :]
    for j in range(ny):
        work[j, :] = fftw_r2hc(work[j, :])
    work.tofile(os.path.join(data_root, 'ref/fwdx.bin'))

    # backward transform in x
    work[:, :] = data[:, :]
    for j in range(ny):
        work[j, :] = fftw_hc2r(work[j, :])
    work.tofile(os.path.join(data_root, 'ref/bwdx.bin'))

    # forward transform in y
    work[:, :] = data[:, :]
    for i in range(nx):
        work[:, i] = fftw_r2hc(work[:, i])
    work.tofile(os.path.join(data_root, 'ref/fwdy.bin'))


    # backward transform in y
    work[:, :] = data[:, :]
    for i in range(nx):
        work[:, i] = fftw_hc2r(work[:, i])
    work.tofile(os.path.join(data_root, 'ref/bwdy.bin'))
    