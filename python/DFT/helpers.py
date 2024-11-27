import numpy as np

def fftw_r2hc(x):
    n = len(x)
    x_hat = np.fft.rfft(x)
    x_tilde = np.zeros_like(x)
    x_tilde[:n//2+1] = x_hat[:n//2+1].real
    x_tilde[n//2+1:] = x_hat[n//2 - 1:0:-1].imag
    return x_tilde

def fftw_hc2r(x_tilde):
    n = len(x_tilde)
    x_hat = np.zeros(n//2 + 1, dtype=complex)
    x_hat[:n//2+1].real = x_tilde[:n//2+1]
    x_hat[n//2 - 1:0:-1].imag = x_tilde[n//2+1:]

    x = np.fft.irfft(x_hat, norm='forward')

    return x

