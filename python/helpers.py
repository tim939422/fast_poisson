import numpy as np

def twoside_stretch_grid(N, H, beta):
    z = np.zeros(N + 1)
    for k in range(N + 1):
        z[k] = H/2*(1 - np.tanh(beta*(1 - 2*k/N))/np.tanh(beta))

    return z
    
def staggered_twoside_stretched_grid(N, H, beta):
    zf = np.zeros(N + 2)
    zc = np.zeros(N + 2)

    # face coordinate
    zf[:N + 1] = twoside_stretch_grid(N, H, beta)
    zf[N + 1] = 2*zf[N] - zf[N - 1]

    # cell coordinate
    for k in range(1, N + 2):
        zc[k] = 0.5*(zf[k - 1] + zf[k])
    zc[0] = 2*zf[0] - zc[1]

    return zf, zc

def staggered_uniform_grid(N, L):
    xf = np.zeros(N + 2); xc = np.zeros(N + 2)

    # face coordinate
    xf[:-1] = np.linspace(0, L, N + 1)
    xf[N + 1] = 2*xf[N] - xf[N - 1]

    # cell coordinate
    xc[1:] = (xf[:-1] + xf[1:])*0.5
    xc[0] = 2*xf[0] - xc[1]

    return xf, xc

def staggered_metric(N, zf, zc):
    dzf = np.zeros(N + 2) # metric dz/dzeta at face coordinate
    dzc = np.zeros(N + 2) # metric dz/dzeta at cell coordinate
    
    for i in range(N + 1):
        dzf[i] = zc[i + 1] - zc[i]
    dzf[N + 1] = dzf[N]

    for i in range(1, N + 2):
        dzc[i] = zf[i] - zf[i - 1]
    dzc[0] = dzc[1]

    return dzf, dzc
    


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


'''
[ b1    c1     0     ...       0     ] [ u1   ]   = [ r1   ]
[ a2    b2    c2     ...       0     ] [ u2   ]     [ r2   ]
[  0    a3    b3    c3     ...       ] [ u3   ]     [ r3   ]
[ ...   ...   ...   ...     ...      ] [ ...  ]     [ ...  ]
[  0     0   aN-1  bN-1   cN-1       ] [ uN-1 ]     [ rN-1 ]
[  0     0     0    aN     bN        ] [ uN   ]     [ rN   ]
'''
# This is a modified version of numerical recipes (p 43)
def tridag(a, b, c, r, n):
    eps = np.finfo(float).eps
    gam = np.zeros(n)

    bet  = b[0]
    r[0] = r[0]/bet

    # Decomposition and forward substitution
    for j in range(1, n):
        gam[j] = c[j - 1]/bet
        bet = b[j] - a[j]*gam[j] + eps # for singular case in N-N BC
        r[j] = (r[j] - a[j]*r[j - 1])/bet

    # Backsubstitution
    for j in range(n - 2, -1, -1):
        r[j] = r[j] - gam[j + 1]*r[j + 1]


def create_laplacian(N, dzf, dzc):
    a = np.zeros(N + 2)
    b = np.zeros(N + 2)
    c = np.zeros(N + 2)
    for i in range(1, N + 1):
        a[i] = 1/(dzf[i - 1]*dzc[i])
        c[i] = 1/(dzf[i]*dzc[i])
        b[i] = -(a[i] + c[i])

    # Neumann BC
    b[1] = b[1] + a[1]
    b[N] = b[N] + c[N]
    
    return a[1:-1], b[1:-1], c[1:-1]