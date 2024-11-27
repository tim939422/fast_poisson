import numpy as np

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


def twoside_stretch_grid(N, L, beta):
    xf = np.zeros(N + 2)
    xc = np.zeros(N + 2)
    dxc = np.zeros(N + 2) # cell width for pressure \Delta x
    dxf = np.zeros(N + 2) # cell width for u \widetilde{\Delta} x
    for i in range(N + 1):
        xf[i] = L*0.5*(1 - np.tanh(beta*(1 - 2*i/N))/np.tanh(beta))
    xf[N + 1] = 2*xf[N] - xf[N - 1]
    for i in range(1, N + 2):
        xc[i] = 0.5*(xf[i - 1] + xf[i])
    xc[0] = 2*xf[0] - xc[1]

    # metric
    # dxc
    for i in range(1, N + 2):
        dxc[i] = xf[i] - xf[i - 1]
    
    dxc[0] = dxc[1]
    
    # dxf
    for i in range(N + 1):
        dxf[i] = xc[i + 1] - xc[i]
    dxf[N + 1] = dxf[N]
    
    return xf,xc,dxf, dxc

def create_laplacian(N, dxf, dxc):
    a = np.zeros(N + 2)
    b = np.zeros(N + 2)
    c = np.zeros(N + 2)
    for i in range(1, N + 1):
        a[i] = 1/(dxf[i - 1]*dxc[i])
        c[i] = 1/(dxf[i]*dxc[i])
        b[i] = -(a[i] + c[i])

    # Neumann
    b[1] = b[1] + a[1]
    b[N] = b[N] + c[N]
    return a, b, c
