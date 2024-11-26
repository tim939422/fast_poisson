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
