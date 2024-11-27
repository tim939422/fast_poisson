def gradpx_2d(nx, ny, dx, phi, work):
    for j in range(1, ny + 1):
        for i in range(nx + 1):
            work[j, i] = (phi[j, i + 1] - phi[j, i])/dx

    return work

def gradpy_2d(nx, ny, dyf, phi, work):
    for j in range(ny+1):
        for i in range(1, nx + 1):
            work[j, i] = (phi[j + 1, i] - phi[j, i])/dyf[j]

    return work

def gradpx_3d(nx, ny, nz, dx, phi, work):
    for k in range(1, nz + 1):
        for j in range(1, ny + 1):
            for i in range(nx + 1):
                work[k, j, i] = (phi[k, j, i + 1] - phi[k, j, i])/dx
    return work

def gradpy_3d(nx, ny, nz, dy, phi, work):
    for k in range(1, nz + 1):
        for j in range(ny + 1):
            for i in range(1, nx + 1):
                work[k, j, i] = (phi[k, j + 1, i] - phi[k, j, i])/dy
    return work

def gradpz_3d(nx, ny, nz, dzf, phi, work):
    for k in range(nz + 1):
        for j in range(1, ny + 1):
            for i in range(1, nx + 1):
                work[k, j, i] = (phi[k + 1, j, i] - phi[k, j, i])/dzf[k]
    return work


