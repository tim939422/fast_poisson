# DFT of complex data
For complex data of size $N$,
$$
\begin{equation}
    u_0,\dots,u_{N-1}
\end{equation}
$$
we can define the Discrete Fourier Transform (DFT) pair
$$
\begin{align}
    \operatorname{DFT}:& \quad \hat{u}_k=\frac{1}{N} \sum_{j=0}^{N-1}u_j\exp\left(-2\pi\sqrt{-1}\frac{kj}{N}\right)\\
    \operatorname{iDFT}:& \quad u_j = \sum_{k=0}^{N-1}\hat{u}_k\exp\left(2\pi\sqrt{-1}\frac{kj}{N}\right)
\end{align}
$$
This is exactly what `FFTW3` or `Numpy` computes except the normalization factor
$1/N$. In this way, we can interpret the wavenumber/frequency as (we only consider even $N$)
$$
\boldsymbol{\kappa} = (0,1,\dots,\frac{N}{2}-1,-\frac{N}{2},\dots, -1)
$$
# DFT of real data

