We are trying to solve
$$
\begin{equation}
    \frac{\partial^2 p}{\partial x^2}+\frac{\partial^2 p}{\partial y^2}=f(x,y)
\end{equation}
$$
with second order central difference
$$
\begin{equation}
    \frac{p_{i-1,j} - 2p_{i,j} + p_{i+1,j}}{(\Delta x)^2} + \frac{p_{i,j-1} - 2p_{i,j} + p_{i,j+1}}{(\Delta y)^2} = f_{i,j}
\end{equation}
$$
After a proper transform in $x$-direction, we can have
$$
\begin{equation}
    \frac{(\lambda_x)_i}{(\Delta x)^2} \hat{p}_{i,j} + \frac{\hat{p}_{i,j-1} - 2\hat{p}_{i,j} + \hat{p}_{i,j+1}}{(\Delta y)^2} = \hat{f}_{i,j}
\end{equation}
$$
And this is a tridiagonal system as
$$
\begin{equation}
    \frac{1}{(\Delta y)^2} \hat{p}_{i,j-1} + \left[-\frac{2}{(\Delta y)^2} +\frac{(\lambda_x)_i}{(\Delta x)^2} \right]\hat{p}_{i,j} + \frac{1}{(\Delta y)^2} \hat{p}_{i,j+1}=\hat{f}_{i,j}
\end{equation}
$$

We are going to test the following combination of BCs and test functions

|Case | $x$-BC | $y$-BC| $p$|
|-----|--------|-------|--|
|Taylor Green| P-P | P-P | $\sin 2\pi x\sin 2\pi y$|
| Steady Channel | P-P | N-N | $\sin 2\pi x\cos\pi y$|
| Cavity | N-N | N-N| $\cos\pi x\cos\pi y$|
| Developing Channel | N-D | N-N| $\cos \frac{3\pi x}{2}\cos\pi y$|
