We need to solve a 1D Poisson equation
$$
\begin{equation}
    \frac{d^2p}{dx^2}=f(x),
\end{equation}
$$
with different boundary condition on an interval $[0,1]$. Now, we consider the following combinations
which exist in the pressure Poisson solver in incompressible NS calculations

| Boundary Condition| $p(x)$| $p^{\prime}(x)$| $f(x)$|
|-------------------|-------|----------------|-------|
| Periodic| $\sin2\pi x$ | $2\pi\cos 2\pi x$ | $-(2\pi)^2\sin2\pi x$ |
|Neumann-Neumann| $\cos\pi x$|$-\pi\sin \pi x$| $-\pi^2\cos\pi x$|
|Neumann-Dirichlet| $\cos \frac{3\pi x}{2}$ | $-\frac{3\pi}{2}\sin \frac{3\pi x}{2}$ | $-\left(\frac{3\pi}{2}\right)^2\cos \frac{3\pi x}{2}$|