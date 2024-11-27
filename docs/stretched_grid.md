In wall turbulence, we must have a stretched grid in the wall-normal direction
to resolve near-wall structures. In this tutorial, we consider a channel flow
with wall-normal direction along the $z$-axis. A double stretched grid can be
obtained on an interval $[0, H]$ by a hyperbolic tangent function
$$
\begin{equation}
    z_k=\frac{H}{2}\left[1 - \frac{\tanh\beta \left(1 - \frac{2k}{N}\right)}{\tanh \beta}\right], \quad k=0,\dots, N,
\end{equation}
$$
with $N$ segments (i.e., $N + 1$ points), and $\beta$ is the stretching parameter.
A function
```python
def twoside_stretch_grid(N, H, beta):
    z = np.zeros(N + 1)
    for k in range(N + 1):
        z[k] = H/2*(1 - np.tanh(beta*(1 - 2*k/N))/np.tanh(beta))

    return z
```
is available to generate this grid. We show an example with
$$
N = 16, \quad H = 2 \text{ and } \beta = 1.2
$$
as shown below.

We can have a more stretched grid with larger $\beta$.

