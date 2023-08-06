# finite-vis
Visualizer for finite-difference method solver written in go using OpenGL and GLFW. Mainly written as an exploration into finite difference methods, and an excuse to learn go.

![Alt text](docs/heat_bar.gif?raw=true)
Example: Copper bar initially at 0K, both ends are held to 300K and 400K respectively.

<img src="https://latex.codecogs.com/gif.latex?\frac{\partial u}{\partial t} \approx \frac{u^{k+1}_{n} - u^{k}_{n}}{\Delta t}" /> 

# Heat Equation
Differential equation describing how heat diffuses through a medium. Equations of a similar form can be used for general diffusion problems.

$\frac{\partial u}{\partial t} = \alpha\nabla^2u $

where

$u(t, x)$ is the temperature at point $x$ at time $t$

In a cartesian plane, the differential equation can be rewritten in terms of the $x, y, z$ components:

$\frac{\partial u}{\partial t} = \alpha(\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} + \frac{\partial^2 u}{\partial z^2}) $

Simplifying to a 1-D case:

$\frac{\partial u}{\partial t} = \alpha\frac{\partial^2 u}{\partial x^2} $

## Finite Difference

The derivatives can be approximated by discretizing the time and space of the problem into finite chunks, and calculating approximations at each chunk.

For a cell $n$ at time $k$ we can formulate the equation using a few differenct schemes:

### Forward Difference in Time - Central Difference in Space
$\approx$

$u^{k+1}_{n}$

$u_{n}^{k+1}$

$\frac{\partial u}{\partial t} \approx \frac{u^{k+1}_{n} - u^{k}_{n}}{\Delta t}$

$\frac{\partial^2 u}{\partial x^2} \approx \frac{u^{k}_{n-1} - 2u^{k}_{n} + u^{k}_{n+1}}{\Delta x^2}$

Which can be combined and rearranged to get:

$r = \alpha \frac{\Delta t}{\Delta x^2}$

$u^{k+1}_{n} = ru^{k}_{n-1} + (1 - 2r)u^{k}_{n} + ru^{k}_{n+1}$

And is numerically stable if: $r \le 0.5$

With boundary conditions:

$ u^{0}_{0} = T_1; u^{0}_{N} = T_2; u^{0}_{n} = U_0(n) $

Where $U_0(n)$ is some initial distribution.
$u_0$ and $u_N$ are held at fixed temperatures throughout.

### Backward Difference in Time - Central Difference in Space
$\frac{\partial u}{\partial t} \approx \frac{u^{k+1}_{n} - u^{k}_{n}}{\Delta t}$

$\frac{\partial^2 u}{\partial x^2} \approx \frac{u^{k+1}_{n-1} - 2u^{k+1}_{n} + u^{k+1}_{n+1}}{\Delta x^2}$

Which can be combined and rearranged to get:

$r = \alpha \frac{\Delta t}{\Delta x^2}$

$u^{k}_{n} = -ru^{k+1}_{n-1} + (1 + 2r)u^{k+1}_{n} - ru^{k+1}_{n+1}$

Which is always numerically stable, but $r$ needs to be low for physical accuracy, so the same condition is used.
This cannot be solved directly for each directly, instead a system of equations needs to be solved for each timestep.
Taking the same boundary conditions as before, we get two additional equations:

$u^k_0 = u^{k+1}_0$ 

$u^k_N = u^{k+1}_N$ 

Putting this into matrix form yields:

$\begin{bmatrix}
1 & 0 & 0 & 0 & 0 &... & 0 & 0 & 0 & 0\\
-r & (1+2r) & -r & 0 & 0 &... & 0 & 0 & 0 & 0\\
0 & -r & (1+2r) & -r & 0 &... & 0 & 0 & 0 & 0\\
0 & 0 & -r & (1+2r) & -r &... & 0 & 0 & 0 & 0\\
\vdots & \vdots & \vdots & \vdots & \vdots & \ddots & \vdots & \vdots & \vdots & \vdots\\
 0 & 0 & 0 & 0 & 0 &...& -r & (1+2r) & -r & 0\\
 0 & 0 & 0 & 0 & 0 &...& 0 & -r & (1+2r) & -r\\
 0 & 0 & 0 & 0 & 0 &...& 0 & 0 & 0 & 1\\
\end{bmatrix}
\begin{bmatrix}u_0\\u_1\\u_2\\u_3\\\vdots\\u_{N-2}\\u_{N-1}\\u_N\\\end{bmatrix}^{k+1} = 
\begin{bmatrix}u_0\\u_1\\u_2\\u_3\\\vdots\\u_{N-2}\\u_{N-1}\\u_N\\\end{bmatrix}^{k}$

Which you can write as 

$\bold{C}\bold{u}^{k+1}=\bold{u}^{k}$

To solve in a naive way would to calculate the inverse of $\bold{C}$ and perform a matrix multiplication:

$\bold{u}^{k+1}=\bold{C}^{-1}\bold{u}^{k}$

A better method would be Gausian elimination, but further inspection of the full matrix reveals that it is a tridiagonal matrix, which can be solved by the tridiagonal matix algorithm (also known as the Thomas algorithm) in $O(n)$ operations instead of $O(n^3)$. This requires the matrix to be of the form:

$\begin{bmatrix}
b_1 & c_1 & & & 0\\
a_2 & b_2 & c_2 & &\\
& a_3 & b_3 & \ddots &\\
& & \ddots & \ddots & c_{n-1}\\
0 & & & a_n & b_n\\
\end{bmatrix}
\begin{bmatrix}x_1\\x_2\\x_3\\\vdots\\x_n\end{bmatrix} = 
\begin{bmatrix}d_1\\d_2\\d_3\\\vdots\\d_n\end{bmatrix}$

The method followed in this code calculates new coefficients, denoted with primes.

$c'_i = 
    \begin{cases}
        \frac{c_i}{b_i}, & i = 1\\
        \frac{c_i}{b_i - a_{i} c'_{i-1}}, & i = 2,3,\ldots,n\\
    \end{cases}$

$d'_i = 
    \begin{cases}
        \frac{d_i}{b_i}, & i = 1\\
        \frac{d_i - a_id'_{i-1}}{b_i - a_{i} c'_{i-1}}, & i = 2,3,\ldots,n\\
    \end{cases}$

Then substituting back, we get the solution:

$x_i = 
    \begin{cases}
        d'_i, & i = n\\
        d'_i - c'_ix_{i+1}, & i = n-1,n-2,\ldots,1\\
    \end{cases}$

This method is used exactly, where the $\bold{x}$ vector is $\bold{u}^{k+1}$ and the $\bold{d}$ vector is $\bold{u}^{k}$


### Central Difference in Time - Central Difference in Space (Crank-Nicolson Method)
$\frac{u^{k+1}_{n} - u^{k}_{n}}{\Delta t} = \frac{\alpha}{2}[\frac{u^{k+1}_{n-1} - 2u^{k+1}_{n} + u^{k+1}_{n+1}}{\Delta x^2} + \frac{u^{k}_{n-1} - 2u^{k}_{n} + u^{k}_{n+1}}{\Delta x^2}]$

Which can be combined and rearranged to get:

$r = \alpha \frac{\Delta t}{\Delta x^2}$

$-u^{k+1}_{n-1} + 2\frac{(1 + r)}{r}u^{k+1}_{n} - u^{k+1}_{n+1} = u^{k}_{n-1} + 2\frac{(1 - r)}{r}u^{k}_{n} + u^{k}_{n+1}$

Taking the same boundary conditions as before, we get two additional equations:

$u^k_0 = u^{k+1}_0$ 

$u^k_N = u^{k+1}_N$ 

In matrix form:

$\begin{bmatrix}
1 & 0 & 0 & 0 & 0 &... & 0 & 0 & 0 & 0\\
-1 & \frac{2(1+r)}{r} & -1 & 0 & 0 &... & 0 & 0 & 0 & 0\\
0 & -1 & \frac{2(1+r)}{r} & -1 & 0 &... & 0 & 0 & 0 & 0\\
0 & 0 & -1 & \frac{2(1+r)}{r} & -1 &... & 0 & 0 & 0 & 0\\
\vdots & \vdots & \vdots & \vdots & \vdots & \ddots & \vdots & \vdots & \vdots & \vdots\\
 0 & 0 & 0 & 0 & 0 &...& -1 & \frac{2(1+r)}{r} & -1 & 0\\
 0 & 0 & 0 & 0 & 0 &...& 0 & -1 & \frac{2(1+r)}{r} & -1\\
 0 & 0 & 0 & 0 & 0 &...& 0 & 0 & 0 & 1\\
\end{bmatrix}
\begin{bmatrix}u_0\\u_1\\u_2\\u_3\\\vdots\\u_{N-2}\\u_{N-1}\\u_N\\\end{bmatrix}^{k+1} = 
\begin{bmatrix}
    u_0\\
    u_0 + \frac{2}{r}(1-r)u_1 + u_2\\
    u_1 + \frac{2}{r}(1-r)u_2 + u_3\\
    u_2 + \frac{2}{r}(1-r)u_3 + u_4\\
    \vdots\\
    u_{N-3} + \frac{2}{r}(1-r)u_{N-2} + u_{N-1}\\
    u_{N-2} + \frac{2}{r}(1-r)u_{N-1} + u_N\\
    u_N\\
\end{bmatrix}^{k}$

Which is also tridiagonal, and can be solved the same way as before.