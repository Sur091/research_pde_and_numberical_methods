# research_pde_and_numberical_methods
Examples for solving pde. Undergraduate research.


### Research


### Heat Equation

The Heat Equation Solution approximates the solution to the parabolic partial differential equation (PDE) given by:
```math
\frac{\partial u}{\partial t}(x,t) - \alpha^2 \frac{\partial^2 u}{\partial x^2}(x,t) = 0, \quad 0 < x < l,\quad  0 < t < T,
```
subject to the boundary conditions
```math
u(0, t)= u(l, t) = 0, \quad 0 < t < T,
```
and the initial conditions
```math
u(x, 0) = f(x), \quad 0 \leq x \leq l.
```
