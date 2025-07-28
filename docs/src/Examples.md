# Examples

## Sod Shock Tube

The [Sod shock tube](https://en.wikipedia.org/wiki/Sod_shock_tube) is a canonical solution to the Euler equations involving a rightwards travelling shock wave and a leftwards travelling expansion wave. Its initial condition is given by:

```math
\left(
\begin{aligned}
\rho_L \\
P_L \\
u_L
\end{aligned}
\right) =
\left(
\begin{aligned}
1.0 \\
1.0 \\
0.0
\end{aligned}
\right),~
\left(
\begin{aligned}
\rho_R \\
P_R \\
u_R
\end{aligned}
\right) =
\left(
\begin{aligned}
0.125 \\
0.1 \\
0.0
\end{aligned}
\right)
```
where ``\rho``, ``P``, and ``u`` are the density, pressure, and velocity, respectively, the domain is ``x=[0,1]``, and the subscripts `L` and `R` refer to states to the left and right of an interface located at ``x=0.5``. Conveniently, an exact solution for this configuration can be found, allowing for validation of solvers for the Euler equations. This example will focus on setting up an `Euler1D` simulation to simulate the Sod shock tube. Comparing to the exact solution is left as a separate exercise.

## 
