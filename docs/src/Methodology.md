# Methodology

!!! note
    The documentation on this page is mainly provided for informational purposes. While it is always good to understand what a numerical code is doing before relying on the answers it gives, this information is not strictly necessary to run Euler1D.

## Overview

This package is a solver for the Euler equations in one dimension on a Lagrangian mesh with an ideal gas equation of state. Thermodynamic quantities (e.g., internal energy, density, pressure) are located at zone centers, and kinematic quantities (e.g., position, velocity) are located at zone edges. As this code is Lagrangian, the mass within a zone is constant and is stored at zone centers. Zone centers are indexed with integers (e.g., ``i``, ``i+1``, ``i+2``, etc, where ``i`` is an integer), and zone edges are fractionally indexed (e.g., ``i-1/2``, ``i+1/2``, etc)

Spatial discretization is performed using central differences, making the solution second order in smooth regions though the solution may degrade to first order in the neighborhood of discontinuities such as shocks. Further details on the governing equations, including their spatial discretization, is shown in [Equations of Motion](@ref equations-of-motion-theory).

Time discretization is performed using a first-order explicit Euler integration. For a parameter ``f`` and the value of the right hand side of the governing equation for ``f``, ``\mathrm{RHS}`` (such as would be found using the equations shown in [Equations of Motion](@ref equations-of-motion-theory)), the time integration from time index ``n`` to ``n+1`` can be written as
```math
\frac{ \partial f }{\partial t} = \frac{ f^{n+1} - f^{n} }{ t^{n+1} - t^n } = \mathrm{RHS}^n .
```
Recognizing that ``t^{n+1} - t^n = \Delta t`` is the timestep size, this equation can be rearranged to obtain the new state at the next time index as:
```math
\boxed{ f^{n+1} = f^n + \Delta t \cdot \mathrm{RHS}^n }
```

## [Equations of Motion](@id equations-of-motion-theory)

### Momentum Equation

The momentum equation solved by this package is given by equation 3 in von Neumann and Richtmyer [1]:
```math
\frac{\partial u}{\partial t} = - \frac{1}{\rho_0} \frac{\partial P}{\partial x}
```
where ``u`` is the velocity, ``P`` is the pressure, ``\rho_0`` is the initial density, and ``t`` and ``x`` are time and position, respectively.

#### Numerical Discretization

Semi-discretizing the momentum equation at the zone interface with index ``i+1/2``:
```math
\frac{\partial u_{i+1/2}}{\partial t} = - \frac{1}{\rho_0} \frac{ P_{i+1} - P_i }{ x_{i+1} - x_i } .
```
Internally, the code tracks the location of zone edges (e.g., ``x_{i+1/2}``), and infers the position of a zone center from the location of the left and right zone edges as ``x_i = ( x_{i+1/2} + x_{i-1/2} ) / 2``. Therefore, the above expression can be written in terms of the locations of the zone edge locations as
```math
\frac{\partial u_{i+1/2}}{\partial t} = - \frac{1}{\rho_0} \frac{ P_{i+1} - P_i }{ \frac{1}{2} (x_{i+3/2} + x_{i+1/2} ) - \frac{1}{2} ( x_{i+1/2} + x_{i-1/2} ) } .
```
Noting that ``x_{i+1/2} - x_{i-1/2} = \Delta x_i``, the length of the ``i``-th zone, this can be further simplified to obtain
```math
\frac{\partial u_{i+1/2}}{\partial t} = - \frac{1}{\rho_0} \frac{ P_{i+1} - P_i }{ \frac{1}{2} (\Delta x_{i+1} + \Delta x_{i}) } .
```
Finally, recognizing that ``\rho_{0,i} \Delta x_i = m_i`` is the (constant) mass in zone ``i``, the discretized equation is obtained as
```math
\boxed{ \frac{\partial u_{i+1/2}}{\partial t} = - \frac{1}{\bar{m}_{i+1/2}} ( P_{i+1} - P_i ) }
```
where ``\bar{m}_{i+1/2} = \frac{1}{2} ( m_{i+1} + m_{i} )`` is the average mass of the zones to the left and right of the interface.

### Energy Equation

The equation for internal energy per unit mass is given in equation 4 of von Neumann and Richtmyer [1] as
```math
\frac{\partial e}{\partial t} = - \frac{P}{\rho_0} \frac{\partial u}{\partial x}
```
where ``e`` is the internal energy per unit mass.

#### Numerical Discretization

A similar procedure as the momentum equation is undertaken for the discretization of the energy equation. Semi-discretizing the energy equation for a zone with index ``i`` yields
```math
\frac{\partial e}{\partial t} = - \frac{P_i}{\rho_0} \frac{u_{i+1/2} - u_{i-1/2}}{x_{i+1/2} - x_{i-1/2}} .
```
Noting that ``x_{i+1/2} - x_{i-1/2} = \Delta x_i`` is the length of the ``i``-th zone, and recalling that ``\rho_0 \Delta x_i = m_i``, this expression can be simplified to
```math
\boxed{\frac{\partial e}{\partial t} = - \frac{P_i}{m_i} (u_{i+1/2} - u_{i-1/2} ) }
```

## Equations of State

This code directly evolves momentum (or, more precisely, the velocity of zone edges and thus zone length) and internal energy per unit mass within each zone. Additionally, by virtue of using a Lagrangian method, mass within each zone is constant. All other quantities are derived from these quantities assuming an ideal gas equation of state. This section enumerates the derived quantities and the mathematical formulas used to compute them.

### Density

Density is derived from the size of a zone at each timestep according to
```math
\rho_i = \frac{ m_i }{ \Delta x_i }
```
where ``\rho_i`` is the density of the ``i``-th zone, ``m_i`` is the mass contained in the zone, and ``\Delta x_i`` is the length of the zone. The mass initially inside each zone is computed at the time of simulation initialization. Given a function for density and a known zone length, the mass within each zone is initialized as
```math
m_i = \rho_i \Delta x_i .
```
The function used to compute density in a given zone is [`EOS_Density()`](@ref).

### Pressure

Pressure is derived from the density and internal energy per unit mass according to
```math
P_i = (\gamma_i - 1) \rho_i e_i
```
where ``P_i`` is the pressure in the ``i``-th zone, ``\gamma`` is the the ratio of specific heats for the gas in the zone, ``\rho_i`` is the density, and ``e_i`` is the internal energy per unit mass. The function used to compute pressure in a zone is [`EOS_Pressure()`](@ref), and internally the density is computed using [`EOS_Density()`](@ref).

### Speed of Sound

The speed of sound, ``c``, within a zone is calculated as
```math
c_i = \sqrt{ \gamma_i \frac{P_i}{\rho_i}} .
```
Speed of sound is computed using the function [`EOS_SpeedOfSound()`](@ref). Internally, the pressure in this expression is found using [`EOS_Pressure()`](@ref) and the density in this expression is found using [`EOS_Density()`](@ref). 
!!! note
    It is important to highlight that the equations shown in this section are written using derived quantities (e.g., pressure being a function of density) in order to be similar to existing literature and thus be more understandable. However, this does not exactly align with what this package does internally. Internally, any derived quantity is evaluated directly from evolved quantities and the appropriate equations of state, *not* from any of the derived fields.

## [Artificial Dissipation](@id artificial-dissipation-theory)

Numerical solutions often form oscillations in the neighborhood of solution discontinuities, such as shock waves. To combat this, an artificial viscosity term is added to the pressure field, and an artificial conductivity term is added to the energy equation. These terms act to spread discontinuities over several zones and thus dampen these oscillations. This section describes the mathematical form of each of these terms.

### Artificial Viscosity

Von Neumann and Richtmyer [1] initially introduced an artificial viscosity term that is added to the pressure field while solving the governing equations, ``P_t = P_a + q_v``, where ``P`` is the pressure, ``P_a`` is the actual (computed from the equation of state) pressure, and ``q_v`` is an artifical addition. Mathematically, the expression used in [1] is
```math
q_{vnr} = - C_v \rho (c\Delta x)^2 \frac{\partial u}{\partial x} \cdot \left| \frac{\partial u}{\partial x} \right|
```
where ``q_{vnr}`` is the artificial von Neumann-Richtmyer viscosity, ``C_v`` is an ``O(1)`` dimensionless coefficient (and not to be confused with ``\hat{c}_v``, the dimensionless heat capacity at constant volume), and ``c`` is the local speed of sound. Later, Landshoff [2] introduced a different definition that was instead linear in the velocity gradient:
```math
q_l = C_v \rho \Delta x c \left| \frac{\partial u}{\partial x} \right|
```
where ``q_l`` is the Landshoff artificial viscosity. Wilkins [3] combined these two definitions into a single expression,
```math
\boxed{ q_v = c_{vnr} + q_l = - C_v \rho (c\Delta x)^2 \frac{\partial u}{\partial x} \cdot \left| \frac{\partial u}{\partial x} \right| + C_v \rho \Delta x c \left| \frac{\partial u}{\partial x} \right| }
```
This is the expression used for artificial viscosity in this package and is computed by the function [`artificial_viscosity()`](@ref). Numerically, this artificial viscosity is zone-centered to be aligned with the pressure field. Derivatives of velocity are computed identically to the method used in the computation of pressure, ``\left( \partial u/\partial x \right)_i = ( u_{i+1/2} - u_{i-1/2} ) / \Delta x_i``.

!!! info
    The same tuning coefficient, ``C_v``, is used for both the von Neumann-Richtmyer and Landshoff components of artificial viscosity.

### Artificial Conductivity

The artificial viscosity is based on velocity gradients, and so will not be active in regions where the solution is continuous in velocity but discontinuous in internal energy or density, such as would be the case in the neighborhood of a contact surface. To address this, this package additionally adds artificial diffusion of internal energy to smooth discontinuities in the internal energy field. Mathematically, this is treated as a [Fickian diffusion](https://en.wikipedia.org/wiki/Fick%27s_laws_of_diffusion), where the flux of internal energy across a zone edge is
```math
J_e = \kappa \frac{\partial e}{\partial x}
```
where ``J_e`` is the flux of internal energy across a zone edge and ``\kappa`` is a diffusion coefficient, modeled as
```math
\kappa = C_k c_m \Delta x
```
where ``C_k`` is an ``O(1)`` dimensionless coefficient and ``c_m`` is a characteristic velocity. By convention, a positive value of ``J_e`` represents energy flux in the positive x direction, and a negative value of ``J_e`` represents flux in the negative x direction. As ``J_e`` represents a flux and is thus located on zone edges, this package computes ``c_m`` at a zone edge as
```math
c_{m,i+1/2} = \mathrm{max}( \bar{c}_{i+1/2} \pm u_{i+1/2}, \bar{c}_{i+1/2} )
```
where ``\bar{c}_{i+1/2} = \frac{1}{2} ( c_i + c_{i+1} )`` is the average speed of sound of the zones to the left and right of the interface. 

The flux of internal energy is added to the energy equation as
```math
\frac{\partial e}{\partial t} = \cdots + q_e
```
where ``q_e`` is
```math
q_e = \frac{\partial}{\partial x} \left( \kappa \frac{\partial e}{\partial x} \right) = \frac{\partial}{\partial x}\left( J_e \right) .
```
Numerically, this term is discretized using the same procedure as the rest of the energy equation to obtain
```math
q_{e,i} = \frac{ J_{e,i+1/2} - J_{e,i-1/2} }{x_{i+1/2} - x_{i-1/2}} = \frac{1}{\Delta x_i} ( J_{e,i+1/2} - J_{e,i-1/2} )
```
Internally, this term is computed using [`artificial_conductivity()`](@ref).

!!! caution
    Artificial conductivity is useful for suppressing oscillations near contact surfaces and stabilizing problems with multiple shock waves, but can also result in material interfaces and contact surfaces becoming unphysically diffuse over long time periods. Generally this term should be small compared to the artificial viscosity, ``C_k \ll C_v``.

## References

* [1] von Neumann, J. and Richtmyer, R. D. "A method for the numerical calculation of hydrodynamic shocks". J. Appl. Phys. (21) pp 232-237 (1950) DOI: [10.1063/1.1699639](https://doi.org/10.1063/1.1699639)
* [2] Landshoff, R., "A Numerical Method for Treating Fluid Flow in the Presence of Shocks", Los Alamos Scientific Laboratory Report LA-1930 (1955) [Link](https://apps.dtic.mil/sti/tr/pdf/ADA382679.pdf)
* [3] Wilkins, M. L., "Use of Artificial Viscosity in Multidimensional Fluid Dynamic Calculations". J. Comp. Phys. (36) pp 281-303 (1980) DOI: [10.1016/0021-9991(80)90161-8](https://doi.org/10.1016/0021-9991(80)90161-8)
