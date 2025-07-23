# Euler1D

This is a package that solves the one-dimensional compressible Euler equations on a Lagrangian (moving) mesh and assuming an ideal gas equation of state. 

The primary contribution of this package is intended to be as a demonstration of how a hydrodynamics code can be implemented entirely in Julia, including problem setup, the simulation code itself, and results postprocessing. In this way, there is no need to maintain a separate parser to interpret "input decks" for configuring the intial problem state. The language used to set up the problem is the same as the one used to perform the simulation. As long as the language supports it, it can be used with no modification to this package. Ideally, this can all be done while maintaining relative parity in execution time with similar codes written in compiled languages (C, C++, Fortran, etc).

This README is broken into several sections:

* `Installation`: How to install this package to be used with Julia
* `Usage`: Information on how to use this package, including simulation setup, execution, and processing of results
* `Methodology`: Description of the solution method used in this package. Includes the equations of motion and artificial terms.
* `Roadmap`: Ideas for additional functionality that might be useful or nice to have. These are not necessarily in progress, and so may also serve as good suggestions for anyone interested in contributing.
* `References`: Citation information for references used in developing this package

# Installation

Installation of this package is straightforward. At the Julia REPL, type

```julia
julia> # Press ] key
pkg> add https://www.github.com/fergu/Euler1D.git
```

or, alternatively

```julia
julia> using Pkg
julia> Pkg.add("https://www.github.com/fergu/Euler1D.git")
```

## A note on data analysis and output

There are numerous packages in the Julia ecosystem to handle plotting (e.g. `Plots.jl`, `Makie.jl`), writing to disk (e.g. `CSV.jl`, `DataFrames.jl`, `HDF5.jl`), and any other data analysis routine that might be desired. The people developing these packages (and others like them) are immensely talented, and one would be well-advised to utilize their hard work if a package with suitable functionality exists.

For this reason, this package actively avoids making assumptions about what a user will want to do with the simulation data after the simulation is complete. This means that there is no built-in functionality to plot simulation data, write it to disk, or perform any other sort of postprocessing. Instead, all simulation outputs are simply elementary Julia types (e.g. `Vector{Float64}`), and can be manipulated using the usual built-in Julia methods, or by using any suitable package that provides the required functionality. 

# Usage

At its core, using `Euler1D.jl` relies on a single Julia file to describe simulation setup, running of the simulation, and processing of outputs. Each of these aspects is detailed in the following subsections. Some complete examples of a simulation setup can be found in the `examples` subdirectory.

## Initialization

A simulation is described by a dictionary of type `Dict{String, Any}` that contains key, value pairs describing various parts of the simulation. The minimal set of parameters can be obtained via the function `DefaultSimulationParameters()`. The dictionary returned by this function can then be modified to match the desired simulation setup. See `? DefaultSimulationParameters()` for more information on the various settings, their meanings, and their default values.

Of particular note, four keys in the dictionary returned by `DefaultSimulationParameters()` are initially set to `nothing` and must be user supplied. These are `init_density_function`, `init_velocity_function`, `init_pressure_function` and `init_gamma_function`, and refer to `Function`s that describe the initial density, velocity, pressure, and ratio of specific heats (gamma) as a function of position, `x`. These should be set to functions with the signature `MyExampleFunction( x::Float64 )`.  For example, consider the following functions for `density`, `velocity`, `pressure`, and `gamma` to set up Sod's shock tube initial condition:

```julia
function density( x::Float64 )
    if ( x < 0.5 )
        return 1.0
    else
        return 0.125
    end
end

function velocity( x::Float64 )
    return 0.0
end

function pressure( x::Float64 )
    if ( x < 0.5 )
        return 1.0
    else
        return 0.1
    end
end

function gamma( x::Float64 )
    return 1.4
end
```

These functions would then be assigned to the simulation by modified the appropriate dictionary key. As an example, if using `DefaultSimulationParameters()` to provide a default dictionary, these functions can be assigned using the following commands:
```julia
julia> sim_dict = DefaultSimulationParameters()
julia> sim_dict["init_density_function"] = density
julia> sim_dict["init_velocity_function"] = velocity
julia> sim_dict["init_pressure_function"] = pressure
julia> sim_dict["init_gamma_function"] = gamma
```

Once all settings have been configured to match the problem configuration, the simulation can be initialized using the function `InitializeSimulation( parameters::Dict{String, Any} )`. Continuing from the example configuration in the previous paragraph, one might call:
```julia
julia> init_state = InitializeSimulation( sim_dict )
```

At this point, `init_state` is a `Simulation{T}` type (where `T` is some concrete subtype of AbstractFloat, most likely `Float64` unless you pick something else for some reason) representing the initial state of the simulation. This struct contains information on the simulation state, and is what is used to actually perform the simulation. For more detail, see the next section.

## Running the simulation

The simulation can then be advanced to the final time using the `AdvanceToTime` function, where `exact=true` tells the code to adjust the final timestep to end as close as possible to the final time:
```julia
end_state = AdvanceToTime( init_state, init_state.tₑ; exact=true )
```

If you want to have more complex output to, for example, get the problem state every 0.01 seconds, you might consider a slightly more complex approach:
```julia
function simulate( init_state::Simulation{T} ) where { T <: AbstractFloat }
    plot_Δt = 0.01 # Frequency to output plots
    next_plot = init_state.t₀ + plot_Δt
    new_state = deepcopy(init_state) # This is our output. We want to use a copy of the input.

    while ( new_state.time.x < end_time )
        new_state = AdvanceToTime( new_state, next_plot; exact=true ) # Advance to t=next_plot
        # Do something with the data, such as plot, save to file, etc
        # ...
        # Then set our next time to stop and plot
        next_plot = next_plot + plot_Δt
    end
    return new_state # Return new_state in case we want to do any final manipulations
end

end_state = simulate( init_state )
```

There are also more complex options to advance by a single timestep or a fixed number of timesteps. These might be useful if you're having issues with problem stability, for example.

## Timestep controls

By default, the simulation utilizes adaptive timestepping. The timestep is computed at every cycle based on the time it would take an acoustic wave moving at the local speed of sound, `c`, to propagate across a zone of length `Δx`. The minimum of this value over the entire domain is then multiplied by the supplied `CFL` to determine the timestep for the next cycle. This approach is chosen as the default to aid in stability as the timestep will be reduced to accomodate regions with densely-clustered zones.

If you would prefer to have a fixed timestep, each of the time-advancing functions (`AdvanceToTime`, `AdvanceOneCycle`, `AdvanceNCycles`) support an optional `Δt` argument. If supplied, this value will be utilized in place of the adaptive timestep, and the function used to determine the timestep will be bypassed. Be aware, however, that this will likely require decreasing the CFL number or increasing artificial viscosity/conductivity to maintain stability. For more information, see the documentation of each function.

## Artificial Viscosity and Conductivity

This code implements an artificial viscosity term that is added to the pressure field, and an artificial conductivity term that is added to the energy evolution equation. These terms act to smooth out oscillations in the solution in the neighborhood of shock waves and contact surfaces, respectively. More detail on these expressions is described in the `Methodology` section, below.

Each of these terms introduces a tuneable O(1) coefficient that impacts how strongly each term is applied. Default values were chosen to hopefully be applicable to a range of problems, but if you find these values to be insufficient they can be overridden using the fields:

```
artificial_viscosity_coefficient
```
for the artificial viscosity and
```
artificial_conductivity_coefficient
```
for the artificial conductivity.

# Methodology

Much of the methodology in this package is based upon the work by von Neumann and Richtmyer (1950) [1], and augmentations made by Landshoff (1955) [2] and Wilkins (1980) [3]. The first subsection describes the equations of motion that are solved and the general solution process. The second subsection outlines the form of the artificial viscosity and conductivity terms that are used to smooth oscillations in the solution.

## Equations of Motion and Solution Process

This package solves the 1D Euler equations in Lagrangian form. That is, the equations are solved on a mesh that moves with the flow. The form of the equations is given by Equations 3 - 5 in [1]. Within this packge, the spatial evolution of the equations is treated separately from the temporal evolution. That is, the equations are re-organized into a form where the time derivative is on the left hand side of the equation, and the spatial derivatives are on the right.

The equations are solved on a staggered mesh, with thermodynamic variables (energy, pressure, density) located at zone centers, and dynamic variables (velocity) located at zone edges.

Practically, the equations are not evolved continuously, but rather over fixed time intervals. Each time interval is referred to as a "cycle". The amount of time that elapses in a given cycle is controlled by the time step, and the methodology for how the time step is determined is described in the Timestep Controls section, above.

Within each cycle, a series of steps take place in order to advance the solution to the next time:
1. The right hand side of the momentum equation (Eq. 3 in [1]) is evaluated using the values of parameters at the current time
1. The right hand size of the energy equation (Eq. 4 in [1]) is evaluated using the values of parameters at the current time
1. The momentum and energy equations are advanced to the new time (`t + Δt`) using a first-order forward Euler integration
1. The edges of each zone are moved according to an average of the velocity at that zone edge from the current time (`t`) and new time (`t + Δt`), times the time step.
1. The new density of each zone is determined by taking the mass of each zone (assumed to be fixed due to the Lagrangian reference frame) divided by the zone size at the new time.
1. The pressure in each zone is similarly updated using the internal energy and density at the new time.
1. New values of artificial viscosity and conductivity are calculated from the complete state at the new time.
1. The cycle is complete

## Artificial Viscosity and Conductivity

Of particular note, this package introduces two "artificial" terms to the equations of motion:

1. An artificial viscosity to damp oscillations in the neighborhood of shock waves (Introduced by [1] and augmented by [2] and [3]).
1. An artificial conductivity to damp oscillations in the neighborhood of contact discontinuities (e.g., jumps in density/internal energy without a corresponding jump in pressure or velocity). This is implemented as a diffusion of internal energy.

# Roadmap

This section describes ideas for future functionality that could be added. Underneath each bullet point below are a few notes on additional functionality that might need to be implemented to support that idea.

Just because something is included here does not mean it is being actively worked on, and just because something isn't included doesn't mean it isn't of interest. If you are interested in contributing to any of these ideas or one of your own, feel free to open an issue or a pull request.

- Add examples
- Add tests
- Support for web documentation using `Documenter.jl`
- Register this package with the Julia package index to simplify installation (After V1.0)
- Support for RANS models to describe, e.g., mixing across a material interface
  - Will require treatment of mixtures, including updating calls to the equation of state to calculate mixture properties
  - Will probably require tracking of individual species mass fraction fields, which is not currently supported.
- Support for an adaptable initial distribution of grid points, with the hope that there is less difference between the timestep requirements of the smallest and largest zones.
  - This could, for example, use the supplied initial condition functions to compute a local speed of sound and adjust the location of the next zone based on that number.
  - The downside is that information on, for example, the zone that corresponds to a given initial location is a bit harder to track
  - This is potentially addressed by introducing passive tracers that move with the local velocity
- Support for passive tracer particles
  - Particles that move exactly with the local flow velocity. This is potentially useful for things like interface tracking given the Lagrangian nature of the solution.
  
# References

* [1] von Neumann, J. and Richtmyer, R. D. "A method for the numerical calculation of hydrodynamic shocks". J. Appl. Phys. (21) pp 232-237 (1950)
* [2] Landshoff, R., "A Numerical Method for Treating Fluid Flow in the Presence of Shocks", Los Alamos Scientific Laboratory Report LA-1930 (1955)
* [3] Wilkins, M. L., "Use of Artificial Viscosity in Multidimensional Fluid Dynamic Calculations". J. Comp. Phys. (36) pp 281-303 (1980)
