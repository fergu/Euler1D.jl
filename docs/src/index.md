```@meta
CurrentModule=Euler1D
```
# Euler1D

This is a package that solves the one-dimensional compressible Euler equations on a Lagrangian (moving) mesh and assuming an ideal gas equation of state. 

The primary contribution of this package is intended to be as a demonstration of how a hydrodynamics code can be implemented entirely in Julia, including problem setup, the simulation code itself, and results postprocessing. In this way, there is no need to maintain a separate parser to interpret "input decks" for configuring the intial problem state. The language used to set up the problem is the same as the one used to perform the simulation. As long as the language supports it, it can be used with no modification to this package. Ideally, this can all be done while maintaining relative parity in execution time with similar codes written in compiled languages (C, C++, Fortran, etc).

This README is broken into several sections:

## Installation

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

## Usage

At its core, using `Euler1D.jl` relies on a single Julia file to describe simulation setup, running of the simulation, and processing of outputs. The [Examples](@ref) page shows how to set up and run a simulation, and post-process its outputs. The `examples` subdirectory also contains ready-to-run examples as well.

## Roadmap

This section describes ideas for future functionality that could be added. Underneath each bullet point below are a few notes on additional functionality that might need to be implemented to support that idea.

Just because something is included here does not mean it is being actively worked on, and just because something isn't included doesn't mean it isn't of interest. If you are interested in contributing to any of these ideas or one of your own, feel free to open an issue or a pull request.

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
  
