module Euler1D

include("Types.jl")                 # Definition of types used to hold things like simulation data
include("Configuration.jl")         # Definition of functions responsible for configuring the simulation
include("EquationOfState.jl")       # Definition of functions that compute quantities from the equation of state
include("SpatialDerivatives.jl")    # Definition of functions to calculate spatial derivatives
include("Timestepping.jl")          # Definition of functions for timestepping.
include("ArtificialDissipation.jl") # Definition of the artificial viscosity and conductivity terms, responsible for stabilizing the solution near shocks and contact surfaces
include("EulerEquations.jl")        # Definition of functions that implement the right hand sides of the Euler equations 

"""
    cycle!( output::Simulation{T}, input::Simulation{T}, Δt::T ) where { T <: AbstractFloat }

Perform one timestep (cycle) of the solution process.

# Returns
`nothing`. Modifies the `output` argument in-place.

# Parameters
- output: A `Simulation{T}` where the state at the next timestep will be stored. Modified in-place by this function.
- input: A `Simulation{T}` representing the state at the beginning of the timestep.
- Δt: The size of the timestep. That is, the amount of time the solution is advanced in one cycle.

# Notes
- This function is generally not intended to be called directly. Instead, use the various time-stepping functions defined in `Timestepping.jl` such as `AdvanceOneCycle()`, `AdvanceNCycles()`, or `AdvanceToTime()`.
- Time integration in this function is performed using a first-order forward Euler scheme.

# Side Effects
- Modifies all of the time-dependent fields (e.g. velocity, energy, and EOS-derived fields) of the `output` argument in-place while leaving constant fields (e.g., mass, gamma) unchanged.
"""
function cycle!( output::Simulation{T}, input::Simulation{T}, Δt::T ) where { T <: AbstractFloat }
    # First compute the RHS of the momentum (velocity) and energy updates
    Momentum!( output, input )
    Energy!( output, input )
    # We now need to integrate these right hand sides to get the new values of velocity and internal energy
    # Currently this is just a basic forward Euler integration
    output.velocity .= input.velocity .+ Δt .* output.∂u∂t
    output.intenergy .= input.intenergy .+ Δt .* output.∂e∂t
    # Next we update the position of each cell edge
    for i in range( 1, output.nedges )
        output.zone_edge[i] = input.zone_edge[i] + Δt * 0.5 * ( input.velocity[i] + output.velocity[i] )
    end
    # Re-compute the location of the zone centers
    for i in 1:output.nzones
        output.zone_center[i] = 0.5 * ( output.zone_edge[i] + output.zone_edge[i+1] )
    end
    # Now update the zone sizes
    for i in range( 1, output.nzones )
        output.zone_length[i] = output.zone_edge[i+1] - output.zone_edge[i]
    end
    # Now update derived quantities using the equation of state
    EquationOfState!( output )
    # Compute the new artificial viscosity, but only if the coefficient is nonzero
    !iszero( output.Cᵥ ) && artificial_viscosity!( output )
    # as well as the new artificial conductivity, also only if the coefficient is nonzero
    !iszero( output.Cₖ ) && artificial_conductivity!( output )
    # Last, update the current simulation time
    output.time.x = output.time.x + Δt
    output.Δt.x = Δt
    output.cycles.x = output.cycles.x + 1
    if ( output.cycles.x > output.max_cycles )
        error("Current cycle ($(output.cycles.x)) exceeds maximum number of cycles ($(output.max_cycles)). Exiting.")
    end
end

end
