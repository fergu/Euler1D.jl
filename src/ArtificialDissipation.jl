"""
    artificial_viscosity( Cᵥ::T, c::T, ρ::T, Δx::T, u₋::T, u₊::T ) where { T <: AbstractFloat }

Compute an artificial viscosity within a zone. 

# Returns
A value of type `T` representing the value of the artificial viscosity.

# Arguments
- `Cᵥ`: An O(1) coefficient to control the strength of the artificial viscosity. (Unit: ⋅)
- `c`: The speed of sound in the zone. (Unit: m/s)
- `ρ`: The density of the zone. (Unit: kg/m³)
- `Δx`: The length of the zone. (Unit: m)
- `u`: The velocity of the zone boundaries, with superscripts - and + referring to the left and right boundaries of the zone, respectively. (Unit: m/s)

# Notes
- This artificial viscosity is based on the method described by Wilkins (1980), which in turn relies upon the methods of Von Neumann and Richtmyer (1950) and Landschoff (1955). This functions by adding the artificial viscosity computed by this function to the pressure field during the [`Momentum!()`](@ref) and [`Energy!()`](@ref) updates.
- The calculation of artificial viscosity requires a gradient of velocity, ∂u/∂x. This is computed using [`∂∂x_ZoneEdgeToZoneCenter()`](@ref). See the documentation of that function for details on the numerical method used for gradient calculation.
- This function returns zero if `∂u/∂x > 0`, which will be the case for regions where the flow is expanding. This is done to restrict artificial viscosity only to regions of compression.
"""
function artificial_viscosity( Cᵥ::T, c::T, ρ::T, Δx::T, u₋::T, u₊::T ) where { T <: AbstractFloat }
    ∂u∂x = ( u₊ - u₋ ) / Δx
    if ( ∂u∂x < 0.0 ) # Only apply viscosity in regions of compression so we restrict this to shocks
        # This is the method described by Wilkins (1980) which is in turn a combination of the methods of Von Neumann and Richtmyer (1950) and Landschoff (1955)
        return -ρ * ( Cᵥ * Δx ).^2 * ∂u∂x * abs( ∂u∂x ) + 0.5 * Cᵥ * ρ * c * Δx * abs( ∂u∂x )
    else
        return 0.0
    end
end

"""
    artificial_viscosity!( state::Simulation{T} ) where { T <: AbstractFloat }

Compute the value of artificial viscosity within every zone.

# Returns
`nothing`. The input `state` is modified by this function.

# Arguments
- `state`: A `Simulation{T}` representing the current problem state that will be used to compute the artificial viscosity term.

# Notes
- This function iterates through every zone interface and calls [`artificial_viscosity()`](@ref) to update the `viscosity` state variable. See the documentation for [`artificial_viscosity()`](@ref) for details on how artificial viscosity is calculated.

# Side Effects
- Modifies the values in the `state.viscosity` vector in-place.
"""
function artificial_viscosity!( state::Simulation{T} ) where { T <: AbstractFloat }
    for i in 1:state.nzones
        state.viscosity[i] = artificial_viscosity( state.Cᵥ, state.speedofsound[i], state.density[i], state.zone_length[i], state.velocity[i], state.velocity[i+1] )
    end
end

"""
    artificial_conductivity( Cₖ::T, u::T, c₋::T, e₋::T, Δx₋::T, c₊::T, e₊::T, Δx₊::T ) where { T <: AbstractFloat }

Compute an artificial flux of internal energy across a zone interface.

# Returns
A value of type `T` representing the artifical flux of internal energy across a zone boundary.

# Arguments
- `Cₖ`: An O(1) coefficient to control the strength of the artificial conductivity. (Unit: ⋅)
- `u`: The velocity of the zone interface. (Unit: m/s)
- `c`: The speed of sound, where superscript - and + refer to the zones to the left and right of the zone interface, respectively. (Unit: m/s)
- `e`: The internal energy per unit mass. Superscript - and + refer to the zones to the left and right of the zone interface, respectively. (Unit: m²/s²)
- `Δx`: Length of the zone. Superscript - and + refer to the zones to the left and right of the zone interface, respectively. (Unit: m)

# Notes
The artificial conductivity is modeled as a Fickian diffusivity. That is, the flux of energy across a zone boundary, fₑ, is described by
    fₑ = -κ ∂e/∂x
where e is the internal energy per unit mass and κ is the (artificial) conductivity coefficient.
The gradient of internal energy is treated with a simple forward finite difference. The artificial conductivity coefficient is modeled as
    κ = Cₖ * cₘ * Δx
where
- Cₖ an O(1) coefficient.
- cₘ  is a characteristic velocity taken to be max(c̄±u, c̄), where c̄=0.5*(c₋+c₊) is the average speed of sound of the two adjacent zones and u is the velocity of the zone interface. 
- Δx is the distance between the zone centers
By convention, this leads to a positive flux if energy is diffusing in the positive x direction, and negative if it is diffusing in the negative x direction.
"""
function artificial_conductivity( Cₖ::T, u::T, c₋::T, e₋::T, Δx₋::T, c₊::T, e₊::T, Δx₊::T ) where { T <: AbstractFloat }
    ∇e = 2.0 * ( e₊ - e₋ ) / ( Δx₋ + Δx₊ )
    c̄ = 0.5 * ( c₋ + c₊ ) # The average speed of sound between the two zones
    cₘ = max( c̄ + u, c̄, c̄ + u ) # A diffusion speed, set as the maximum wave propagation speed
    κ = Cₖ * cₘ * 0.5 * ( Δx₋ + Δx₊ ) # Artificial conduction coefficient
    return -κ * ∇e
end

"""
    artificial_conductivity!( state::Simulation{T} ) where { T <: AbstractFloat }

Compute the artificial flux of energy across each zone interface in the simulation. 

# Returns
`nothing`. Modifies the `state` argument in-place. 

# Arguments
- `state`: A `Simulation{T}` object representing the problem state.

# Notes
This function calls [`artificial_conductivity()`](@ref) at each zone interface in the problem. See the documentation of [`artificial_conductivity()`](@ref) for further detail on how artificial conductivity is calculated.

# Side Effects
- Modifies the values in the `state.energy_flux` vector in-place.
"""
function artificial_conductivity!( state::Simulation{T} ) where { T <: AbstractFloat }
    for i in 1:state.nzones-1
        state.energy_flux[i+1] = artificial_conductivity( state.Cₖ, state.velocity[i+1], state.speedofsound[i], state.intenergy[i], state.zone_length[i], state.speedofsound[i+1], state.intenergy[i+1], state.zone_length[i+1] )
    end
end

export artificial_viscosity, artificial_conductivity
