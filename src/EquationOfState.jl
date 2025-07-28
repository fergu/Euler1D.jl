"""
    EOS_Density( mass::T, Δx::T ) where { T <: AbstractFloat } = mass / Δx

Compute the density of the fluid in a zone.

# Returns
A scalar of type `T` representing the density of a zone.

# Arguments
- `mass`: The total mass contained within the zone. (Unit: kg)
- `Δx`: The size of the zone. (Unit: m)

# Notes
Density is calculated as:
    ρ = mass / Δx
As this is the physical definition of density, this calculation does not assume any particular equation of state.
"""
EOS_Density( mass::T, Δx::T ) where { T <: AbstractFloat } = mass / Δx

"""
    EOS_Density!( simulation::Simulation{T} ) where { T <: AbstractFloat }

Update all zone densities.

# Returns
`nothing`. Modifies the `simulation` input in-place. See the Side Effects section for further detail.

# Arguments
- `simulation`: A `Simulation{T}` representing the simulation state.

# Notes
Calls [`EOS_Density()`](@ref) to compute the density in each zone. See the documentation of that function for further details.

# Side Effects
- Modifies the values in the `simulation.density` vector in-place
"""
function EOS_Density!( simulation::Simulation{T} ) where { T <: AbstractFloat }
    @inbounds for i in 1:simulation.nzones
        simulation.density[i] = EOS_Density( simulation.mass[i], simulation.zone_length[i] )
    end
end

"""
    EOS_Pressure( γ::T, ρ::T, e::T ) where { T <: AbstractFloat }
    EOS_Pressure( γ::T, mass::T, Δx::T, e::T ) where { T <: AbstractFloat }

Compute the pressure of a zone using an ideal gas equation of state.

# Returns
A scalar of type `T` representing the pressure within the zone.

# Arguments
- `γ`: The ratio of specific heats of the fluid in the zone. (Unit: ⋅)
- `ρ`: The density of the zone. (Unit: kg/m³)
- `mass`: The total mass contained within the zone. (Unit: kg)
- `Δx`: The size of the zone. (Unit: m)
- `e`: The internal energy per unit mass of the zone. (Unit: m²/s²)

# Notes
The pressure is calculated as:
    P = ( γ - 1 ) * ρ * e
The four-parameter version of this function computes density using [`EOS_Density()`](@ref). See the documentation for that function for further details.
"""
EOS_Pressure( γ::T, ρ::T, e::T ) where { T <: AbstractFloat } = ( γ - 1.0 ) * ρ * e
EOS_Pressure( γ::T, mass::T, Δx::T, e::T ) where { T <: AbstractFloat } = EOS_Pressure( γ, EOS_Density( mass, Δx ), e )

"""
    EOS_Pressure!( simulation::Simulation{T} ) where { T <: AbstractFloat }

Update all zone pressures using an ideal gas equation of state.

# Returns
`nothing`. Modifies the `simulation` input in-place. See the Side Effects section for further detail.

# Arguments
- `simulation`: A `Simulation{T}` representing the simulation state.

# Notes
Calls [`EOS_Pressure()`](@ref) to compute the pressure in each zone. See the documentation of that function for further details.

# Side Effects
Modifies the values in the `simulation.pressure` vector in-place.
"""
function EOS_Pressure!( simulation::Simulation{T} ) where { T <: AbstractFloat }
    @inbounds for i in 1:simulation.nzones
        simulation.pressure[i] = EOS_Pressure( simulation.gamma[i], simulation.mass[i], simulation.zone_length[i], simulation.intenergy[i] )
    end
end

"""
    EOS_SpeedOfSound( γ::T, P::T, ρ::T ) where { T <: AbstractFloat }
    EOS_SpeedOfSound( γ::T, e::T, mass::T, Δx::T ) where { T <: AbstractFloat }

Compute the speed of sound in a zone using an ideal gas equation of state.

# Returns
A scalar of type `T` representing the speed of sound in the zone.

# Arguments
- `γ`: The ratio of specific heats in the zone. (Unit: ⋅)
- `e`: The internal energy per unit mass in the zone. (Unit: m²/s²)
- `P`: The pressure in the zone. (Unit: kg/(m⋅s²))
- `ρ`: The density of the fluid in the zone. (Unit: kg/m³)
- `mass`: The mass contained in the zone. (Unit: kg)
- `Δx`: The length of the zone. (Unit: m)

# Notes
The speed of sound is calculated as:
    c = √( γ * P / ρ )
The four-parameter version of this function computes density using [`EOS_Density()`](@ref) and pressure using [`EOS_Pressure()`](@ref). See the documentation of those functions for further details.
"""
EOS_SpeedOfSound( γ::T, P::T, ρ::T ) where { T <: AbstractFloat } = sqrt( γ * P / ρ )
EOS_SpeedOfSound( γ::T, e::T, mass::T, Δx::T ) where { T <: AbstractFloat } = EOS_SpeedOfSound( γ, EOS_Pressure( γ, mass, Δx, e ), EOS_Density( mass, Δx ) )

"""
    EOS_SpeedOfSound!( simulation::Simulation{T} ) where { T <: AbstractFloat }

Update all speed of sound values using an ideal gas equation of state.

# Returns
`nothing`. Modifies the `simulation` input in-place. See the Side Effects section for further detail.

# Arguments
- `simulation`: A `Simulation{T}` representing the simulation state.

# Notes
Calls [`EOS_SpeedOfSound()`](@ref) to calculate the speed of sound in each zone. See the documentation of that function for further details.

# Side Effects
Modifies the values in the `simulation.speedofsound` vector in-place.
"""
function EOS_SpeedOfSound!( simulation::Simulation{T} ) where { T <: AbstractFloat }
    @inbounds for i in 1:simulation.nzones
        simulation.speedofsound[i] = EOS_SpeedOfSound( simulation.gamma[i], simulation.intenergy[i], simulation.mass[i], simulation.zone_length[i] )
    end
end

"""
    EquationOfState!( simulation::Simulation{T} ) where { T <: AbstractFloat }

Update all densities, pressures, and speeds of sound using the equation of state.

# Returns
`nothing`. Modifies the `simulation` input in-place. See the Side Effects section for further detail.

# Arguments
- `simulation`: A `Simulation{T}` representing the simulation state.

# Notes
Calls [`EOS_Density!()`](@ref), [`EOS_Pressure!()`](@ref), and [`EOS_SpeedOfSound!()`](@ref) to compute the density, pressure, and speed of sound in each zone respectively. See the documentation of those functions for further details.

# Side Effects
- Modifies the values in the `simulation.density` vector in-place.
- Modifies the values in the `simulation.pressure` vector in-place.
- Modifies the values in the `simulation.speedofsound` vector in-place.
"""
function EquationOfState!( simulation::Simulation{T} ) where { T <: AbstractFloat }
    # Update the solution state using the equation of state together with the state variables
    # Compute the density of each zone
    EOS_Density!( simulation )
    # Compute pressure of each zone
    EOS_Pressure!( simulation )
    # Compute the speed of sound of each zone
    EOS_SpeedOfSound!( simulation )
end

# Only the pure EOS calculation functions are exported to aid in problem setup.
# The versions that modify the simulation state should only be needed internally, so they are not exported.
export EOS_Density, EOS_Pressure, EOS_SpeedOfSound
