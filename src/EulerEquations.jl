"""
    Momentum( m₋::T, P₋::T, m₊::T, P₊::T ) where { T <: AbstractFloat }

Compute the right hand side of the momentum equation at a zone interface. See the Notes section for information on the equation being solved.

# Returns
A scalar of type `T` representing the rate of change of velocity of a zone interface over time.

# Arguments
- `m`: The mass of the zone. (Unit: kg)
- `P`: The pressure of the zone (Unit: kg/(m⋅s²))

For each of these parameters, a - subscript refers to the zone to the left of the zone interface and a + subscript refers to the zone to the right of the zone interface.

# Notes
The governing equation solved in this function is:
    ∂u/∂t = (1/ρ₀) * ∂P/∂x
Through the numerical disretization, this reduces to
    ∂u/∂t = 1/m̄ ( P₊ - P₋ )
If using artificial viscosity per the method of Von Neumann and Richtmyer (1950), the artificial viscosity term should be added to the pressure field.
"""
function Momentum( m₋::T, P₋::T, m₊::T, P₊::T ) where { T <: AbstractFloat }
    m̄ = 0.5 * ( m₋ + m₊ )
    return - ( P₊ - P₋ ) / m̄ 
end

"""
    Momentum!( output::Simulation{T}, input::Simulation{T} ) where { T <: AbstractFloat }

Compute the right hand side of the momentum equation for every zone interface. 

# Returns
`nothing`. The `output` argument is modified in-place.

# Arguments
- `output`: A `Simulation{T}` object representing the output problem state. This argument is modified by the function.
- `input`: A `Simulation{T}` object representing the input problem state.

# Notes
This function calls [`Momentum()`](@ref) to perform the actual calculation. See the documentation of that function for more information.
The artificial viscosity term is added to the pressure field during the call to [`Momentum()`](@ref), per the artificial viscosity method of Von Neumann and Richtmyer (1950).

# Side Effects
- Modifies the values contained in `output.∂u∂t` in-place.
"""
function Momentum!( output::Simulation{T}, input::Simulation{T} ) where { T <: AbstractFloat }
    for i in range( 2, input.nedges - 1 )
        output.momentum_rhs[i] = Momentum( input.mass[i-1], input.pressure[i-1] + input.viscosity[i-1], input.mass[i], input.pressure[i] + input.viscosity[i] )
    end
end

"""
    Energy( ρ::T, P::T, Δx::T, u₋::T, u₊::T, q₋::T, q₊::T ) where { T <: AbstractFloat }

Compute the right hand side of the energy equation in a zone. See the Notes section for the specific equations that are solved.

# Returns
A scalar of type `T` representing the rate of change of internal energy inside the zone.

# Arguments
- `m`: The mass of the zone. (Unit: kg)
- `P`: The pressure of the zone. (Unit: kg/(m⋅s²))
- `Δx`: The length of the zone. (Unit: m)
- `u`: The velocity of the zone edges on the (-): left and (+): right of the zone. (Unit: m/s)
- `q`: The (artificial) flux of internal energy across the (-): left and (+): right zone edges. (Unit: m³/s³) 

# Notes:
The governing equation solved in this function is
    ∂e/∂t = - ( P / ρ₀ ) * ∂u/∂x + ∑q
Through discretization, this becomes
    ∂eᵢ/∂t = - ( Pᵢ / mᵢ ) * ( u₊ - u₋ ) + ∑q 
"""
function Energy( m::T, P::T, Δx::T, u₋::T, u₊::T, q₋::T, q₊::T ) where { T <: AbstractFloat }
    return -( P / m ) * ( u₊ - u₋ ) + ( q₋ - q₊ ) / Δx
end

"""
    Energy!( output::Simulation{T}, input::Simulation{T} ) where { T <: AbstractFloat }

Compute the right hand side of the energy equation for every zone. 

# Returns
`nothing`. The `output` argument is modified in-place.

# Arguments
- `output`: A `Simulation{T}` object representing the output problem state. This argument is modified by this function.
- `input`: A `Simulation{T}` object representing the input problem state.

# Notes
- This function calls [`Energy()`](@ref) to perform the actual calculation. See the documentation of that function for more information.
- The artificial viscosity term is added to the pressure field during the call to [`Energy()`](@ref), per the artificial viscosity method of Von Neumann and Richtmyer (1950).
"""
function Energy!( output::Simulation{T}, input::Simulation{T} ) where { T <: AbstractFloat }
    for i in range( 1, input.nzones )
        output.energy_rhs[i] = Energy( input.mass[i], input.pressure[i] + input.viscosity[i], input.zone_length[i], input.velocity[i], input.velocity[i+1], input.energy_flux[i], input.energy_flux[i+1] )
    end
end
