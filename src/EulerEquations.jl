"""
    Momentum( ρ₋::T, P₋::T, Δx₋::T, ρ₊::T, P₊::T, Δx₊::T ) where { T <: AbstractFloat }

Compute the right hand side of the momentum equation at a cell interface. See the Notes section for information on the equation being solved.

# Returns
A scalar of type `T` representing the rate of change of velocity of a cell interface over time.

# Arguments
- ρ : The density of the cell. [ kg/m³ ]
- P : The pressure of the cell [ kg/(m⋅s²) ]
- Δx: The length of the cell. [ m ]

For each of these parameters, a - subscript refers to the cell to the left of the cell interface and a + subscript refers to the cell to the right of the cell interface.

# Notes
The governing equation solved in this function is:
    ∂u/∂t = (1/ρ) * ∂P/∂x
Here, ρ is taken to be the average of the densities in the cells to the left and right of the interface.
The pressure term is augmented with an artifical viscosity per the method of Von Neumann and Richtmyer (1950). This artificial viscosity is added in the `Momentum!()` function. See the documentation of that function for further details.
"""
function Momentum( ρ₋::T, P₋::T, Δx₋::T, ρ₊::T, P₊::T, Δx₊::T ) where { T <: AbstractFloat }
    ∂p∂x = ∂∂x_CellCenterToCellEdge( P₋, P₊, Δx₋, Δx₊ ) # Derivative of pressure, evaluated at the cell boundary
    ρ̄ = 0.5 * ( ρ₋ + ρ₊ ) # Compute an average density from the cells to the left and the right of the boundary
    return - ∂p∂x / ρ̄ 
end

"""
    Momentum!( output::Simulation{T}, input::Simulation{T} ) where { T <: AbstractFloat }

Compute the right hand side of the momentum equation for every cell interface. 

# Returns
`nothing`. The `output` argument is modified in-place.

# Arguments
- `output`: A `Simulation{T}` object representing the output problem state. This argument is modified by the function.
- `input` : A `Simulation{T}` object representing the input problem state.

# Notes
This function calls `Momentum()` to perform the actual calculation. See the documentation of that function for more information.
The artificial viscosity term is added to the pressure field during the call to `Momentum()`, per the artificial viscosity method of Von Neumann and Richtmyer (1950).

# Side Effects
- Modifies the values contained in `output.∂u∂t` in-place.
"""
function Momentum!( output::Simulation{T}, input::Simulation{T} ) where { T <: AbstractFloat }
    for i in range( 2, input.nedges - 1 )
        output.∂u∂t[i] = Momentum( input.density[i-1], input.pressure[i-1] + input.viscosity[i-1], input.Δx[i-1], input.density[i], input.pressure[i] + input.viscosity[i], input.Δx[i] )
    end
end

"""
    Energy( ρ::T, P::T, Δx::T, u₋::T, u₊::T, q₋::T, q₊::T ) where { T <: AbstractFloat }

Compute the right hand side of the energy equation in a cell. See the Notes section for the specific equations that are solved.

# Returns
A scalar of type `T` representing the rate of change of internal energy inside the cell.

# Arguments
- ρ : The density of the cell, [ kg/m³ ]
- P : The pressure of the cell, [ kg/(m⋅s²) ]
- Δx: The length of the cell, [ m ]
- u : The velocity of the cell edges on the (-): left and (+): right of the cell, [ m/s ]
- q : The (artificial) flux of internal energy across the (-): left and (+): right cell edges, [ m³/s³ ] 

# Notes:
The governing equation solved in this function is
    ∂e/∂t = - ( P / ρ ) * ∂u/∂x + ∑q
"""
function Energy( ρ::T, P::T, Δx::T, u₋::T, u₊::T, q₋::T, q₊::T ) where { T <: AbstractFloat }
    ∂u∂x = ∂∂x_CellEdgeToCellCenter( u₋, u₊, Δx )
    return -( P / ρ ) * ∂u∂x + ( q₋ - q₊ ) / Δx
end

"""
    Energy!( output::Simulation{T}, input::Simulation{T} ) where { T <: AbstractFloat }

Compute the right hand side of the energy equation for every cell. 

# Returns
`nothing`. The `output` argument is modified in-place.

# Arguments
- `output`: A `Simulation{T}` object representing the output problem state. This argument is modified by this function.
- `input` : A `Simulation{T}` object representing the input problem state.

# Notes
- This function calls `Energy()` to perform the actual calculation. See the documentation of that function for more information.
- The artificial viscosity term is added to the pressure field during the call to `Energy()`, per the artificial viscosity method of Von Neumann and Richtmyer (1950).
"""
function Energy!( output::Simulation{T}, input::Simulation{T} ) where { T <: AbstractFloat }
    for i in range( 1, input.nzones )
        output.∂e∂t[i] = Energy( input.density[i], input.pressure[i] + input.viscosity[i], input.Δx[i], input.velocity[i], input.velocity[i+1], input.energy_flux[i], input.energy_flux[i+1] )
    end
end
