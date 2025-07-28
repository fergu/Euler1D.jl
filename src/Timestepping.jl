"""
    AdvanceToTime( state::Simulation{T}, stoptime::T, Δt::T; exact::Bool=false ) where { T <: AbstractFloat }

Advance the simulation with initial state `state` to a time specified by `stoptime` with a fixed timestep. 

# Returns
A `Simulation{T}` representing the state at the end of the final cycle.

# Arguments
- `input`: A `Simulation{T}` representing the simulation state at the start of the first cycle.
- `stoptime`: The time to advance to. (Unit: s)
- `Δt`: The timestep size. (Unit: s)
- `exact`: If true, try to stop as close as possible to `stoptime` by adjusting the final timestep size (Default: `false`)

# Notes
- The current simulation time is determined by the `Simulation` field `state.time`. If `state.time > stoptime`, no steps will be taken.
- This function allocates two `deepcopy()`s of the input state and returns the copy corresponding to the final state.
- This function simply calls `AdvanceOneCycle!` repeatedly until the simulation time reaches `stoptime`. The primary advantage to using this function as opposed to [`AdvanceOneCycle()`](@ref) or [`AdvanceNCycles()`](@ref) is that various backing arrays are pre-allocated to improve speed.
- If `exact=true`, the timestep of the final cycle is adjusted so that the time of the final state is as close as possible to `stoptime`.
"""
function AdvanceToTime( state::Simulation{T}, stoptime::T, Δt::T; exact::Bool=false ) where { T <: AbstractFloat }
    # To start, we'll make two copies of the input state that we'll cycle between to avoid allocations in the actual numerical scheme
    # We need two copies as otherwise the algorithm will mutate the input state, which we don't want
    input = deepcopy(state)
    output = deepcopy(state)
    while ( output.time.x <= stoptime )
        # Swap the input and output structures
        input, output = output, input
        # Reduce Δt if we want to end at exactly the specified time
        if ( ( input.time.x + Δt > stoptime ) & ( exact == true ) )
            Δt = stoptime - input.time.x + eps(Float64) # We add machine epsilon to Δt just to ensure the less-than-or-equal-to check passes after this step
        end
        # Advance by one cycle
        AdvanceOneCycle!( output, input, Δt )
    end
    return output
end

"""
    AdvanceToTime( state::Simulation{T}, stoptime::T; exact::Bool=false ) where { T <: AbstractFloat }

Advance the simulation with initial state `state` to a time specified by `stoptime` with a variable timestep. 

# Returns
A `Simulation{T}` representing the state at the end of the final cycle

# Arguments
- `input`: A `Simulation{T}` representing the simulation state at the start of the first cycle.
- `stoptime`: The time to advance to. (Unit: s)
- `exact`: If true, try to stop as close as possible to `stoptime` by adjusting the final timestep size (Default: false)

# Notes
- The current simulation time is determined by the `Simulation` field `state.time`. If `state.time > stoptime`, no steps will be taken.
- This function allocates two `deepcopy()`s of the input state and returns the copy corresponding to the final state.
- This function simply calls `AdvanceOneCycle!()` repeatedly until the simulation time reaches `stoptime`. The primary advantage to using this function as opposed to [`AdvanceOneCycle()`](@ref) or [`AdvanceNCycles()`](@ref) is that various backing arrays are pre-allocated to improve speed.
- The timestep size, Δt, is determined for each cycle based on the minimum time for an acoustic wave to traverse a zone. See the documentation for [`CalculateTimestepSize()`](@ref) for further details. 
- If `exact=true`, the timestep of the final cycle is adjusted so that the time of the final state is as close as possible to `stoptime`.
"""
function AdvanceToTime( state::Simulation{T}, stoptime::T; exact::Bool=false ) where { T <: AbstractFloat }
    # To start, we'll make two copies of the input state that we'll cycle between to avoid allocations in the actual numerical scheme
    # We need two copies as otherwise the algorithm will mutate the input state, which we don't want
    input = deepcopy(state)
    output = deepcopy(state)
    while ( output.time.x <= stoptime )
        # Swap the input and output structures
        input, output = output, input
        # Get the timestep required by the CFL condition
        Δt = CalculateTimestepSize( input )
        # Reduce Δt if we want to end at exactly the specified time
        if ( ( input.time.x + Δt > stoptime ) & ( exact == true ) )
            Δt = stoptime - input.time.x + eps(Float64) # We add machine epsilon to Δt just to ensure the less-than-or-equal-to check passes after this step
        end
        # Advance by one cycle
        AdvanceOneCycle!( output, input, Δt )
    end
    return output
end

"""
    AdvanceOneCycle!( output::Simulation{T}, input::Simulation{T}, Δt::T ) where { T <: AbstractFloat }

Advance the simulation by one cycle with a timestep of Δt.

# Returns
`nothing`. Modifies `output` in-place.

# Arguments
- `output`: A `Simulation{T}` that will represent the output state. This will be modified by the function to represent the simulation state after advancing one cycle.
- `input`: A `Simulation{T}` that represents the simulation state at the start of the cycle.
- `Δt`: The size of the time step. (Unit: s)

# Side Effects
- All fields of `output` are modified in-place.
"""
function AdvanceOneCycle!( output::Simulation{T}, input::Simulation{T}, Δt::T ) where { T <: AbstractFloat }
    cycle!( output, input, Δt )
end

"""
    AdvanceOneCycle!( output::Simulation{T}, input::Simulation{T} ) where { T <: AbstractFloat }

Advance the simulation by one cycle. 

# Returns
`nothing`. Modifies `output` in-place.

# Arguments
- `output`: A `Simulation{T}` that will represent the output state. This will be modified by the function to represent the simulation state after advancing one cycle.
- `input`: A `Simulation{T}` representing the simulation state at the start of the cycle

# Notes
- The timestep size, Δt, is determined based on the minimum time for an acoustic wave to traverse a zone. See [`CalculateTimestepSize()`](@ref) for further details.

# Side Effects
- All fields of `output` are modified in-place.
"""
function AdvanceOneCycle!( output::Simulation{T}, input::Simulation{T} ) where { T <: AbstractFloat }
    Δt = CalculateTimestepSize( input )
    AdvanceOneCycle!( output, input, Δt )
end

"""
    AdvanceOneCycle( state::Simulation{T}, Δt::T ) where { T <: AbstractFloat }

Advance the simulation by one cycle with a timestep of Δt.

# Returns
A `Simulation{T}` representing the state at the end of the cycle

# Arguments
- `input`: A `Simulation{T}` representing the simulation state at the start of the cycle.
- `Δt`: The size of the time step. (Unit: s)

# Notes
- This function allocates a `deepcopy()` of the input state. The copy is modified and returned from this function.
"""
function AdvanceOneCycle( state::Simulation{T}, Δt::T ) where { T <: AbstractFloat }
    output = deepcopy( state )
    AdvanceOneCycle!( output, state, Δt )
    return output
end

"""
    AdvanceOneCycle( state::Simulation{T} ) where { T <: AbstractFloat }

Advance the simulation by one cycle. 

# Returns
- A `Simulation{T}` representing the state at the end of the cycle

# Arguments
- `input`: A `Simulation{T}` representing the simulation state at the start of the cycle

# Notes
- This function allocates a `deepcopy()` of the input state and returns the copy.
- The timestep size, Δt, is determined based on the minimum time for an acoustic wave to traverse a zone. See [`CalculateTimestepSize()`](@ref) for further details.
"""
function AdvanceOneCycle( state::Simulation{T} ) where { T <: AbstractFloat }
    output = deepcopy( state )
    AdvanceOneCycle!( output, state )
    return output
end

"""
    AdvanceNCycles( state::Simulation{T}, ncycles::UInt, Δt::T ) where { T <: AbstractFloat }

Advance the simulation by `ncycles` cycles with a fixed timestep.

# Returns
- A `Simulation{T}` representing the state at the end of the final cycle

# Arguments
- `input`: A `Simulation{T}` representing the simulation state at the start of the first cycle.
- `ncycles`: The number of cycles to advance.
- `Δt`: The size of the time step. (Unit: s)

# Notes
- This function allocates two `deepcopy()`s of the input state and returns the copy corresponding to the final state.
- This function calls [`AdvanceOneCycle!()`](@ref) a total of `ncycles` times to advance the simulation. The primary advantage to using this function as opposed to [`AdvanceOneCycle()`](@ref) if the number of cycles to advance is known is that various backing arrays are pre-allocated to improve speed.
"""
function AdvanceNCycles( state::Simulation{T}, ncycles::UInt, Δt::Float64 ) where { T <: AbstractFloat }
    # To start, we'll make two copies of the input state that we'll cycle between to avoid allocations in the actual numerical scheme
    # We need two copies as otherwise the algorithm will mutate the input state, which we don't want
    input = deepcopy(state)
    output = deepcopy(state)
    # Perform the requested number of cycles
    for cycle in range( 1, ncycles )
        # Swap the input and output structures
        input, output = output, input
        # Advance by one cycle
        AdvanceOneCycle!( output, input, Δt )
    end
    # Return the result
    return output
end

"""
    AdvanceNCycles( state::Simulation{T}, ncycles::UInt ) where { T <: AbstractFloat }

Advance the simulation by `ncycles` cycles with a variable timestep. 

# Returns
- A `Simulation{T}` representing the state at the end of the final cycle

# Arguments
- `input`: A `Simulation{T}` representing the simulation state at the start of the first cycle.
- `ncycles`: The number of cycles to advance. (Unit: ⋅)

# Notes
- This function allocates two `deepcopy()`s of the input state and returns the copy corresponding to the final state.
- This function calls [`AdvanceOneCycle!()`](@ref) a total of `ncycles` times to advance the simulation. The primary advantage to using this function as opposed to [`AdvanceOneCycle()`](@ref) if the number of cycles to advance is known is that various backing arrays are pre-allocated to improve speed.
- The timestep size, Δt, is determined based on the minimum time for an acoustic wave to traverse a zone. See [`CalculateTimestepSize()`](@ref) for further details.
"""
function AdvanceNCycles( state::Simulation{T}, ncycles::UInt ) where { T <: AbstractFloat }
    # To start, we'll make two copies of the input state that we'll cycle between to avoid allocations in the actual numerical scheme
    # We need two copies as otherwise the algorithm will mutate the input state, which we don't want
    input = deepcopy(state)
    output = deepcopy(state)
    # Perform the requested number of cycles
    for cycle in range( 1, ncycles )
        # Swap the input and output structures
        input, output = output, input
        # Get the timestep required by the CFL condition
        Δt = CalculateTimestepSize( input )
        # Advance by one cycle
        AdvanceOneCycle!( output, input, Δt )
    end
    # Return the result
    return output
end

"""
    CalculateTimestepSize( state::Simulation{T} ) where { T <: AbstractFloat }

Compute an automatic timestep size for the next simulation cycle based on the current simulation state. See the Notes section for details on how the timestep is determined.

# Returns
- A scalar of type `T` representing the timestep size for the next cycle based on the current simulation state.

# Arguments
- `state`: A `Simulation{T}` representing the problem state.

# Notes
- For each zone, the local speed of sound is computed according to `c = √( γ P / ρ )`, where γ, P, and ρ are the ratio of specific heats, the pressure, and the density of the gas in that zone.
- The time for an acoustic wave to traverse a zone with length `Δx` is computed as `t = Δx / c`.
- The minimum traversal time for all zones is multiplied by the user-specified CFL number to obtain the timestep size.
- Sanity checking for negative zone sizes and small timesteps is performed to detect problem instability.
"""
function CalculateTimestepSize( state::Simulation{T} ) where { T <: AbstractFloat }
    # CFL = soundspeed * dt / dx, therefore dt = CFL * dx / soundspeed
    # We therefore want the minimum value of dx / soundspeed
    minΔt = 1.0e6
    for i in range( 1, state.nzones )
        Δx = state.zone_length[i]
        if ( Δx <= 0.0 )
            error("Negative zone size detected at x=$(state.x[i]), t=$(state.time.x)")
        end
        c = state.speedofsound[i]
        Δt = Δx / c
        if ( Δt < minΔt )
            minΔt = Δx / c
        end
    end
    if ( state.CFL * minΔt < state.min_Δt )
        error("Δt ($(state.CFL * minΔt)) < minimum Δt ($(state.min_Δt)). Exiting.")
    end
    return state.CFL * minΔt
end

export AdvanceToTime
export AdvanceOneCycle, AdvanceOneCycle!
export AdvanceNCycles
export CalculateTimestepSize
