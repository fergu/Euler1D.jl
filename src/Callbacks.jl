function EvaluateCycleCallbacks( simulation::Simulation{T}, callbacks::Vector{CycleCallback} ) where { T <: AbstractFloat }
    for cb in callbacks
        if ( simulation.cycles.x >= cb.last_called.x + cb.every )
            # Call the callback
            cb.func( simulation )
            # Update the last called cycle
            cb.last_called.x = simulation.cycles.x
        end
    end
end

function EvaluateTimeCallbacks( simulation::Simulation{T}, callbacks::Vector{TimeCallback} ) where { T <: AbstractFloat }
    for cb in callbacks
        if ( simulation.time.x >= cb.times[ cb.next_index.x ] )
            # Call the callback
            cb.func( simulation )
            # Update the last called cycle
            cb.next_index.x = cb.next_index.x + 1
        end
    end
end

function EvaluateTimeDeltaCallbacks( simulation::Simulation{T}, callbacks::Vector{TimeDeltaCallback} ) where { T <: AbstractFloat }
    for cb in callbacks
        if ( simulation.time.x >= cb.last_called.x + cb.dt )
            # Call the callback
            cb.func( simulation )
            # Update the last called cycle
            cb.last_called.x = cb.last_called.x + cb.dt
        end
    end
end

function EvaluateCallbacks( simulation::Simulation{T}, callbacks::SimulationCallback ) where { T <: AbstractFloat }
    EvaluateCycleCallbacks( simulation, callbacks.callback_cycle )
    EvaluateTimeCallbacks( simulation, callbacks.callback_time )
    EvaluateTimeDeltaCallbacks( simulation, callbacks.callback_dt )
end

"""
    ConfigureSimulationCallbacks( simulation::Simulation{T} ) where { T <: AbstractFloat }

Configure a [`SimulationCallback`](@ref) structure for setting up callbacks as part of a simulation run.

# Returns
- A [`SimulationCallback`](@ref) structure

# Parameters
- `simulation`: A `Simulation{T}` describing the simulation that the callbacks are intended for.

# Notes
- Currently, the `simulation` parameter has no effect in this function and the resulting [`SimulationCallback`](@ref) structure can be used by any simulation, not just the one passed as an argument. However, this may change in future version.
- Callback entries can be added using [`RegisterCycleCallback!()`](@ref), [`RegisterTimeCallback!()`](@ref), or [`RegisterTimeDeltaCallback!()`](@ref).
"""
function ConfigureSimulationCallbacks( simulation::Simulation{T} ) where { T <: AbstractFloat }
    # Functions to set up callbacks
    callbacks = SimulationCallback(
                                Vector{CycleCallback}(undef,0),         # Vector of callback functions that are called after every N cycles
                                Vector{TimeCallback}(undef,0),          # Vector of callback functions that are called at a fixed list of times
                                Vector{TimeDeltaCallback}(undef,0)      # Vector of callback functions that are called at a fixed frequency
                               )
end

"""
    RegisterCycleCallback!( callbacks::SimulationCallback, func::Function, N::UInt; initial_cycle::UInt=0 )

Register a function `func` to be executed every `N` cycles of the simulation starting at cycle `initial_cycle`. Callbacks are stored in the `callbacks` input parameter and passed as an argument to the timestepping routines

# Returns
`nothing`, mutates the `callbacks` parameter

# Parameters
- `callbacks`: A `SimulationCallback` that describes the callbacks to be executed
- `func`: A `Function` that should be called by this callback. See Notes for information on expected function signature
- `N`: A `UInt` indicating how many cycles should occur between callbacks. Must be N ≥ 1
- `initial_cycle`: A `UInt` indicating the cycle number to start calling this callback at (Default: 0)

# Notes
- `func` is expected to accept arguments of the form `CallbackFunction( arg::Simulation{T} ) where { T <: AbstractFloat }`. `arg` will be the simulation state at the end of the cycle that the callback is executed on.
- This type of callback is useful for cases where a callback should be called at a fixed number of cycles between each callback.
"""
function RegisterCycleCallback!( callbacks::SimulationCallback, func::Function, N::UInt, initial_cycle::UInt=0 )
    # Check to be sure the requested callback frequency is at most once per cycle (i.e., not zero or negative)
    if ( N < 1 )
        throw("Cycle callback frequency must be ≥ 1 cycle")
    end

    # Create the callback
    #   Subtract N from initial_cycle to ensure the callback is first called at initial_cycle
    cb = CycleCallback( func, N, Ref( initial_cycle - N ) )

    # Add it to the list of callbacks
    push!( callbacks.callback_cycle, cb )
end

# These are just a few overloaded calls to catch cases where either N or initial_cycle are given as signed integer arguments
RegisterCycleCallback!( callbacks::SimulationCallback, func::Function, N::Int, initial_cycle::Union{Int,UInt}=0 ) = RegisterCycleCallback!( callbacks, func, convert(UInt, N ), convert(UInt, initial_cycle) )

"""
    RegisterTimeCallback!( callbacks::SimulationCallback, func::Function, times::Vector{T} ) where { T <: AbstractFloat }

Register a function `func` to be executed at the list of times specified in `times`. Callbacks are stored in the `callbacks` input parameter and passed as an argument to the timestepping routines

# Returns
`nothing`, mutates the `callbacks` parameter

# Parameters
- `callbacks`: A `SimulationCallback` that describes the callbacks to be executed
- `func`: A `Function` that should be called by this callback. See Notes for information on expected function signature
- `times`: A `Vector{T}` indicating the times at which to execute the callback

# Notes
- `func` is expected to accept arguments of the form `CallbackFunction( arg::Simulation{T} ) where { T <: AbstractFloat }`. `arg` will be the simulation state at the end of the cycle that the callback is executed on.
- Callbacks are executed at the end of the first cycle where the simulation time is greater than the next element in `times`. As a result, it is not guaranteed that the callback will be called at exactly the time specified in the `times` vector.
- `times` is sorted into ascending order before being stored. 
- This type of callback is useful if the callback should be called at an irregular series of times. See [`CycleCallback`](@ref) if the callback should be called at a regular number of cycles, or [`TimeDeltaCallback`](@ref) if the callback should be called at a fixed temporal cadence.
"""
function RegisterTimeCallback!( callbacks::SimulationCallback, func::Function, times::Vector{T} ) where { T <: AbstractFloat }
    # Create a copy of `times` and sort it into ascending order
    sorted_times = deepcopy( times )
    sorted_times = sort( times )
    # Append Infinity (Inf) to the end. 
    # We do this because the check to determine if the callback should be called works by comparing the value of `times` at the next index to be called against the current simulation time
    # If the current simulation time > time[index], the callback is run and `index` is set to index+1
    # Therefore, if we don't append a dummy element to the end of the `times` array, we will read off the end of the array after the last element in `times` is read
    # Appending ∞ specifically ensures that there is no chance this dummy element will ever trigger the callback
    push!( sorted_times, convert(T, Inf) )

    # Create the callback
    cb = TimeCallback( func, sorted_times, Ref( UInt(1) ) )

    # Add it to the list of callbacks
    push!( callbacks.callback_time, cb )
end

"""
    RegisterTimeDeltaCallback!( callbacks::SimulationCallback, func::Function, delta::T, initial_time::T=T(0.0) ) where { T <: AbstractFloat }

Register a function `func` to be executed every `delta` seconds starting at `initial_time`. Callbacks are stored in the `callbacks` input parameter and passed as an argument to the timestepping routines

# Returns
`nothing`, mutates the `callbacks` parameter

# Parameters
- `callbacks`: A `SimulationCallback` that describes the callbacks to be executed
- `func`: A `Function` that should be called by this callback. See Notes for information on expected function signature
- `delta`: A scalar `T` indicating how often to execute the callback
- `inital_time`: A scalar `T` indicating the time of the first callback to execute. (Default: 0.0)

# Notes
- `func` is expected to accept arguments of the form `CallbackFunction( arg::Simulation{T} ) where { T <: AbstractFloat }`. `arg` will be the simulation state at the end of the cycle that the callback is executed on.
- Callbacks are executed at the end of the first cycle where the simulation time is greater than the last time the callback was executed plus `delta`. As a result, it is not guaranteed that the callback will be called at exactly `delta` seconds since the last call.
- The time of the next callback is computed relative to the expected time of the current callback, _not_ when the callback was actually called. In other words, if a callback that should be called at `t₀` was actually called at `t₀+ϵ` due to the timestep size, the next callback will be scheduled for `t₀+δ`, *not* `t₀+ϵ+δ`.
- This type of callback is useful for cases where a callback should be called at a fixed temporal spacing between each callback.
"""
function RegisterTimeDeltaCallback!( callbacks::SimulationCallback, func::Function, delta::T, initial_time::T=T(0.0) ) where { T <: AbstractFloat }
    # Create our time delta callback
    #   Subtract `delta` from `init` to ensure the callback is first called at `init`
    #   This is required because the criteria used to decide to executed a callback is whether t > t_last + dt
    #   A small factor epsilon is also subtracted to ensture the greater than test passes at t = init
    cb = TimeDeltaCallback( func, delta, Ref( initial_time - delta - eps(T) ) )
    
    # Add it to the list of callbacks
    push!( callbacks.callback_dt, cb )
end

export ConfigureSimulationCallbacks
export RegisterCycleCallback!
export RegisterTimeCallback!
export RegisterTimeDeltaCallback!
