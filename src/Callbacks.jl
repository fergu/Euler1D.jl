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

function EvaluateCallbacks( simulation::Simulation{T}, callbacks::SimulationCallback ) where { T <: AbstractFloat }
    EvaluateCycleCallbacks( simulation, callbacks.callback_cycle )
end

function ConfigureSimulationCallbacks( simulation::Simulation{T} ) where { T <: AbstractFloat }
    # Functions to set up callbacks
    callbacks = SimulationCallback(
                                Vector{CycleCallback}(undef,0),         # Vector of callback functions that are called after every N cycles
                                Vector{TimeCallback}(undef,0),          # Vector of callback functions that are called at a fixed list of times
                                Vector{TimeDeltaCallback}(undef,0)      # Vector of callback functions that are called at a fixed frequency
                               )
end

"""
    RegisterCycleCallback!( callbacks::SimulationCallback, simulation::Simulation{T}, func::Function, N::UInt ) where { T <: AbstractFloat }

Register a function `func` to be executed every `N` cycles of the simulation defined by `simulation`. Callbacks are stored in the `callbacks` input parameter

# Returns
`nothing`, mutates the `callbacks` parameter

# Parameters
- `callbacks` : A `SimulationCallback` that describes the callbacks to be executed
- `simulation` : A `Simulation{T}` used to initialize the callback
- `func` : A `Function` that should be called by this callback. See Notes for information on expected function signature
- `N` : A `UInt` indicating how many cycles should occur between callbacks. Must be N ≥ 1

# Notes
- `func` is expected to accept arguments of the form `CallbackFunction( arg::Simulation{T} ) where { T <: AbstractFloat }`. `arg` will be the simulation state at the end of the cycle that the callback is executed on.
- The `simulation` argument is used to initalize the callback, primarily setting the in
"""
function RegisterCycleCallback!( callbacks::SimulationCallback, func::Function, N::UInt; simulation::Union{Simulation{T},Nothing}=nothing ) where { T <: AbstractFloat }
    # Check to be sure the requested callback frequency is at most once per cycle (i.e., not zero or negative)
    if ( N < 1 )
        throw("Cycle callback frequency must be ≥ 1 cycle")
    end

    # Determine the initial cycle for the callback
    # This is set to the current simulation cycle if `simulation` is not `nothing`, otherwise 0
    initial_cycle = isnothing( simulation ) ? UInt( 0 ) : simulation.cycles.x

    # Create the callback
    cb = CycleCallback( func, N, Ref( initial_cycle ) )

    # Add it to the list of callbacks
    push!( callbacks.callback_cycle, cb )
end

RegisterCycleCallback!( callbacks::SimulationCallback, func::Function, N::Int; simulation::Union{Simulation{T},Nothing}=nothing ) where { T <: AbstractFloat } = RegisterCycleCallback!( callbacks, func, UInt( N ); simulation=simulation )

export ConfigureSimulationCallbacks
export RegisterCycleCallback!
