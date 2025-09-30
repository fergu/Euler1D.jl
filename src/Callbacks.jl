function EvaluateCycleCallbacks( simulation::Simulation{T}, callbacks::Vector{CycleCallback} ) where { T <: AbstractFloat }
    # Per-cycle callbacks
    # FIXME: This currently misbehaves. Since the output vector becomes the input vector on the next step for any of the timestepping routines that pre-allocate storage
    #       The changes made here become inputs on the next step - in other words, callbacks end up happening twice
    #       e.g., when cb.last_called.x is updated from cycle 0 to cycle 10, the updated struct then becomes an input on the next cycle, 
    #       but since we call callbacks based on the *output* state, we end up using old data (because the output state on the next cycle was the input state on this one, which was not updated)
    #       This can probably be addressed by somehow updating the output structs based on input data, but I need to think about it
    #
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
    function RegisterCycleCallback!( input::Simulation{T}, func::Function, N::UInt ) where { T <: AbstractFloat }

Register a function `func` to be executed every `N` cycles of the simulation defined by `input`

# Returns
`nothing`, mutates the `input` parameter

# Parameters
- `input` : A `Simulation{T}` to register the callback to
- `func` : A `Function` that should be called by this callback. See Notes for information on expected function signature
- `N` : A `UInt` indicating how many cycles should occur between callbacks. Must be N ≥ 1

# Notes
- `func` is expected to accept arguments of the form `CallbackFunction( arg::Simulation{T} ) where { T <: AbstractFloat }`. `arg` will be the simulation state at the end of the cycle that the callback is executed on.
"""
function RegisterCycleCallback!( callbacks::SimulationCallback, simulation::Simulation{T}, func::Function, N::UInt ) where { T <: AbstractFloat }
    # Check to be sure the requested callback frequency is at most once per cycle (i.e., not zero or negative)
    if ( N < 1 )
        throw("Cycle callback frequency must be ≥ 1 cycle")
    end

    # Create the callback
    cb = CycleCallback( func, N, Ref( simulation.cycles.x ) )

    # Add it to the list of callbacks
    push!( callbacks.callback_cycle, cb )
end

RegisterCycleCallback!( callbacks::SimulationCallback, simulation::Simulation{T}, func::Function, N::Int ) where { T <: AbstractFloat } = RegisterCycleCallback!( callbacks, simulation, func, UInt( N ) )

export ConfigureSimulationCallbacks
export RegisterCycleCallback!
