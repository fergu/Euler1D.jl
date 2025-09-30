using Test
using Euler1D


# Now set up a simulation with a quiescent state. Nothing should happen, we just want to run for a known number of cycles
function TestCycleCallbacks()
    # These are two functions that increment a counter by one every time they are called
    global every_cycle_callback_count = 0
    function cycle_count_every_one( state::Simulation{T} ) where { T <: AbstractFloat }
        global every_cycle_callback_count += 1
    end

    global every_ten_cycles_callback_count = 0
    function cycle_count_every_ten( state::Simulation{T} ) where { T <: AbstractFloat }
        global every_ten_cycles_callback_count += 1
    end

    # Ratio of specific heats
    function γ( x::Float64 )
        return 1.4
    end
    
    # Density
    function ρ( x::Float64 )
        return 1.0
    end
    
    # Velocity
    function u( x::Float64 )
        return 0.0
    end
    
    # Pressure
    function P( x::Float64 )
        return 1.0
    end

    # Use as many default values as possible, only modifying parameters that are required for this specific case
    simulation_parameters = DefaultSimulationParameters() # Get a default set of parameters
    simulation_parameters["start_position"] = 0.0
    simulation_parameters["end_position"] = 1.0
    simulation_parameters["init_density_function"] = ρ # Set the function used to initialize density as the "ρ" function, above
    simulation_parameters["init_velocity_function"] = u # Set the function used to initialize velocity as the "u" function, above
    simulation_parameters["init_pressure_function"] = P # Set the function used to initialize pressure as the "P" function, above
    simulation_parameters["init_gamma_function"] = γ # Set the function used to initialize the ratio of specific heats as the "γ" function, above
    simulation_parameters["number_of_zones"] = 6000 # Set the number of zones in the problem
    init_state = InitializeSimulation( simulation_parameters ) # Initialize the simulation. This sets up the various arrays and such needed for the simulation and returns them in the `init_state` variable

    # Set up and register the callbacks
    callbacks = ConfigureSimulationCallbacks( init_state )
    RegisterCycleCallback!( callbacks, init_state, cycle_count_every_one, 1 )
    RegisterCycleCallback!( callbacks, init_state, cycle_count_every_ten, 10 )

    final_state = AdvanceNCycles( init_state, 1000; callbacks=callbacks ) # Advance 1000 cycles

    # First test that the final cycle count from the simulation state is equal to the number of times the callback was called
    # This tests that the callbacks are firing as often as expected
    @test final_state.cycles.x == every_cycle_callback_count

    # Second test that the final cycle count divided by 10 is equal to the number of times the second callback (which should only fire every 10 cycles) was called
    # This tests that callbacks that should only fire periodically are correctly doing so
    @test floor( final_state.cycles.x / 10 ) == every_ten_cycles_callback_count
end

TestCycleCallbacks()

