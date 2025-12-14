using Test
using Euler1D

# Now set up a simulation with a quiescent state. Nothing should happen, we just want to run for a known amount of time
function TestTimeDeltaCallbacks()
    # These are two functions that increment a counter by one every time they are called
    global dt_001_count = 0
    function dt_001_test( state::Simulation{T} ) where { T <: AbstractFloat }
        global dt_001_count += 1
    end

    global dt_005_count = 0
    function dt_005_test( state::Simulation{T} ) where { T <: AbstractFloat }
        global dt_005_count += 1
    end

    global dt_001_init_05_count = 0
    function dt_001_init_05_test( state::Simulation{T} ) where { T <: AbstractFloat }
        global dt_001_init_05_count += 1
    end

    global dt_005_init_05_count = 0
    function dt_005_init_05_test( state::Simulation{T} ) where { T <: AbstractFloat }
        global dt_005_init_05_count += 1
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

    # A test that should be called every 0.01 seconds starting at 0.0, for 101 total calls (the call at t=0 adds an extra 1 call)
    RegisterTimeDeltaCallback!( callbacks, dt_001_test, 0.01 )
    # A test that should be called every 0.05 seconds starting at 0.0, for 21 total calls (the call at t=0 adds an extra 1 call)
    RegisterTimeDeltaCallback!( callbacks, dt_005_test, 0.05 )
    # A test that should be called every 0.01 seconds starting at 0.5, for 51 total calls (the call at t=0.5 adds an extra 1 call)
    RegisterTimeDeltaCallback!( callbacks, dt_001_init_05_test, 0.01, 0.5 )
    # A test that should be called every 0.05 seconds starting at 0.5, for 11 total calls (the call at t=0.5 adds an extra 1 call)
    RegisterTimeDeltaCallback!( callbacks, dt_005_init_05_test, 0.05, 0.5 )

    final_state = AdvanceToTime( init_state, 1.0; callbacks=callbacks ) # Advance to t=1.0

    # Now test each callback to make sure they were called the correct number of times
    @test dt_001_count == 101 # Called every 0.01 seconds starting at 0.0 == 101 calls
    @test dt_005_count == 21  # Called every 0.05 seconds starting at 0.0 == 21 calls
    @test dt_001_init_05_count == 51 # Called every 0.01 seconds starting at 0.5 == 51 calls
    @test dt_005_init_05_count == 11 # Called every 0.05 seconds starting at 0.5 == 11 calls
end

TestTimeDeltaCallbacks()

