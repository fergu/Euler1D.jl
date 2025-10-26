using Test
using Euler1D

# Now set up a simulation with a quiescent state. Nothing should happen, we just want to run for a known amount of time
function TestTimeCallbacks()
    # These are two functions that increment a counter by one every time they are called
    global ordered_callback_count = 0
    function ordered_time_test( state::Simulation{T} ) where { T <: AbstractFloat }
        global ordered_callback_count += 1
    end

    global unordered_callback_count = 0
    function unordered_time_test( state::Simulation{T} ) where { T <: AbstractFloat }
        global unordered_callback_count += 1
    end

    global callall_callback_count = 0
    global callall_actual_callback_times = Vector{Float64}()
    function callall_time_test( state::Simulation{T} ) where { T <: AbstractFloat }
        push!( callall_actual_callback_times, state.time.x )
        global callall_callback_count += 1
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
    # These are two vectors of callback times - one initially in order and one out of order
    callback_times_sorted = Vector{Float64}([0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 1.1])
    callback_times_random = Vector{Float64}([1.1, 0.1, 0.25, 0.05, 0.5, 0.9, 0.75]) 
    # Also a vector where all elements should be called
    callback_times_allcall = Vector{Float64}([0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95])

    RegisterTimeCallback!( callbacks, ordered_time_test, callback_times_sorted )
    RegisterTimeCallback!( callbacks, unordered_time_test, callback_times_random )
    RegisterTimeCallback!( callbacks, callall_time_test, callback_times_allcall )

    final_state = AdvanceToTime( init_state, 1.0; callbacks=callbacks ) # Advance 1000 cycles

    # First, test that the ordered callback was called the correct numbr of times, which should be the length of the vector minus 1
    @test ordered_callback_count == length( callback_times_sorted ) - 1

    # Second, test that the random callback list was called the correct number of times, which should be the length of the vector minus 1
    @test unordered_callback_count == length( callback_times_random ) - 1

    # Third, check that the callback that should have called every element in the vector was called the right number of times
    @test callall_callback_count == length( callback_times_allcall )

    # Also check that all elements that were stored in the callall_actual_callback_times vector are approximately equal to what we wanted
    # We do not expect exactly equal as we do not try to control the timestep to match the requested time, so exactly how close we get depends on the simulation state
    # The chosen absolute tolerance of 1e-4 (4 digits) was determined heuristically to be about as precise as possible for this quiescent test case which will probably have an abnormally large timestep
    # Realistic cases may be more precise than this, but at least this test should be sensitive to any unexpected changes
    @test all( isapprox.( callback_times_allcall, callall_actual_callback_times; atol=1e-4 ) )
end

TestTimeCallbacks()

