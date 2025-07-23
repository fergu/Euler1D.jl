using Euler1D

### Functions to set up the 1D simulation

# Ratio of specific heats
function γ( x::Float64 )
    return 1.4
end

# Density
function ρ( x::Float64 )
    if ( x < 0.5 )
        return 1.0
    else
        return 0.125
    end
end

# Velocity
function u( x::Float64 )
    return 0.0
end

# Pressure
function P( x::Float64 )
    if ( x < 0.5 )
        return 1.0
    else
        return 0.1
    end
end

simulation_parameters = DefaultSimulationParameters() # Get a default set of parameters
simulation_parameters["end_time"] = 0.1 # Set the end time to 0.1 seconds
simulation_parameters["init_density_function"] = ρ # Set the function used to initialize density as the "ρ" function, above
simulation_parameters["init_velocity_function"] = u # Set the function used to initialize velocity as the "u" function, above
simulation_parameters["init_pressure_function"] = P # Set the function used to initialize pressure as the "P" function, above
simulation_parameters["init_gamma_function"] = γ # Set the function used to initialize the ratio of specific heats as the "γ" function, above
simulation_parameters["artificial_conductivity_coefficient"] = 1.0e-2 # Set the artificial conductivity coefficient. This value is the same as the default value, but is included here for demonstration
simulation_parameters["artificial_viscosity_coefficient"] = 1.0e0 # Set the artificial viscosity coefficient. This value is the same as the default value, but is included here for demonstration
simulation_parameters["number_of_zones"] = 6000 # Set the number of zones in the problem
init_state = InitializeSimulation( simulation_parameters ) # Initialize the simulation. This sets up the various arrays and such needed for the simulation and returns them in the `init_state` variable
final_state = AdvanceToTime( init_state, init_state.t₁; exact=true ) # Run the simulation to the final time

# At this point you can explore the final problem state
println( final_state )
# Or look at various arrays, etc
# For example, you could get the profile of velocity using
# println( final_state.velocity )
# Sicne velocity is located at zone edges, the corresponding x locations as
# println( final_state.zone_edge )
# You could then use Plots.jl to plot this data, DataFrames to save it to a file, and so on.
# Similarly, you can get the profile of density using
# println( final_state.density ) 
# Since density is zone-centered, the x positions are
# println( final_state.zone_center )
# And, as above, you can use these values to plot the data, save it to file, or so on
