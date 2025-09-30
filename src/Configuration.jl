"""
    InitializeSimulation( parameters::Dict{String, Any} )

Initialize a simulation with the parameters provided by `parameters`. A default set of `parameters` can be returned from the [`DefaultSimulationParameters()`](@ref) function.

# Returns
A struct of type `Simulation{T}` that describes the current simulation state. The type `T` is determined by the type of the members of `parameters`.

# Notes
- This function performs some basic sanity checking on input parameters such as ensuring that `end_time > start_time` or `end_position > start_position`, and will raise an error if any of these checks fail.
- A warning will be raised if any keys in the dictionary are unused, but this will not halt execution.
- This function will not warn about any potential issues with problem configuration or potential stability issues.
"""
function InitializeSimulation( parameters::Dict{String, Any} )
    # We'll be using mutating functions on the parameters dictionary, so we make a copy to avoid mutating the user's dictionary
    internal_parameters = deepcopy( parameters ) 

    # Get the number of grid zones
    number_of_zones = zero(UInt)
    try
        input_number_of_zones = pop!( internal_parameters, "number_of_zones" )
        number_of_zones = UInt( input_number_of_zones )
    catch
        error("Failed to interpret input 'number_of_zones' as an unsigned integer. Check that supplied value ($input_number_of_zones) is an integer > 0" )
    end

    # Read the timestep and cycle limits
    min_Δt = pop!( internal_parameters, "minimum_timestep" )

    max_cycles = zero(UInt)
    try
        input_max_cycles = pop!( internal_parameters, "maximum_cycles" )
        max_cycles = UInt( input_max_cycles )
    catch
        error("Failed to interpret input 'max_cycles' as an unsigned integer. Check that the supplied value ($input_max_cycles) is an integer > 0" )
    end

    # Get the CFL number
    CFL = pop!( internal_parameters, "CFL" )

    # Get the artificial viscosity and conductivity coefficients
    artificial_viscosity = pop!( internal_parameters, "artificial_viscosity_coefficient" )
    artificial_conductivity = pop!( internal_parameters, "artificial_conductivity_coefficient" )

    # Generate the grid
    x₀ = pop!( internal_parameters, "start_position" )
    x₁ = pop!( internal_parameters, "end_position" )
    if ( x₁ ≤ x₀ )
        error("Parameter `end_position` (=$x₁) must be greater than `start_position` (=$x₀)")
    end
    x = collect( range( start=x₀, stop=x₁, length=number_of_zones+1 ) ) # Need to add one to number_of_zones to get number of edges
    Δx = diff( x )
    xₘ = 0.5 .* ( x[2:end] .+ x[1:end-1] ) # Position of zone centers. This will be passed to the initialization functions for zone-centered quantities

    # Get the starting and ending times
    start_time = pop!( internal_parameters, "start_time" )

    # Get the initial condition functions
    density = pop!( internal_parameters, "init_density_function" )
    velocity = pop!( internal_parameters, "init_velocity_function" )
    pressure = pop!( internal_parameters, "init_pressure_function" )
    gamma = pop!( internal_parameters, "init_gamma_function" )

    # Check that all of the supplied initialization functions exist (e.g., that they are not `nothing`)
    # Throw an error if one was not supplied since we won't be able to continue
    if ( any( [ density, velocity, pressure, gamma ] .== nothing ) )
        error("Some initialization functions were not supplied. Double check that you have supplied initialization functions for density, velocity, pressure, and the ratio of specific heats.")
    end

    # Now double check that all keys have been used and emit a warning if not to warn the user that one of their keys wasn't used 
    if ( !isempty( internal_parameters ) )
        @warn "Not all supplied parameters were used in problem initialization. Unused parameters: $(keys(internal_parameters))"
    end

    # Initialize the simulation struct
    simulation = Simulation( 
                            number_of_zones,                        # Number of zones
                            number_of_zones + 1,                    # Number of edges
                            CFL,                                    # CFL Number
                            start_time,                             # Starting time
                            artificial_viscosity,                   # Coefficient for artificial viscosity
                            artificial_conductivity,                # Coefficient for artificial conductivity
                            Ref( start_time ),                      # Current time, initialized to the initial time
                            Ref( min_Δt ),                          # Size of the last timestep, initialized to min_Δt
                            Ref( zero(UInt) ),                      # Number of cycles taken so far, initialized to 0
                            min_Δt,                                 # Minimum timestep size
                            max_cycles,                             # Maximum numer of cycles
                            x,                                      # Vector of zone edge positions
                            xₘ,                                     # Vector of zone center positions
                            Δx,                                     # Vector of zone lengths
                            zeros( number_of_zones ),               # Vector of ratios of specific heats (zone centered)
                            zeros( number_of_zones ),               # Vector of zone masses (zone centered)
                            zeros( number_of_zones ),               # Vector of zone densities (zone centered)
                            zeros( number_of_zones + 1 ),           # Vector of zone edge velocities (edge centered)
                            zeros( number_of_zones ),               # Vector of zone pressures (zone centered)
                            zeros( number_of_zones ),               # Vector of zone internal energy (zone centered)
                            zeros( number_of_zones ),               # Vector of zone speeds of sound (zone centered)
                            zeros( number_of_zones ),               # Vector of artificial viscosities (zone centered)
                            zeros( number_of_zones + 1 ),           # Vector of artificial energy conduction (edge centered)
                            zeros( number_of_zones + 1 ),           # Vector of the acceleration of zone edges (edge centered)
                            zeros( number_of_zones )                # Vector of the time rate of change of internal energy of a zone (zone centered)
                           )

    # Now initialize the simulation fields using the supplied functions
    simulation.gamma .= gamma.( xₘ ) # Ratio of specific heats for each zone
    simulation.mass .= density.( xₘ ) .* Δx # Mass of each zone. Computed from the density function and the initial zone sizes
    simulation.velocity .= velocity.( x ) # Initial velocity
    simulation.intenergy .= pressure.( xₘ ) ./ ( ( gamma.( xₘ ) .- 1.0 ) .* density.( xₘ ) ) # Compute the initial internal energy assuming an ideal gas EOS
    EquationOfState!( simulation ) # Use the equation of state to update 

    return simulation
end

"""
    UpdateSimulationState!( state::Simulation{T}, gamma::Function, density::Function, velocity::Function, pressure::Function ) where { T <: AbstractFloat }

Update the simulation state using the supplied functions to define a new state. See Notes for information about the expected signature of the functions.

# Returns
`nothing`. Updates the `state` input in-place.

# Arguments
- `state`: A `Simulation{T}` representing the simulation state to be updated
- `gamma`: A `Function` that returns the new values of the ratio of specific heats, gamma
- `density`: A `Function` that returns the new values of the density
- `velocity`: A `Function` that returns the new values of velocity
- `pressure`: A `Function` that returns the new values of pressure

# Notes
Unlike the functions used in initial problem setup, the functions supplied to `UpdateSimulationState` have a slightly different expected signature of
    ExampleFunction( x::T, oldValue::T ) where { T <: AbstractFloat }
where `oldValue` will be the current value of the state variable at current position `x` (e.g., the previous value of pressure).

# Side Effects
- Updates the values stored in the vectors for `state.gamma`, `state.mass`, `state.velocity`, `state.intenergy`, and fields derived from the equation of state in-place.

!!! warning
    If the `density` function alters the `state.density` field (that is, it does not just return `oldValue`), mass will be added or removed from a given zone to achieve a specified value of `density` without changing zone size. As a consequence, mass will not be conserved in the system in this case.
"""
function UpdateSimulationState!( state::Simulation{T}, gamma::Function, density::Function, velocity::Function, pressure::Function ) where { T <: AbstractFloat }
    new_density = density.( state.zone_center, state.density )
    new_gamma = gamma.( state.zone_center, state.gamma )
    state.gamma .= new_gamma
    state.mass .= new_density .* state.zone_length
    state.velocity .= velocity.( state.zone_edge, state.velocity )
    state.intenergy .= pressure.( state.zone_center, state.pressure ) ./ ( ( new_gamma .- 1.0 ) .* new_density )
    EquationOfState!( state )
end

"""
    DefaultSimulationParameters()

Return a set of default parameters for a simulation.

# Returns 
A `Dict{String, Any}` with default parameters for a simulation. 

# Notes
- Keys with a value of `nothing` must be set before using this dictionary to initialize a simulation. This is done instead of providing a default function to avoid the possibility of accidentally initializing a simulation with unintended initial conditions.

# Parameters
- `start_time`: The initial time of the simulation. (Unit: s; Default: 0.0)
- `end_time`: The final time of the simulation. (Unit: s; Default: 1.0)
- `start_position`: The position of the left side of the domain. (Unit: m; Default: 0.0)
- `end_position`: The position of the right side of the domain. (Unit: m; Default: 1.0)
- `number_of_zones`: The number of zones to divide the domain into. (Unit: ⋅; Default: 1000)
- `CFL`: The CFL number to use. (Unit: ⋅; Default: 0.2)
- `artificial_viscosity_coefficient`: The coefficient to use for artificial viscosity. (Unit: ⋅; Default: 1e0)
- `artificial_conductivity_coefficient`: The coefficient to use for artificial conductivity. (Unit: ⋅; Default: 1e-2)
- `min_timestep`: The minimum allowable timestep size. Simulation will halt if timestep falls below this value. (Unit: s; Default: 1e-7)
- `max_cycles`: The maximum number of cycles to perform. Simulation will halt if more than this many timesteps are taken. (Unit: ⋅; Default: 1e6)
- `init_density_function`: A `Function` that returns the initial density as a function of position `x`. (Unit: kg/m³; Default: `nothing`, must be user-supplied)
- `init_velocity_function`: A `Function` that returns the initial velocity as a function of position `x`. (Unit: m/s; Default: `nothing`, must be user-supplied)
- `init_pressure_function`: A `Function` that returns the initial pressure as a function of position `x`. (Unit: N/m²; Default: `nothing`, must be user-supplied)
- `init_gamma_function`: A `Function` that returns the ratio of specific heats as a function of position `x`. (Unit: ⋅; Default: `nothing`, must be user-supplied)
"""
function DefaultSimulationParameters()
    return Dict{String, Any}(
                             "start_time" => 0.0,
                             "start_position" => 0.0,
                             "end_position" => 1.0,
                             "number_of_zones" => 1000,
                             "CFL" => 0.2,
                             "artificial_viscosity_coefficient" => 1e0,
                             "artificial_conductivity_coefficient" => 1e-2,
                             "minimum_timestep" => 1e-7,
                             "maximum_cycles" => 1e6,
                             "init_density_function" => nothing,
                             "init_velocity_function" => nothing,
                             "init_pressure_function" => nothing,
                             "init_gamma_function" => nothing,
                            )
end

export InitializeSimulation
export UpdateSimulationState!
export DefaultSimulationParameters
