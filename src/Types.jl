# These are some structs to define callback functions
abstract type AbstractSimulationCallback end

"""
    struct CycleCallback

A structure containing information about a callback intended to be executed every `every` cycles

# Parameters
- `func`: The `Function` to be executed
- `every`: How often (in cycles) to execute this callback
- `last_called`: The last cycle that this callback was called

# Notes
- This type of callback is useful for cases where a callback should be called at a fixed number of cycles between each callback.
- While this structure can be initialized directly and added to the `callback_cycle` entry of a [`SimulationCallback`](@ref) structure, it is recommended to call [`RegisterCycleCallback!()`](@ref) instead.
"""
struct CycleCallback{F} <: AbstractSimulationCallback
    func::F
    every::UInt
    last_called::Base.RefValue{UInt}
end

"""
    struct TimeCallback

A structure containing information about a callback intended to be executed at a fixed list of times given by `times`.

# Parameters
- `func`: The `Function` to be executed
- `times`: The list of times at which to execute this callback
- `next_index`: An index into the `times` vector indicating the next time the callback should be executed.

# Notes
- This type of callback is useful if the callback should be called at an irregular series of times. See [`CycleCallback`](@ref) if the callback should be called at a regular number of cycles, or [`TimeDeltaCallback`](@ref) if the callback should be called at a fixed temporal cadence.
- While this structure can be initialized directly and added to the `callback_time` entry of a [`SimulationCallback`](@ref) structure, it is recommended to call [`RegisterTimeCallback!()`](@ref) instead.
- The entries in `times` are assumed to be sorted in ascending order. [`RegisterTimeCallback!()`](@ref) will handle this automatically, but this will need to be handled manually if creating this structure directly.
"""
struct TimeCallback{F,T} <: AbstractSimulationCallback
    func::F
    times::Vector{T}
    next_index::Base.RefValue{UInt}
end

"""
    struct TimeDeltaCallback

A structure containing information about a callback intended to be executed every `every` seconds

# Parameters
- `func`: The `Function` to be executed
- `every`: How often (in seconds) to execute this callback
- `last_called`: The last time that this callback was called

# Notes
- This type of callback is useful for cases where a callback should be called at a fixed temporal spacing between each callback.
- While this structure can be initialized directly and added to the `callback_dt` entry of a [`SimulationCallback`](@ref) structure, it is recommended to call [`RegisterTimeDeltaCallback!()`](@ref) instead.
"""
struct TimeDeltaCallback{F,T} <: AbstractSimulationCallback
    func::F
    dt::T
    last_called::Base.RefValue{T}
end

"""
    struct SimulationCallback

A structure containing callbacks to be called by a simulation.

# Parameters
- `callback_cycle`: A `Vector` of [`CycleCallback`](@ref)'s that are called based on the number of cycles (timesteps) that have been performed
- `callback_time`: A `Vector` of [`TimeCallback`](@ref)'s that are called at a fixed set of times
- `callback_dt`: A `Vector` of [`TimeDeltaCallback`](@ref)'s that are called at a fixed temporal frequency

# Notes
- This structure is typically initialized using [`ConfigureSimulationCallbacks()`](@ref). See the documentation of that function for further detail.
"""
struct SimulationCallback
    # Vectors of callback functions to be called at various points during the simulation
    callback_cycle::Vector{CycleCallback}         # A vector of functions to be called every N cycles
    callback_time::Vector{TimeCallback}           # A vector of functions to be called at a fixed set of times
    callback_dt::Vector{TimeDeltaCallback}        # A vector of functions to be called at fixed time increments
end

"""
    struct Simulation{T}

A structure containing all the internal variables and arrays used in the simulation.

# Parameters
- `nzones::Int`: The number of zones in the simulation. (Unit: ⋅)
- `nedges::Int`: The number of zone edges in the simulation, equal to `nzones + 1`. (Unit: ⋅)
- `CFL::Float64`: The CFL number to be used when calculating timesteps. (Unit: ⋅)
- `start_time::Float64`: The initial time of the simulation. (Unit: s)
- `viscosity_coefficient::Float64`: The coefficient used to scale artificial viscosity. (Unit: ⋅)
- `conductivity_coefficient::Float64`: The coefficient used to scale artificial conductivity. (Unit: ⋅)
- `time::Base.RefValue{T}`: The current time of the simulation. (Unit: s)
- `dt::Base.RefValue{T}`: The size of the timestep taken in the last cycle. (Unit: s)
- `cycles::Base.RefValue{UInt}`: The number of cycles performed so far. (Unit: ⋅)
- `min_dt::Float64`: The minimum allowable timestep size. Simulation will halt if `Δt` falls below this value. (Unit: s)
- `max_cycles::UInt`: The maximum number of cycles to perform. Simulation will halt if `cycles` exceeds this value. (Unit: ⋅)
- `zone_edge::Vector{T}`: A vector of locations of zone edges. Increases monotonically. (Unit: m)
- `zone_center::Vector{T}`: A vector of locations of zone centers, defined as the midpoint between two zone edges. Increases monotonically. (Unit: m)
- `zone_length::Vector{T}`: A vector of the length of each zone. (Unit: m)
- `gamma::Vector{T}`: A vector of the ratio of specific heats inside each zone. (Unit: ⋅)
- `mass::Vector{T}`: A vector of the mass contained within each zone. Assumed constant. (Unit: kg)
- `density::Vector{T}`: A vector of the density of the fluid within each zone. Computed using the equation of state as mass/Δx. (Unit: kg/m³)
- `velocity::Vector{T}`: A vector of velocities of zone edges. (Unit: m/s)
- `pressure::Vector{T}`: A vector of pressures inside each zone. Computed from the equation of state. (Unit: kg/(m⋅s²))
- `intenergy::Vector{T}`: A vector of the internal energy per unit mass of each zone. (Unit: m²/s²)
- `speedofsound::Vector{T}`: A vector of the speed of sound within each zone. Computed from the equation of state. (Unit: m/s)
- `viscosity::Vector{T}`: A vector of the artificial viscosity within each zone, added to the pressure field. See ArtificialDissipation.jl for more information. (Unit: kg/(m⋅s²))
- `energy_flux::Vector{T}`: A vector of fluxes of internal energy per unit mass across each zone boundary. Added to the energy equation as a diffusion term. See ArtificialDissipation.jl for more information. (Unit: m³/s³)
- `momentum_rhs::Vector{T}`: A vector of the right hand side of the momentum equation at each zone edge. (Unit: m/s²)
- `energy_rhs::Vector{T}`: A vector of the right hand side of the energy equation within each zone. (Unit: m²/s³)
"""
struct Simulation{T}
    # Simulation configuration information (Constant)
    nzones::UInt                        # Number of zones in the simulation [⋅]
    nedges::UInt                        # Number of zone edges in the simulation. Equal to nzones + 1
    CFL::Float64                        # CFL number of the simulation [⋅]
    start_time::Float64                 # Initial time of the simulation [s]
    # Artificial viscosity and conductivity coefficients
    viscosity_coefficient::Float64      # Coefficient for artificial viscosity
    conductivity_coefficient::Float64   # Coefficient for artificial conductivity
    # Simulation state information
    time::Base.RefValue{T}              # Current time [s]
    dt::Base.RefValue{T}                # The timestep size of the last cycle
    cycles::Base.RefValue{UInt}         # The number of cycles performed so far
    min_dt::T                           # The minimum allowable Δt. Simulation will halt if Δt goes below this value
    max_cycles::UInt                    # The maximum allowable number of cycles. Simulation will halt if cycles goes above this value
    # Zone Positions
    zone_edge::Vector{T}                # Position of the zone edges [m]
    zone_center::Vector{T}              # Position of the zone centers [m]
    zone_length::Vector{T}              # Length of the zones [m]
    # Fluid properties (constant as a function of time)
    gamma::Vector{T}                    # Ratio of specific heats of each zone [⋅]
    mass::Vector{T}                     # Mass of each zone [kg]
    # State variables (change as a function of time)
    density::Vector{T}                  # Density of each zone [kg/m³]
    velocity::Vector{T}                 # Velocity at each zone edge [m/s]
    pressure::Vector{T}                 # Pressure in each zone [N/m^2]
    intenergy::Vector{T}                # Internal energy of each zone
    speedofsound::Vector{T}             # Speed of sound within each zone [m/s]
    viscosity::Vector{T}                # Artificial viscosity to be applied to a given zone
    energy_flux::Vector{T}              # Energy flux across a given zone boundary [edge centered]
    # RHS variables. These are integrated to calculate the state variables, above
    momentum_rhs::Vector{T}             # RHS of the momentum equation, ∂u/∂t 
    energy_rhs::Vector{T}               # RHS of the energy equation, ∂e/∂t
end

function Base.show( io::IO, obj::Simulation{T} ) where T
    println("Euler1D simulation")
    println("\tStart time: ", obj.start_time)
    println("\tCurrent time: ", obj.time.x)
    println("\tLast Δt: ", obj.dt.x, " (Min: ", obj.min_dt, ")")
    println("\tCycle count: ", obj.cycles.x, " (Max: ", obj.max_cycles, ")")
    println("\tCFL number: ", obj.CFL)
    println("\tNumber of zones: ", obj.nzones)
    println("\tAvailable fields: ", join( string.( fieldnames( Simulation{T} ) ), ", " ) )
end


export Simulation
export SimulationCallback
export CycleCallback
export TimeCallback
export TimeDeltaCallback
