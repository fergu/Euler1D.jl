"""
    struct Simulation{T}

A structure containing all the internal variables and arrays used in the simulation.

# Parameters
- nzones::Int: The number of cells in the simulation. [ ⋅ ]
- nedges::Int: The number of cell edges in the simulation, equal to nzones + 1. [ ⋅ ]
- CFL::Float64: The CFL number to be used when calculating timesteps. [ ⋅ ]
- t₀::Float64: The initial time of the simulation. [ s ]
- t₁::Float64: The final time of the simulation. [ s ]
- Cᵥ::Float64: The coefficient used to scale artificial viscosity. [ ⋅ ]
- Cₖ::Float64: The coefficient used to scale artificial conductivity. [ ⋅ ]
- time::Base.RefValue{T}: The current time of the simulation. [ s ]
- Δt::Base.RefValue{T}: The size of the timestep taken in the last cycle. [ s ]
- cycles::Base.RefValue{UInt}: The number of cycles performed so far. [ ⋅ ]
- min_Δt::Float64: The minimum allowable timestep size. Simulation will halt if `Δt` falls below this value. [ s ]
- max_cycles::UInt: The maximum number of cycles to perform. Simulation will halt if `cycles` exceeds this value. [ ⋅ ]
- x::Vector{T}: A vector of locations of cell edges. Increases monotonically. [ m ]
- Δx::Vector{T}: A vector of the size of each cell. [ m ]
- gamma::Vector{T}: A vector of the ratio of specific heats inside each cell. [ ⋅ ]
- mass::Vector{T}: A vector of the mass contained within each cell. Assumed constant. [ kg ]
- density::Vector{T}: A vector of the density of the fluid within each cell. Computed using the equation of state as mass/Δx. [ kg/m³ ]
- velocity::Vector{T}: A vector of velocities of cell edges. [ m/s ]
- pressure::Vector{T}: A vector of pressures inside each cell. Computed from the equation of state. [kg/(m⋅s²)]
- intenergy::Vector{T}: A vector of the internal energy per unit mass of each cell. [ m²/s² ]
- speedofsound::Vector{T}: A vector of the speed of sound within each cell. Computed from the equation of state. [ m/s ]
- viscosity::Vector{T}: A vector of the artificial viscosity within each cell, added to the pressure field. See ArtificialDissipation.jl for more information. [ kg/(m⋅s²) ]
- energy_flux::Vector{T}: A vector of fluxes of internal energy per unit mass across each cell boundary. Added to the energy equation as a diffusion term. See ArtificialDissipation.jl for more information. [ m³/s³ ]
- ∂u∂t::Vector{T}: A vector of the right hand side of the momentum equation at each cell edge. [ m/s² ]
- ∂e∂t::Vector{T}: A vector of the right hand side of the energy equation within each cell. [ m²/s³ ]
"""
struct Simulation{T}
    # Simulation configuration information (Constant)
    nzones::UInt                        # Number of cells in the simulation [⋅]
    nedges::UInt                        # Number of cell edges in the simulation. Equal to nzones + 1
    CFL::Float64                        # CFL number of the simulation [⋅]
    t₀::Float64                         # Initial time of the simulation [s]
    t₁::Float64                         # End time of the simulation [s]
    # Artificial viscosity and conductivity coefficients
    Cᵥ::Float64                         # Coefficient for artificial viscosity
    Cₖ::Float64                         # Coefficient for artificial conductivity
    # Simulation state information
    time::Base.RefValue{T}              # Current time [s]
    Δt::Base.RefValue{T}                # The timestep size of the last cycle
    cycles::Base.RefValue{UInt}         # The number of cycles performed so far
    min_Δt::T                           # The minimum allowable Δt. Simulation will halt if Δt goes below this value
    max_cycles::UInt                    # The maximum allowable number of cycles. Simulation will halt if cycles goes above this value
    # Cell Positions
    x::Vector{T}                        # Position of the cell edges [m]
    Δx::Vector{T}                       # Spacing between the cells [m]
    # Fluid properties (constant as a function of time)
    gamma::Vector{T}                    # Ratio of specific heats of each cell [⋅]
    mass::Vector{T}                     # Mass of each cell [kg]
    # State variables (change as a function of time)
    density::Vector{T}                  # Density of each cell [kg/m³]
    velocity::Vector{T}                 # Velocity at each cell edge [m/s]
    pressure::Vector{T}                 # Pressure in each cell [N/m^2]
    intenergy::Vector{T}                # Internal energy of each cell
    speedofsound::Vector{T}             # Speed of sound within each cell [m/s]
    viscosity::Vector{T}                # Artificial viscosity to be applied to a given cell
    energy_flux::Vector{T}              # Energy flux across a given cell boundary [edge centered]
    # RHS variables. These are integrated to calculate the state variables, above
    ∂u∂t::Vector{T}                     # RHS of the momentum equation, ∂u/∂t 
    ∂e∂t::Vector{T}                     # RHS of the energy equation, ∂e/∂t
end

function Base.show( io::IO, obj::Simulation{T} ) where T
    println("Euler1D simulation")
    println("\tStart time: ", obj.t₀)
    println("\tCurrent time: ", obj.time.x)
    println("\tFinal time: ", obj.t₁)
    println("\tLast Δt: ", obj.Δt.x, " (Min: ", obj.min_Δt, ")")
    println("\tCycle count: ", obj.cycles.x, " (Max: ", obj.max_cycles, ")")
    println("\tCFL number: ", obj.CFL)
    println("\tNumber of zones: ", obj.nzones)
    println("\tAvailable fields: ", join( string.( fieldnames( Simulation{T} ) ), ", " ) )
end

export Simulation
