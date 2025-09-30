using Test
using Euler1D


# Test the simulation solution of Sod's shock tube vs the exact solution
function TestSodShockTube()
    """
        SodShockTube( time::T, x::Vector{T}, γ::T, ρ₁::T, u₁::T, P₁::T, ρ₅::T, u₅::T, P₅::T ) where { T <: AbstractFloat }

    Compute the exact solution to Sod's shock tube with the specified initial conditions

    # Arguments
    - time: The time at which to compute the exact solution
    - x : A vector of points indicating the position of zone centers. Used for density and internal energy solutions.
    - xₑ: A vector of points indicating the position of zone edges. Used for the velocity solution.
    - xᵢ: The initial position of the interface between the two regions
    - γ : The ratio of specific heats of the gas
    - ρ : The initial density
    - u : The initial velocity
    - P : The initial pressure

    For density, velocity, and pressure, the subscripts refer to the regions of the solution where the values are defined. See Notes for a schematic of the solution.

    # Returns
    A `Dict{String,Vector{T}}` containing the following keys:
    - "x": The coordinate of each point in the solution
    - "density": The profile of density
    - "velocity": The profile of velocity
    - "intenergy": The profile of internal energy per unit mass

    # Notes
    Schematically, the solution states are
    (1) (Expansion:2) (3) (ContactSurface) (4) (Shock) (5)
    where the integer numbers indicates the index for that section used this function
    """
    function SodShockTubeExact( time::T, x::Vector{T}, xₑ::Vector{T}, xᵢ::T, γ::T, ρ₁::T, u₁::T, P₁::T, ρ₅::T, u₅::T, P₅::T ) where { T <: AbstractFloat }
        # This function references Chapter 7 of
        # Anderson, J, "Modern Compressible Flow with Historical Perspective", 3rd ed (2003)
        # for the equations governing each of the waves

        # Speed of sound
        c₁ = sqrt( γ * P₁ / ρ₁ ) # Left state
        c₅ = sqrt( γ * P₅ / ρ₅ ) # Right state
        
        # A few constants to make things easier
        Γ = ( γ - 1.0 ) / ( γ + 1.0 )
        β = ( γ - 1.0 ) / ( 2 * γ )

        # Now compute the states of the regions where the solution is constant
        # The expansion region (region 2) will need to be treated separately since it varies in space
        
        # State 4 (Between the shock and the contact surface)
        # We will use Eq. 7.94 to implicitly solve for the pressure jump across the shock in terms of the diaphragm pressure ratio
        ϵ = 100.0 # Starting error
        ϵₘ = 1e-6 # Maximum error (difference between pressure ratio estimates)
        # Initial guesses for the pressure ratio
        Pₗ = 1.0 # Initial lower guess is a pressure ratio of 1 (essentially an acoustic wave)
        Pₕ = 50.0 # Initial upper guess is a pressure ratio of 50, essentially a very strong shock
        P₄P₅ = 0.5 * ( Pₗ + Pₕ ) # An initial guess of an average between the upper and lower bounds
        while ( ϵ > ϵₘ )
            # Compute equation 7.94
            residual = P₄P₅ * ( 1.0 - ( ( γ - 1.0 ) * ( c₅ / c₁ ) * ( P₄P₅ - 1.0 ) ) / sqrt( 2.0 * γ * ( 2.0 * γ + ( γ + 1.0 ) * ( P₄P₅ - 1.0 ) ) ) )^( -2.0 * γ / ( γ - 1 ) ) - P₁/P₅
            # Check the residual.
            if ( residual >= 0.0 ) # Shock was too strong, lower our upper estimate
                Pₕ = P₄P₅
            else # Shock was too weak, raise our lower estimate
                Pₗ = P₄P₅
            end
            ϵ = Pₕ - Pₗ # Recompute the error as the difference between the upper and lower guesses
            P₄P₅ = 0.5 * ( Pₗ + Pₕ ) # Compute the next guess as the average of the upper and lower bounds
        end

        # Now we know the pressure jump across the shock. Can use this to solve for the state following the shock
        P₄ = P₄P₅ * P₅ # Use the solution to get the pressure after the shock wave
        u₄ = u₅ + c₅ / γ * ( P₄P₅ - 1.0 ) * sqrt( ( 2.0 * γ / ( γ + 1.0 ) ) / ( P₄P₅ + Γ ) ) # Eq. 7.16
        #ρ₄ = ρ₅ * ( P₄ + Γ * P₅ ) / ( P₅ + Γ * P₄ ) # Rankine-Hugioniot conditions give density after a contact surface
        ρ₄ = ρ₅ * ( ( 1.0 + 1.0/Γ * P₄P₅ ) / ( 1.0 / Γ + P₄P₅ ) ) # Eq 7.11, density ratio across the shock
        c₄ = sqrt( γ * P₄ / ρ₄ )

        # State 3 (Between contact surface and expansion wave)
        P₃ = P₄ # Pressure is constant across a contact surface
        u₃ = u₄ # Velocity is constant across a contact surface
        ρ₃ = ρ₁ * ( P₃ / P₁ )^( 1.0 / γ ) # Contact surface is isentropic, so density comes from isentropic relationships
        c₃ = sqrt( γ * P₃ / ρ₃ )

        # Lastly, compute the locations of the tip of the expansion wave, the contact surface, and the shock, as these will set where in physical space each of the above solutions is located
        # Shock velocity comes from the shock wave pressure ratio and Eq. 7.14 in Anderson
        xₛ = xᵢ + c₅ * sqrt( ( γ + 1.0 ) / ( 2.0 * γ ) * ( P₄P₅ - 1.0 ) + 1.0 ) * time
        # The contact surface has a velocity of u₃ = u₄, which we just solved for
        x₀ = xᵢ + u₃ * time
        # The left end (head) of the expansion fan can be found from the speed of sound of the gas it is propagating into
        xₑₗ = xᵢ - c₁ * time # A minus sign is added here to account for the direction of propagation
        # The right end of the expansion fan can be found using figure 7.15
        xₑᵣ = xᵢ + ( u₃ - c₃ ) * time

        # Now we can fill in the solution
        # First, fill in density and internal energy per unit mass
        ρ₀ = zeros( length( x ) )
        e₀ = zeros( length( x ) )
        for i in 1:length(x)
            thisX = x[i]
            if ( thisX < xₑₗ ) # To the left of the expansion fan (Region 1)
                ρ₀[i] = ρ₁
                e₀[i] = P₁ / ( ( γ - 1.0 ) * ρ₁ )
            elseif ( ( thisX >= xₑₗ ) && ( thisX < xₑᵣ ) ) # Inside of the expansion fan (Region 2)
                # We'll need Equations 7.86, 7.87, and 7.89 from Anderson for this
                u₂ = 2.0 / ( γ + 1.0 ) * ( c₁ + ( thisX - xᵢ ) / time ) # Equation 7.89
                P₂ = P₁ * ( 1.0 - ( γ - 1.0 ) / 2.0 * ( u₂ / c₁ ) )^( 2.0 * γ / ( γ - 1.0 ) ) # Equation 7.86
                ρ₂ = ρ₁ * ( 1.0 - ( γ - 1.0 ) / 2.0 * ( u₂ / c₁ ) )^( 2.0 / ( γ - 1.0 ) ) # Equation 7.87
                ρ₀[i] = ρ₂
                e₀[i] = P₂ / ( ( γ - 1.0 ) * ρ₂ )
            elseif ( ( thisX >= xₑᵣ ) && ( thisX < x₀ ) ) # Between the expansion fan and the contact surface (Region 3)
                ρ₀[i] = ρ₃
                e₀[i] = P₃ / ( ( γ - 1.0 ) * ρ₃ )
            elseif ( ( thisX >= x₀ ) && ( thisX < xₛ ) ) # Between the contact surface and the shock (Region 4)
                ρ₀[i] = ρ₄
                e₀[i] = P₄ / ( ( γ - 1.0 ) * ρ₄ )
            elseif ( thisX >= xₛ ) # To the right of the shock (Region 5)
                ρ₀[i] = ρ₅
                e₀[i] = P₅ / ( ( γ - 1.0 ) * ρ₅ )
            end
        end

        # Then fill in velocity.
        # This is done separately because velocity is edge-centered in the simulation and so we want to evaluate the solution at zone edges
        u₀ = zeros( length( xₑ ) )
        for i in 1:length(xₑ)
            thisX = xₑ[i]
            if ( thisX < xₑₗ ) # To the left of the expansion fan (Region 1)
                u₀[i] = u₁
            elseif ( ( thisX >= xₑₗ ) && ( thisX < xₑᵣ ) ) # Inside of the expansion fan (Region 2)
                # We'll need Equations 7.86, 7.87, and 7.89 from Anderson for this
                u₂ = 2.0 / ( γ + 1.0 ) * ( c₁ + ( thisX - xᵢ ) / time ) # Equation 7.89
                u₀[i] = u₂
            elseif ( ( thisX >= xₑᵣ ) && ( thisX < x₀ ) ) # Between the expansion fan and the contact surface (Region 3)
                u₀[i] = u₃
            elseif ( ( thisX >= x₀ ) && ( thisX < xₛ ) ) # Between the contact surface and the shock (Region 4)
                u₀[i] = u₄
            elseif ( thisX >= xₛ ) # To the right of the shock (Region 5)
                u₀[i] = u₅
            end
        end
        return Dict{String, Vector{T}}(
                                       "density"=>ρ₀,
                                       "velocity"=>u₀,
                                       "intenergy"=>e₀
                                      )
    end

    # Functions to set up the Sod's shock tube simulation
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

    final_time = 0.1

    # Run the simulation to the final state
    final_state = AdvanceToTime( init_state, final_time; exact=true ) # Run the simulation to the final time of 0.1 for comparison to Sod solution

    # Get the exact solution
    sod_exact = SodShockTubeExact( final_time, final_state.zone_center, final_state.zone_edge, 0.5, 1.4, 1.0, 0.0, 1.0, 0.125, 0.0, 0.1 )

    # Now compare the two
    # Errors will be checked using an averaged normalized L₂ error
    # The L₂ norm of a vector f is defined as
    L₂(f) = sqrt( sum( f.^2 ) )
    # The normalized error is then defined as
    #   ϵ₂ = || a - b ||₂ / || b ||₂
    # where a and b are the simulation and exact solutions, respectively
    L₂(a,b) = L₂( a .- b ) / L₂( b )
    # Finally, an average error is found by dividing the normalized L₂ error by the number of points
    L₂(a,b,N) = L₂(a,b) / N
    
    # Our testing threshold will be set as an average error of less than 0.1% from the exact solution
    ϵₘ = 1e-3
    # Check Density
    @test L₂( final_state.density, sod_exact["density"], final_state.nzones ) ≤ ϵₘ
    # Check Internal Energy
    @test L₂( final_state.intenergy, sod_exact["intenergy"], final_state.nzones ) ≤ ϵₘ
    # Check Velocity
    @test L₂( final_state.velocity, sod_exact["velocity"], final_state.nedges) ≤ ϵₘ
end

TestSodShockTube()
