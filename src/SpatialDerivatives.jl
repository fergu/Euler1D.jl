"""
    ∂∂x_ZoneCenterToZoneEdge( f₋::T, f₊::T, Δx₋::T, Δx₊::T ) where { T <: AbstractFloat }

Compute the derivative of a zone-centered variable, f, evaluated at the zone edge.

# Returns
A scalar of type `T` representing the value of ∂f/∂x evaluated at the boundary between two adjacent zones.

# Arguments
- `f`: A zone-centered variable located in the zone to the (-) left or (+) right of the zone boundary at which the derivative is evaluated.
- `Δx`: The length of the zone to the (-) left or (+) right of the zone boundary at which the derivative is evaluated.

# Notes
The derivative is calculated using a first-order forward difference. Functionally,

    ∂f/∂x[i+1/2] = ( f[i+1] - f[i] ) / ( xₘ[i+1] - xₘ[i] )

where f is the zone-centered variable and xₘ is the coordinate of the center of the zone
Since this code is Lagrangian, the location of zone edges is what is tracked, with the zone center computed as 

    xₘ[i] = 0.5 * ( x[i-1/2] + x[i+1/2] )

Therefore, the denominator above can be shown to be equivalent to

    xₘ[i+1] - xₘ[i] = 0.5 * ( x[i+3/2] + x[i+1/2] ) - 0.5 * ( x[i+1/2] + x[i-1/2] )

                    = 0.5 * ( x[i+3/2] - x[i+1/2] ) + 0.5 * ( x[i+1/2] - x[i-1/2] )

                    = 0.5 * Δx[i+1] + 0.5 * Δx[i]

                    = 0.5 * ( Δx[i+1] + Δx[i] )

and thus the derivative can be expressed as

    ∂f/∂x = 2.0 * ( f[i+1] - f[i] ) / ( Δx[i+1] + Δx[i] )
"""
function ∂∂x_ZoneCenterToZoneEdge( f₋::T, f₊::T, Δx₋::T, Δx₊::T ) where { T <: AbstractFloat }
    return 2.0 * ( f₊ - f₋ ) / ( Δx₊ + Δx₋ )
end

"""
    ∂∂x_ZoneCenterToZoneCenter( f₋::T, f::T, f₊::T, Δx₋::T, Δx::T, Δx₊::T ) where { T <: AbstractFloat }

Compute the derivative of a zone-centered variable, f, evaluated at the zone center.

# Returns
A scalar of type `T` representing the value of ∂f/∂x evaluated at the center of a zone.

# Arguments
- `f`: The value of the zone centered variable to the (-) left or (+) right of the (0) current zone.
- `Δx`: The length of the zone to the (-) left or (+) right of the (0) current zone.

# Notes
This derivative is calculated by first linearly interpolating the zone-centered variables to the left and right zone boundaries of the current zone. These interpolated values are then used to evaluate the derivative at the zone center.
The interpolated value at the left hand side of the zone is (where xₘ is the midpoint of the zone):

    f[i-1/2]    = f[i] + ∂f/∂x[i-1/2] * ( x[i-1/2] - xₘ[i] )

                = f[i] + ∂f/∂x[i-1/2] * ( x[i-1/2] - 0.5 * ( x[i+1/2] + x[i-1/2] ) )

                = f[i] + ∂f/∂x[i-1/2] * 0.5 * ( x[i-1/2] - x[i+1/2] )

                = f[i] - 0.5 * ∂f/∂x[i-1/2] * Δx[i]

And similarly for the right hand side of the zone:

    f[i+1/2]    = f[i] + ∂f/∂x[i+1/2] * ( x[i+1/2] - xₘ[i] )

                = f[i] + ∂f/∂x[i+1/2] * ( x[i+1/2] - 0.5 * ( x[i+1/2] + x[i-1/2] ) )

                = f[i] + ∂f/∂x[i-1/2] * 0.5 * ( x[i+1/2] - x[i-1/2] )

                = f[i] + 0.5 * ∂f/∂x[i-1/2] * Δx[i]

In both cases, `∂f/∂x[i+1/2]` (or `i-1/2`) can be evaluated using [`∂∂x_ZoneEdgeToZoneEdge()`](@ref).
The derivative is then:

    ∂f/∂x[i]    = ( f[i+1/2] - f[i-1/2] ) / ( x[i+1/2] - x[i-1/2] )

!!! warning
    This function is not currently used by this package and is mainly included for completeness. As such, it is not well-tested and should be used with caution.
"""
function ∂∂x_ZoneCenterToZoneCenter( f₋::T, f₀::T, f₊::T, Δx₋::T, Δx₀::T, Δx₊::T ) where { T <: AbstractFloat }
    ∂f∂x₋ = ∂f∂x_ZoneEdge( f₋, f₀, Δx₋, Δx₀ )
    f₋ = f₀ - 0.5 * ∂f∂x₋ * Δx
    ∂f∂x₊ = ∂f∂x_ZoneEdge( f₀, f₊, Δx₀, Δx₊ )
    f₊ = f₀ + 0.5 * ∂f∂x₊ * Δx
    return ( f₊ - f₋ ) / Δx
end

"""
    ∂∂x_ZoneEdgeToZoneCenter( f₋::T, f₊::T, Δx::T ) where { T <: AbstractFloat }

Evaluate the derivative of an edge-centered variable, f, at the center of the zone bounded by the two edges.

# Returns
A scalar of type `T` representing the value of ∂f/∂x evaluated at the center of a zone.

# Arguments
- `f`: The value of the edge-centered variable at the (-) left or (+) right edge of the zone.
- `Δx`: The length of the zone.

# Notes
This is computed using a simple forward difference with the values at the left and right edges of the zone.

    ∂f/∂x[i]    = ( f[i+1/2] - f[i-1/2] ) / ( x[i+1/2] - x[i-1/2] )

                = ( f[i+1/2] - f[i-1/2] ) / Δx[i]
"""
function ∂∂x_ZoneEdgeToZoneCenter( f₋::T, f₊::T, Δx::T ) where { T <: AbstractFloat }
    return ( f₊ - f₋ ) / Δx
end

"""
    ∂∂x_ZoneEdgeToZoneEdge( f₋::T, f::T, f₊::T, Δx₋::T, Δx₊::T ) where { T <: AbstractFloat }

Computes the derivative of an edge-centered variable, f, at the zone edge.

# Returns
A scalar of type `T` representing the value of ∂f/∂x

# Arguments
- `f`: The edge-centered variable to be evaluated. Subscripts refer to the zone edge to the (-) left or (+) right of the central (0) zone edge at which the derivative will be evaluated.
- `Δx`: The size of the zones to the (-) left or (+) right of the central zone edge.

# Notes
This function uses the values of the edge-centered variables to the (-) left and (+) right of the (0) central zone edge.  These left and right values are used to interpolate the edge centered values to the center of the zone to the (-,0) left or (0,+) right of the central zone edge.
These interpolation steps are performed according to:

    f[i+1]  = f[i+1/2] + ∂f/∂x[i+1] * ( xₘ[i+1] - x[i+1/2] )

            = f[i+1/2] + ∂f/∂x[i+1] * ( 0.5 * ( x[i+3/2] + x[i+1/2] ) - x[i+1/2] )

            = f[i+1/2] + ∂f/∂x[i+1] * 0.5 * ( x[i+3/2] - x[i+1/2] )

            = f[i+1/2] + 0.5 * ∂f/∂x[i+1] * Δx[i+1]

Similarly for the evaluation to the left:

    f[i]    = f[i+1/2] + ∂f/∂x[i] * ( xₘ[i] - x[i+1/2] )

            = f[i+1/2] + ∂f/∂x[i] * ( 0.5 * ( x[i+1/2] + x[i-1/2] ) - x[i+1/2] )

            = f[i+1/2] + ∂f/∂x[i] * 0.5 * ( x[i-1/2] - x[i+1/2] )

            = f[i+1/2] - 0.5 * ∂f/∂x[i] * Δx[i]

The calculation of `∂f/∂x` is performed using [`∂∂x_ZoneEdgeToZoneCenter()`](@ref).
The derivative evaluated at the central zone edge is therefore

    ∂f/dx[i+1/2]    = ( f[i+1] - f[i] ) / ( xₘ[i+1] - xₘ[i] )

                    = ( f[i+1] - f[i] ) / ( 0.5 * ( x[i+3/2] + x[i+1/2] ) - 0.5 * ( x[i+1/2] + x[i-1/2] ) )

                    = 2.0 * ( f[i+1] - f[i] ) / ( ( x[i+3/2] - x[i+1/2] ) + ( x[i+1/2] - x[i-1/2] ) )

                    = 2.0 * ( f[i+1] - f[i] ) / ( Δx[i+1] + Δx[i] )

!!! warning
    This function is not currently used by this package and is mainly included for completeness. As such, it is not well-tested and should be used with caution.
"""
function ∂∂x_ZoneEdgeToZoneEdge( f₋::T, f₀::T, f₊::T, Δx₋::T, Δx₊::T ) where { T <: AbstractFloat }
    ∂f∂x₋ = ∂∂x_ZoneEdgeToZoneCenter( f₋, f₀, Δx₋ ) # Derivative evaluated at the center of the zone to the left of the boundary
    f₋ = f₀ - 0.5 * ∂f∂x₋ * Δx₋ # The zone-edge value interpolated to the center of the zone to the left of the boundary
    ∂f∂x₊ = ∂∂x_ZoneEdgeToZoneCenter( f₀, f₊, Δx₊ ) # Derivative evaluated at the center of the zone to the right of the boundary
    f₊ = f₀ + 0.5 * ∂f∂x₊ * Δx₊ # The zone-edge value interpolated to the center of the zone to the right of the boundary
    return 2.0 * ( f₊ - f₋ ) / ( Δx₊ + Δx₋ )
end
