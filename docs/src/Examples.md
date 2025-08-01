# Examples

## Sod Shock Tube

The [Sod shock tube](https://en.wikipedia.org/wiki/Sod_shock_tube) is a case involving a rightwards travelling shock wave and a leftwards travelling expansion wave. Its initial condition is given by:

```math
\left(
\begin{aligned}
\rho_L \\
P_L \\
u_L
\end{aligned}
\right) =
\left(
\begin{aligned}
1.0 \\
1.0 \\
0.0
\end{aligned}
\right),~
\left(
\begin{aligned}
\rho_R \\
P_R \\
u_R
\end{aligned}
\right) =
\left(
\begin{aligned}
0.125 \\
0.1 \\
0.0
\end{aligned}
\right)
```

where ``\rho``, ``P``, and ``u`` are the density, pressure, and velocity, the domain is in the range ``x=[0,1]``, and the subscripts `L` and `R` refer to states to the left and right of an interface initially located at ``x=0.5``. While an exact solution for this configuration can be found, this example will focus on setting up an `Euler1D` simulation to simulate the Sod shock tube.

A simulation is initialized using a set of dictionary key:value pairs describing various aspects of the simulation, including the functions used to describe the initial profiles of density, velocity, pressure, and the ratio of specific heats. A default set of parameters can be obtained using the function [`DefaultSimulationParameters()`][@ref]:

```@repl sodsim
using Euler1D
simulation_parameters = DefaultSimulationParameters()
```
You may notice that there are several parameters with a value of `nothing`. These must be supplied before the simulation can be initialized since it is difficult to make reasonable default assumptions. Of particular note are the keys `init_density_function`, `init_velocity_function`, `init_pressure_function`, and `init_gamma_function`. These refer to functions that describe the inital profiles of density, velocity, pressure, and the ratio of specific heats, respectively. A set of functions to implement the Sod shock tube initial conditions looks like:
```@repl sodsim
function init_density( x )
    if ( x < 0.5 )
        return 1.0
    else
        return 0.125
    end
end;

function init_velocity( x )
    return 0.0
end;

function init_pressure( x )
    if ( x < 0.5 )
        return 1.0
    else
        return 0.1
    end
end;

function init_gamma( x )
    return 1.4
end;
```
!!! caution
    The `x` argument in these functions is the location of where a given variable will be stored. Due to the fact that some variables are zone-centered while others are edge-centered, it should not be assumed that `x` will have the same value in all functions. Additionally, due to how Julia performs vectorized computations, these functions should not assume anything about the order in which they are executed. i.e., that the `n`-th function call will always have a certain value of `x`.

The `simulation_parameters` dictionary can now be modified to describe the desired initial condition. In this case, this would look like:
```@repl sodsim
simulation_parameters["init_density_function"] = init_density;
simulation_parameters["init_velocity_function"] = init_velocity;
simulation_parameters["init_pressure_function"] = init_pressure;
simulation_parameters["init_gamma_function"] = init_gamma;
```
The `simulation_parameters` dictionary has many other keys that can be modified. See [`DefaultSimulationParameters()`](@ref) for a list of parameters that can be set.

Once all desired parameters have been set, the simulation can be initialized:
```@repl sodsim
init_state = InitializeSimulation( simulation_parameters )
println( init_state ) # hide
```
The `init_state` variable is a [`Simulation`](@ref) type that holds information about the simulation. 
!!! note
    The functions that describe the initial conditions are not called until [`InitializeSimulation()`](@ref) is called, at which time they are called at every point in the domain. For these reasons, it is good to be sure that the functions are not computationally heavy.

A common time to examine the Sod shock tube solution is `t=0.1`. The [`AdvanceToTime()`](@ref) function can be used to advance the simulation to this time:
```@repl sodsim
end_state = AdvanceToTime( init_state, 0.1; exact=true )
```
where `0.1` is the time to advance to and the `exact=true` keyword tells the simulation to modify the final time step to stop as close as possible to the final time. If `exact=false` were supplied instead (or if the `exact` argument was omitted), the simulation would stop as soon as the current time is greater than the specified stopping time, but no modification of the timestep would be made and so the actual stopping time may differ from the specified stopping time. The size of the difference will depend on the timestep size.

!!! tip
    This example wrote out the functions used for initial conditions in full for instructive purposes. However, there is no requirement that the functions be defined in this way. For example, for simple initial conditions, a function could be defined using [ternary operators](https://docs.julialang.org/en/v1/manual/control-flow/):
    ```@repl
    init_density( x ) = x < 0.5 ? 1.0 : 0.125; # A function that varies in space
    init_gamma( x ) = 1.4; # A function with a constant value
    ```
    or even using [anonymous functions](https://docs.julialang.org/en/v1/manual/functions/#man-anonymous-functions):
    ```@repl
    simulation_parameters = Dict{String,Any}() # hide
    simulation_parameters["init_density_function"] = (x) -> x < 0.5 ? 1.0 : 0.125; # A function that varies in space
    simulation_parameters["init_gamma_function"] = (x) -> 1.4; # A function with a constant value
    ```
    There is no requirement or advantage to using one method over any other. These are simply alternative ways of defining a function to describe the initial condition.

## Plotting Results

This section outlines a few examples of how simulation data can be post-processed to visualize results or perform advanced processing.

!!! note
    There are numerous packages in the Julia ecosystem to handle plotting (e.g. [Plots.jl](https://docs.juliaplots.org/stable/), [Makie.jl](https://docs.makie.org/stable/)), writing to disk (e.g. [CSV.jl](https://csv.juliadata.org/stable/), [DataFrames.jl](https://dataframes.juliadata.org/stable/), [HDF5.jl](https://juliaio.github.io/HDF5.jl/stable/)), and any other data analysis routine that might be desired. The people developing these packages (and others like them) are immensely talented, and one would be well-advised to utilize their hard work if a package with suitable functionality exists. 
    For this reason, this package actively avoids making assumptions about what a user will want to do with the simulation data after the simulation is complete. This means that there is no built-in functionality to plot simulation data, write it to disk, or perform any other sort of postprocessing. Instead, all simulation outputs are simply elementary Julia types (e.g. `Vector{Float64}`), and can be manipulated using the usual built-in Julia methods, or by using any suitable package that provides the required functionality. 

### Line Plot

The simplest case of post-processing data is to just plot profiles of the solution state at the final time in the simulation. Fortunately, all necessary fields can be obtained from the `Simulation` variable in the form of elementary Julia types. For example, to plot the profile of density from the simulation performed in [Sod Shock Tube](@ref) using the [Plots.jl](https://docs.juliaplots.org/stable/) library, one can do
```@repl sodsim
using Plots
plot( end_state.zone_center, end_state.density );
savefig( "sod_density.svg" ) # hide
```
![](sod_density.svg)

Similar plots can be made for other zone-centered quantities such as internal energy or pressure. 

!!! tip
    If plotting edge-centered quantities such as velocity, the `end_state.zone_edge` variable should be used for the spatial coordinate.

### X-T diagrams

It is often useful to look at the evolution of the problem solution as a function of space and time rather than at just a single time instant. X-T diagrams are a useful way to perform this visualization. Gathering the data required to generate this type of plot is a little bit more complicated than the previous examples, however.

The simulation can be set up identically to previous examples: 
```@repl xtsim
using Euler1D # hide
simulation_parameters = DefaultSimulationParameters();
# We'll define our initialization functions using anonymous functions to keep things short
simulation_parameters["init_density_function"] = (x) -> x < 0.5 ? 1.0 : 0.125;
simulation_parameters["init_velocity_function"] = (x) -> 0.0;
simulation_parameters["init_pressure_function"] = (x) -> x < 0.5 ? 1.0 : 0.1;
simulation_parameters["init_gamma_function"] = (x) -> 1.4;
init_state = InitializeSimulation( simulation_parameters )
```
In order to generate an X-T diagram, we will need to record the simulation state at multiple times in the simulation. To do this, let's define a new function that will be responsible for advancing the simulation through time: 
```@repl xtsim
using Interpolations
using Plots

function run_simulation( init_state::Simulation{T}, end_time::T ) where { T <: AbstractFloat }
    new_state = deepcopy( init_state ) # We'll need a copy of our initial state to modify inside the while loop
    # Information about plotting
    plot_dt = 0.001 # How often to stop and record data
    next_plot = 0.0 # The next time to stop and record data
    # Set up some matricies to hold our plotting data
    plot_times = zeros( 0 ) # The times at which we've recorded data
    plot_positions = init_state.zone_center # Where to plot the data at in space. We'll use the initial grid
    plot_J₊ = zeros( 0, init_state.nzones ) # A matrix of the positive Riemann invariant values in every zone in the simulation
    plot_J₋ = zeros( 0, init_state.nzones ) # A matrix of the negative Riemann invariant values at every zone in the simulation
    while ( new_state.time.x <= end_time ) # Keep looping until we've reached our end time
        new_state = AdvanceToTime( new_state, next_plot; exact=true ) # Advance the solution to the next time to record data
        # Now record data about the simulation state
        plot_times = vcat( plot_times, new_state.time.x ) # Record our current time, appending it to the plot_times vector
        # Now compute the Riemann invariants. 
        # For this we'll need to interpolate the velocities to zone centers, which we'll do with a simple average
        uₘ = 0.5 .* ( new_state.velocity[1:end-1] .+ new_state.velocity[2:end] )
        # Now compute the Riemann invariants
        J₊ = uₘ .+ ( 2.0 .* new_state.speedofsound ) ./ ( new_state.gamma .- 1.0 ) # Positive Riemann invariant
        J₋ = uₘ .- ( 2.0 .* new_state.speedofsound ) ./ ( new_state.gamma .- 1.0 ) # Negative Riemann invariant
        # Most contour plot functions assume data is on a cartesian grid 
        # However, because the grid zones in the simulation are moving, we'll need to interpolate the data back onto a cartesian grid.
        # For simplicity, we will use the initial simulation grid as the grid for plotting
        # We can use Interpolations.jl for this
        J₊_interp = interpolate( ( new_state.zone_center, ), J₊, Gridded(Linear()) ) # Set up a linear interpolation of our data
        J₊_extrap = extrapolate( J₊_interp, Line() ) # Linearly extrapolate if needed.
        plot_J₊ = vcat( plot_J₊, J₊_extrap( plot_positions )' ) # Interpolate the simulation data back onto the initial grid and append it to our matrix of data

        J₋_interp = interpolate( ( new_state.zone_center, ), J₋, Gridded(Linear()) ) # Set up a linear interpolation of our data
        J₋_extrap = extrapolate( J₋_interp, Line() ) # Linearly extrapolate if needed
        plot_J₋ = vcat( plot_J₋, J₋_extrap( plot_positions )' ) # Interpolate the simulation data back onto the initial grid and append it to our matrix of data
        next_plot = next_plot + plot_dt
    end
    # Now that the simulation is done, we can save our plot
    # First compute the minimum and maximum Riemann invariant values. This helps set the contour levels for our plot
    minJ = min( minimum( plot_J₊ ), minimum( plot_J₋ ) ) 
    maxJ = max( maximum( plot_J₊ ), maximum( plot_J₊ ) )
    # Now plot the Riemann invariants as contours
    p = contour( plot_positions, plot_times, plot_J₊; c=:black, levels=minJ:0.2:maxJ, cbar=false ) # Positive invariants
    contour!( p, plot_positions, plot_times, plot_J₋; c=:black, levels=minJ:0.2:maxJ ) # Negative invariants
    savefig( p, "sod_xt.svg" )
    return new_state
end;
```
There's a lot going on in this function, but fortunately most of it is straightforward and is well explained by the comments. However, there's a few things worth highlighting. In particular, X-T diagrams can be made using quantities such as pressure or density, but [Riemann Invariants](https://en.wikipedia.org/wiki/Riemann_invariant) are a particularly powerful way to visualize the movement of waves within the domain. The Euler equations have positive and negative Riemann invariants, corresponding to rightwards- and leftwards-moving waves, respectively. These are computed using the following lines:
```julia
# Now compute the Riemann invariants. 
# For this we'll need to interpolate the velocities to zone centers, which we'll do with a simple average
uₘ = 0.5 .* ( new_state.velocity[1:end-1] .+ new_state.velocity[2:end] )
# Now compute the Riemann invariants
J₊ = uₘ .+ ( 2.0 .* new_state.speedofsound ) ./ ( new_state.gamma .- 1.0 ) # Positive Riemann invariant
J₋ = uₘ .- ( 2.0 .* new_state.speedofsound ) ./ ( new_state.gamma .- 1.0 ) # Negative Riemann invariant
```
These lines first translate the velocity (which is edge-centered) to the zone centers, and then uses that average velocity along with the speed of sound and ratio of specific heats to compute the positive and negative Riemann invariants as a function of space.

Now that we have the Riemann invariants, we need to store them for plotting. However, most contour plotting functions (which we will use to plot our X-T diagram) assume that the data is stored on a rectangular grid. In this context, this would mean that the `x` positions are the same for all time. However, because this is a Lagrangian code, the `x` positions of the data are not constant as a function of time, so we need to massage the data into a format the contour plotting function can accept. This is done by interpolating the data at a given time back onto the initial grid using the `Interpolations.jl` package and the following lines: 
```julia
J₊_interp = interpolate( ( new_state.zone_center, ), J₊, Gridded(Linear()) ) # Set up a linear interpolation of our data
J₊_extrap = extrapolate( J₊_interp, Line() ) # Linearly extrapolate if needed.
plot_J₊ = vcat( plot_J₊, J₊_extrap( plot_positions )' ) # Interpolate the simulation data back onto the initial grid and append it to our matrix of data
        
J₋_interp = interpolate( ( new_state.zone_center, ), J₋, Gridded(Linear()) ) # Set up a linear interpolation of our data
J₋_extrap = extrapolate( J₋_interp, Line() ) # Linearly extrapolate if needed
plot_J₋ = vcat( plot_J₋, J₋_extrap( plot_positions )' ) # Interpolate the simulation data back onto the initial grid and append it to our matrix of data
```
The `interpolate` and `extrapolate` functions set up this interpolation, and then the `vcat` lines actually perform the interpolation and append the result to the matricies we're using to store our data.

Now that the simulation is done, we'll configure the contour plot to actually produce the X-T diagram. First, we want to compute the minimum and maximum values of the Riemann invariants. We do this so that the contour plots for the positive and negative invariants use the same contour levels, which helps make sure lines are connected to each other and aids in visualization.
```julia
# First compute the minimum and maximum Riemann invariant values. This helps set the contour levels for our plot
minJ = min( minimum( plot_J₊ ), minimum( plot_J₋ ) ) 
maxJ = max( maximum( plot_J₊ ), maximum( plot_J₊ ) )
```
We can now actually create the contour plot:
```julia
# Now plot the Riemann invariants as contours
p = contour( plot_positions, plot_times, plot_J₊; c=:black, levels=minJ:0.2:maxJ, cbar=false ) # Positive invariants
contour!( p, plot_positions, plot_times, plot_J₋; c=:black, levels=minJ:0.2:maxJ ) # Negative invariants
savefig( p, "sod_xt.svg" )
```
These lines are responsible for actually creating the contour plot. The first line plots the positive (rightwards-moving) Riemann invariants, and the second line adds a plot for the negative (leftwards-moving) Riemann invariants. The optional arguments are `c`, which sets the line color, `levels`, which sets the contour levels to plot, and `cbar`, which tells the plotting routine not to plot a colorbar as it isn't useful in these contexts. Finally, `savefig` saves the plot to a file called `sod_xt.svg`.

!!! tip
    The plot colors were set to black as this makes a uniform looking plot. However, you might try plotting the two with different colors (say, `c=:red` in one plot) as this will show how the positive and negative Riemann invariants correspond to left- and right-moving waves.

!!! note
    The `levels` argument, particularly the step size of `0.2`, was tuned to produce a good looking plot for this case. You will probably need to adjust this for other configurations.
    
With everything set up we can run our simulation and take a look at our plot. For illustration purposes, the final time in the simulation will be a later time in order to better show the movement of waves in the final plot.
```@repl xtsim
final_state = run_simulation( init_state, 1.0 )
```
![](sod_xt.svg)

and there we go! You can see the leftwards moving expansion wave and the rightwards moving shock and contact surface. Reflections off of the end walls are also visible, as are interactions between different waves. 

!!! tip
    You might notice that things like the material interface look a little wider than might be expected. Try playing around with the `artificial_viscosity_coefficient` and `artificial_conductivity_coefficient` values to see what effect these parameters have.

## Modifying Problem State

There may be cases where it is desirable to stop and modify the problem state part way through a simulation. For example, to add another shock wave. This can be done using the [`UpdateSimulationState!()`](@ref) function. This function takes four arguments that are `Function`s describing the new problem density, velocity, pressure, and ratio of specific heats. However, in this case the functions have a slightly different signature:
```julia
function my_new_state( x, oldValue )
    ...
end
```
Notice that there is now an `oldValue` argument. This will hold the current value of the variable at the position `x`. If you don't want to modify anything, you can just return `oldValue` from this function.

As a trivial example, let's say we want to modify the example from the [Line Plot](@ref) section to add a sinusoidal profile to the existing density, set the velocity to zero, and leave everything else untouched. Our functions might look like
```@repl sodsim
function update_density( x, oldValue )
    return 0.05 * sin( 20 * pi * x ) + 0.05 + oldValue # Need to add a value to make sure density doesn't go to zero
end;

function update_velocity( x, oldValue )
    return 0.0
end;

function update_pressure( x, oldValue )
    return oldValue
end;

function update_gamma( x, oldValue )
    return oldValue
end;
```

With these functions we can now update the problem state:
```@repl sodsim
UpdateSimulationState!( end_state, update_gamma, update_density, update_velocity, update_pressure )
```
and if we plot density again, we'll see the density field has been updated:
```@repl sodsim
plot( end_state.zone_center, end_state.density );
savefig( "sod_density_updated.svg" ) # hide
```
![](sod_density_updated.svg)

and the velocity field is now zero everywhere:
```@repl sodsim
plot( end_state.zone_edge, end_state.velocity );
savefig( "sod_velocity_updated.svg" ) # hide
```
![](sod_velocity_updated.svg)

!!! caution
    The functions to update the simulation state were chosen to be illustrative, and it's likely the simulation would be unstable following this change. Care should be taken to ensure the updated state makes sense.
