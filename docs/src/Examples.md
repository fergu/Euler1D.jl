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

## Callbacks

Often it may be desirable to have the simulation stop at certain points in order to perform some processing, modify simulation state, or write the current simulation state to a file, amongst other things. To facilitate this, Euler1D supports callbacks that can be called based on various criteria. In particular, Euler1D supports callbacks that can be called at a regular cadence based on the number of cycles (timesteps) that have occurred ([`CycleCallback`](@ref)), at a regular cadence based on an elapsed time ([`TimeDeltaCallback`](@ref)), or at a fixed list of (potentially irregularly-spaced) times ([`TimeCallback`](@ref)). 

Setting up callbacks is straightforward. Considering the Sod Shock Tube example above, a structure containing information about callbacks can be created using [`ConfigureSimulationCallbacks()`](@ref):
```@repl sodsim
callbacks = ConfigureSimulationCallbacks(init_state)
```
From here, one can define functions that should be called at various defined points during the simulation. For illustrative purposes, we will define three callback functions here:
```@repl sodsim
# A callback that will be called at a fixed list of times
function MyTimeCallback( state::Simulation{T} ) where { T <: AbstractFloat }
    println("In MyTimeCallback: Current time is $(state.time.x)")
end;

# A callback that will be called at a regular number of cycles
function MyCycleCallback( state::Simulation{T} ) where { T <: AbstractFloat }
    println("In MyCycleCallback: Current cycle is $(state.cycles.x)")
end;

# A callback that will be called at a regular temporal cadence
function MyTimeDeltaCallback( state::Simulation{T} ) where { T <: AbstractFloat }
    println("In MyTimeDeltaCallback: Current time is $(state.time.x)")
end;
```

These three functions can then be assigned to the `callbacks` structure:
```@repl sodsim
RegisterTimeCallback!( callbacks, MyTimeCallback, Vector{Float64}([0.023, 0.067, 0.096]) ); # Register our time callback to run at t=0.023, 0.067, and 0.096
RegisterCycleCallback!( callbacks, MyCycleCallback, 250 ); # Register our cycle callback to run every 250 cycles
RegisterTimeDeltaCallback!( callbacks, MyTimeDeltaCallback, 0.03 ); # Register our time delta callback to run every 0.03 seconds

# We can also register the same callback multiple times.
RegisterTimeDeltaCallback!( callbacks, MyTimeDeltaCallback, 0.009, 0.05 ); # Let's register the time delta callback again, but starting at t=0.05 and with a increment of 0.009 seconds
RegisterCycleCallback!( callbacks, MyCycleCallback, 20, 1200 ); # We can also do the same thing with cycle-based callbacks. Here, we'll start at the 1200th cycle and print every 20 cycles
```

We can then run the simulation and pass our callbacks as an argument:
```@repl sodsim
end_state = AdvanceToTime( init_state, 0.1; exact=true, callbacks=callbacks )
```
As you can see, the output from each of the callbacks is visible. Here, a rather simple example of just printing to the screen was used, but other options such as saving to disk, making plots, or other manipulations are possible.

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

It is often useful to look at the evolution of the problem solution as a function of space and time rather than at just a single time instant. X-T (or space-time) diagrams are a useful way to perform this visualization. Gathering the data required to generate this type of plot is a little bit more complicated than the previous examples, however, so this section will walk through an example of how this might be done.

X-T diagrams can be made using quantities such as pressure or density, but [Riemann invariants](https://en.wikipedia.org/wiki/Riemann_invariant) are a particularly powerful way to visualize the movement of waves within the domain. The Euler equations have positive and negative Riemann invariants, corresponding to rightwards- and leftwards-moving waves, respectively:
```math
J_\pm = u \pm \frac{2c}{\gamma - 1}
```
where ``J_\pm`` is the positive and negative Riemann invariant, ``c`` is the speed of sound, and ``\gamma`` is the ratio of specific heats. The value of the Riemann invariant is a constant along a wave, and so visualizing isocontours of ``J_\pm`` as a function of space and time will correspond to the various waves moving around the domain.

As the above equation makes clear, we will need to combine multiple aspects of the simulation solution in order to compute the Riemann invariants. Additionally, the `contour` function from `Plots.jl` that we will use to plot the result assumes that the data it is provided lies on a regular grid, while `Euler1D` simulations are on a Lagrangian mesh, so we will need to interpolate our data to a regular grid. Fortunately, we can use callbacks to perform both of these tasks as the simulation runs.

!!! note
    The need to interpolate the data onto a regular grid is a specific limitation of the `contour` function in `Plots.jl`. Other plotting routines may exist that do not have this requirement, in which case the interpolation step will not be required.

To start, the simulation can be set up identically to the Sod shock tube used in previous examples: 
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
Now we need to construct our callback to compute the Riemann invariants and interpolate it onto a regular grid. For this, we'll use the grid at ``t=0`` as our reference grid:
```@repl xtsim
plot_positions = init_state.zone_center; # The x positions for plotting
```
We will also need some arrays to store the resulting Riemann invariant values into. For this example, we'll run the simulation to ``t=1.0``, while stopping to plot every 0.001 seconds:
```@repl xtsim
sim_end = 1.0; # Run the simulation to t=1.0 seconds
xt_dt = 0.001; # Stop to gather data for the X-T diagram every 0.001 seconds
N_callbacks = Int( div( sim_end, xt_dt ) ) + 1; # This is how many times the callback should be called. We add 1 to account for the first callback at t=0
xt_callback_count = 0; # How many times our plotting callback has been called. We'll use this to know where to save the Riemann invariant values for each callback
```
Using this, we can construct the arrays we'll use the save the Riemann invariants:
```@repl xtsim
plot_J₊ = zeros( N_callbacks, init_state.nzones ); # A matrix of the positive Riemann invariant values in every zone in the simulation
plot_J₋ = zeros( N_callbacks, init_state.nzones ); # A matrix of the negative Riemann invariant values at every zone in the simulation
plot_times = zeros( N_callbacks ); # The times at which the simulation data was saved. 
```
Now, finally, we can set up our actual callback routine
```@repl xtsim
using Interpolations

function xt_callback( state::Simulation{T} ) where { T <: AbstractFloat }
    global plot_times[begin + xt_callback_count] = state.time.x # Record our current time

    # Now compute the Riemann invariants
    # For this we'll need to interpolate the velocities to zone centers, which we'll do with a simple average
    uₘ = 0.5 .* ( state.velocity[1:end-1] .+ state.velocity[2:end] )

    # Now compute the Riemann invariants
    J₊ = uₘ .+ ( 2.0 .* state.speedofsound ) ./ ( state.gamma .- 1.0 ) # Positive Riemann invariant
    J₋ = uₘ .- ( 2.0 .* state.speedofsound ) ./ ( state.gamma .- 1.0 ) # Negative Riemann invariant

    # Most contour plot functions assume data is on a cartesian grid 
    # However, because the grid zones in the simulation are moving, we'll need to interpolate the data back onto a cartesian grid.
    # We can use Interpolations.jl for this

    # First the positive (right-moving) invariant
    J₊_interp = interpolate( ( state.zone_center, ), J₊, Gridded(Linear()) ) # Set up a linear interpolation of our data
    J₊_extrap = extrapolate( J₊_interp, Line() ) # Linearly extrapolate if needed.
    global plot_J₊[begin + xt_callback_count, :] = J₊_extrap( plot_positions ) # Interpolate the simulation data back onto the initial grid and add it to our matrix of data

    # Now the negative (left-moving) one
    J₋_interp = interpolate( ( state.zone_center, ), J₋, Gridded(Linear()) ) # Set up a linear interpolation of our data
    J₋_extrap = extrapolate( J₋_interp, Line() ) # Linearly extrapolate if needed
    global plot_J₋[begin + xt_callback_count, :] = J₋_extrap( plot_positions ) # Interpolate the simulation data back onto the initial grid and add it to our matrix of data

    # Finally, increment the counter that tracks the number of times this callback has been called
    global xt_callback_count += 1
end;
```
Finally, we can add this callback to our simulation:
```@repl xtsim
callbacks = ConfigureSimulationCallbacks(init_state);
RegisterTimeDeltaCallback!( callbacks, xt_callback, xt_dt ); # Register our X-T plotting callback
```
And now, we can run the simulation:
```@repl xtsim
end_state = AdvanceToTime( init_state, sim_end; exact=true, callbacks=callbacks )
```
With the simulation complete, we can now plot our X-T diagram. First, compute the minimum and maximum Riemann invariant values, as this will help us set the limits for our plot:
```@repl xtsim
minJ = min( minimum( plot_J₊ ), minimum( plot_J₋ ) );
maxJ = max( maximum( plot_J₊ ), maximum( plot_J₊ ) );
```
These values can then be used to create the contour plot
```@repl xtsim
using Plots
p = contour( plot_positions, plot_times, plot_J₊; c=:black, levels=minJ:0.2:maxJ, cbar=false ); # Positive invariants
contour!( p, plot_positions, plot_times, plot_J₋; c=:black, levels=minJ:0.2:maxJ ); # Negative invariants
savefig( p, "sod_xt.svg" );
```
These lines are responsible for actually creating the contour plot. The first line plots the positive (rightwards-moving) Riemann invariants, and the second line adds a plot for the negative (leftwards-moving) Riemann invariants. The optional arguments are `c`, which sets the line color, `levels`, which sets the contour levels to plot, and `cbar`, which tells the plotting routine not to plot a colorbar as it isn't useful in these contexts. Finally, `savefig` saves the plot to a file called `sod_xt.svg`. Running these lines gives us our X-T diagram:

![](sod_xt.svg)

And, voila! You can see the leftwards moving expansion wave and the rightwards moving shock and contact surface. Reflections off of the end walls are also visible, as are interactions between different waves. 

!!! note
    The `levels` argument, particularly the step size of `0.2`, was tuned to produce a good looking plot for this case. You will probably need to adjust this for other configurations.

!!! tip
    The plot colors were set to black as this makes a uniform looking plot. However, you might try plotting the two with different colors (say, `c=:red` in one plot) as this will show how the positive and negative Riemann invariants correspond to left- and right-moving waves.

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

Functions like these could, for example, be set up to run within a [`TimeCallback`](@ref) in order to update the problem state at a known time in order to, for example, add another shock wave.

!!! caution
    The functions to update the simulation state in this example were chosen to be illustrative, and it's likely the simulation would be unstable following this change. Care should be taken to ensure the updated state makes sense.
