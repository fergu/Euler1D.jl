# Function Reference

## Problem Configuration
```@docs
    DefaultSimulationParameters
    InitializeSimulation
    UpdateSimulationState!
```

## Timestepping
```@docs
    AdvanceToTime
    AdvanceOneCycle
    AdvanceOneCycle!
    AdvanceNCycles
    CalculateTimestepSize
```


## Equation of State
```@docs
    EOS_Density
    EOS_Pressure
    EOS_SpeedOfSound
```

## Artificial Dissipation
```@docs
    artificial_viscosity
    artificial_conductivity
```

## Types
```@docs
    Simulation
```

## Governing Equations
!!! note
    These functions are not intended to be called directly as part of simulation setup. However, as their functionality is central to this package, they are documented here for reference.
    
```@docs
    Euler1D.Momentum
    Euler1D.Energy
```

!!! warning
    These functions deal with calculating the rate of change of energy and momentum over the entire simulation domain. They should not be called as part of simulation configuration, but are documented here for completeness.

```@docs
    Euler1D.Momentum!
    Euler1D.Energy!
```

## Internal Functions

!!! warning
    These functions are internal to `Euler1D.jl` and are not intended for use in setting up simulations. Their function signatures may change at any time and without notice. The documentation of these functions is primarily included to support documentation of user-facing functionality.

```@docs
    Euler1D.∂∂x_ZoneCenterToZoneEdge
    Euler1D.∂∂x_ZoneCenterToZoneCenter
    Euler1D.∂∂x_ZoneEdgeToZoneCenter
    Euler1D.∂∂x_ZoneEdgeToZoneEdge
```
