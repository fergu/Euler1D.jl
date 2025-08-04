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
