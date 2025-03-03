```@meta
CurrentModule = VortexStepMethod
```
## Enumerations
```@docs
Model
AeroModel
PanelDistribution
InitialGammaDistribution
SolverStatus
```

## Basic Vectors
```@docs
MVec3
PosVector
VelVector
```

## Aerodynamic data
```@docs
AeroData
```

## Wing Geometry, Panel and Aerodynamics
A body is constructed of one or more abstract wings. An abstract wing can be a Wing or a KiteWing. 
A Wing/ KiteWing has one or more sections.
```@docs
Section
Section(LE_point::Vector{Float64}, TE_point::Vector{Float64}, aero_model=nothing, aero_data=nothing)
Wing
Wing(n_panels::Int; spanwise_panel_distribution::PanelDistribution=LINEAR,
     spanwise_direction::PosVector=MVec3([0.0, 1.0, 0.0]))
KiteWing
KiteWing(obj_path, dat_path; alpha=0.0, crease_frac=0.75, wind_vel=10., mass=1.0, 
         n_panels=54, n_sections=n_panels+1, spanwise_panel_distribution=UNCHANGED, 
         spanwise_direction=[0.0, 1.0, 0.0])
BodyAerodynamics
```

## The Solver and its results
```@docs
Solver
Result
```
