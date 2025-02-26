module VortexStepMethod

using LinearAlgebra
using StaticArrays
using Logging
using Statistics
using Colors
using DelimitedFiles
using ControlPlots
using Measures
using LaTeXStrings
using NonlinearSolve
using Interpolations
using Interpolations: Extrapolation
using Serialization
using SharedArrays

# Export public interface
export Wing, Section, KiteWing
export BodyAerodynamics
export Solver, solve
export calculate_results, solve_circulation_distribution
export add_section!, set_va!
export calculate_span, calculate_projected_area
export plot_wing, plot_circulation_distribution, plot_geometry, plot_distribution, plot_polars
export show_plot, save_plot, menu
export Model, VSM, LLT
export PanelDistribution, LINEAR, COSINE, COSINE_VAN_GARREL, SPLIT_PROVIDED, UNCHANGED

"""
   const MVec3    = MVector{3, Float64}

Basic 3-dimensional vector, stack allocated, mutable.
"""
const MVec3    = MVector{3, Float64}

"""
   const PosVector=Union{MVec3, Vector, SizedVector{3, Float64, Vector{Float64}}}

Position vector, either a `MVec3` or a `Vector` for use in function signatures.
"""
const PosVector=Union{MVec3, Vector}

"""
   const VelVector=Union{MVec3, Vector, SizedVector{3, Float64, Vector{Float64}}}

Velocity vector, either a `MVec3` or a `Vector` for use in function signatures.
"""
const VelVector=Union{MVec3, Vector}

"""
   Model `VSM` `LLT`

Enumeration of the implemented model types.

# Elements
- VSM: Vortex Step Method
- LLT: Lifting Line Theory
"""
@enum Model VSM LLT

"""
   PanelDistribution `LINEAR` `COSINE` `COSINE_VAN_GARREL` `SPLIT_PROVIDED` `UNCHANGED`

Enumeration of the implemented model types.

# Elements
- LINEAR               # Linear distribution
- COSINE               # Cosine distribution
- `COSINE_VAN_GARREL`  # van Garrel cosine distribution
- `SPLIT_PROVIDED`     # Split provided sections
- UNCHANGED            # Keep original sections
"""
@enum PanelDistribution begin
   LINEAR             # Linear distribution
   COSINE             # Cosine distribution
   COSINE_VAN_GARREL  # van Garrel cosine distribution
   SPLIT_PROVIDED     # Split provided sections
   UNCHANGED          # Keep original sections
end

abstract type AbstractWing end

function menu()
   Main.include("examples/menu.jl")
end

# Include core functionality
include("wing_geometry.jl")
include("kite_geometry.jl")
include("filament.jl")
include("panel.jl")
include("wake.jl")
include("body_aerodynamics.jl")
include("solver.jl")

# include plotting
include("plotting.jl")

end # module