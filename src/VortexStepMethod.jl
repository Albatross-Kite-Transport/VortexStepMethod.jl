module VortexStepMethod

using LinearAlgebra
using StaticArrays
using Logging
using Statistics
using Colors
using DelimitedFiles
using Measures
using LaTeXStrings
using NonlinearSolve
using Interpolations
using Interpolations: Extrapolation
using Parameters
using Serialization
using SharedArrays
using NonlinearSolve

# Export public interface
export Wing, Section, KiteWing
export BodyAerodynamics
export Solver, solve
export calculate_results, solve_circulation_distribution
export add_section!, set_va!
export calculate_span, calculate_projected_area
export menu
export Model, VSM, LLT
export AeroModel, LEI_AIRFOIL_BREUKELS, POLAR_DATA, INVISCID
export PanelDistribution, LINEAR, COSINE, COSINE_VAN_GARREL, SPLIT_PROVIDED, UNCHANGED
export InitialGammaDistribution, ELLIPTIC, ZEROS

export plot_geometry, plot_distribution, plot_circulation_distribution, plot_geometry, plot_polars, plot_foil_polars, save_plot, show_plot

# the following functions are defined in ext/VortexStepMethodExt.jl
function plot_geometry end
function plot_distribution end
function plot_circulation_distribution end
function plot_geometry end
function plot_polars end
function plot_foil_polars end
function save_plot end
function show_plot end

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
   AeroModel `LEI_AIRFOIL_BREUKELS` `POLAR_DATA` `INVISCID`

Enumeration of the implemented aerodynamic models.

# Elements
- `LEI_AIRFOIL_BREUKELS`: Polynom approximation for leading edge inflatable kites
- `POLAR_DATA`: Polar data (lookup tables with interpolation)
- INVISCID
"""
@enum AeroModel begin
   LEI_AIRFOIL_BREUKELS
   POLAR_DATA
   INVISCID
end

"""
   PanelDistribution `LINEAR` `COSINE` `COSINE_VAN_GARREL` `SPLIT_PROVIDED` `UNCHANGED`

Enumeration of the implemented panel distributions.

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

"""
   InitialGammaDistribution ELLIPTIC ZEROS

Enumeration of the implemented initial gamma distributions.

# Elements
- ELLIPTIC
- ZEROS
"""
@enum InitialGammaDistribution ELLIPTIC ZEROS

abstract type AbstractWing end

function menu()
   Main.include("examples/menu.jl")
end

# Include core functionality
include("wing_geometry.jl")
include("kite_geometry.jl")
include("filament.jl")
include("panel.jl")
include("body_aerodynamics.jl")
include("wake.jl")
include("solver.jl")

end # module