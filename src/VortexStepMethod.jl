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

# Export public interface
export Wing, Section, KiteWing
export BodyAerodynamics
export Solver, solve, solve_base
export calculate_results, solve_circulation_distribution
export add_section!, set_va!
export calculate_span, calculate_projected_area
export menu, MVec3
export Model, VSM, LLT
export AeroModel, LEI_AIRFOIL_BREUKELS, POLAR_VECTORS, POLAR_MATRICES, INVISCID
export PanelDistribution, LINEAR, COSINE, COSINE_VAN_GARREL, SPLIT_PROVIDED, UNCHANGED
export InitialGammaDistribution, ELLIPTIC, ZEROS

export plot_geometry, plot_distribution, plot_circulation_distribution, plot_geometry, plot_polars, save_plot, show_plot

# the following functions are defined in ext/VortexStepMethodExt.jl
function plot_geometry end
function plot_distribution end
function plot_circulation_distribution end
function plot_geometry end
function plot_polars end
function save_plot end
function show_plot end

"""
   const MVec3    = MVector{3, Float64}

Basic 3-dimensional vector, stack allocated, mutable.
"""
const MVec3    = MVector{3, Float64}

"""
   const PosVector=Union{MVec3, Vector}

Position vector, either a `MVec3` or a `Vector` for use in function signatures.
"""
const PosVector=Union{MVec3, Vector}

"""
   const VelVector=Union{MVec3, Vector}

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
   AeroModel `LEI_AIRFOIL_BREUKELS` `POLAR_VECTORS` `POLAR_MATRICES` `INVISCID`

Enumeration of the implemented aerodynamic models. See also: [AeroData](@ref)

# Elements
- `LEI_AIRFOIL_BREUKELS`: Polynom approximation for leading edge inflatable kites
- `POLAR_VECTORS`: Polar vectors as function of alpha (lookup tables with interpolation)
- `POLAR_MATRICES`: Polar matrices as function of alpha and beta (lookup tables with interpolation)
- INVISCID

where `alpha` is the angle of attack, `beta` is trailing edge angle.
"""
@enum AeroModel begin
   LEI_AIRFOIL_BREUKELS
   POLAR_VECTORS
   POLAR_MATRICES
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

"""
   SolverStatus FEASIBLE INFEASIBLE FAILURE

Enumeration to report back the validity of the result of the solve! function.
Used in the [Result](@ref) struct.

# Elements
- FEASIBLE: The gamma distribution is physically feasible
- INFEASIBLE: The gamma distribution is physically infeasible
- FAILURE: The result did not converge within the maximal number of iterations
"""
@enum SolverStatus FEASIBLE INFEASIBLE FAILURE

abstract type AbstractWing end

"""
    AeroData= Union{
        Nothing,
        NTuple{2, Float64},
        Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}},
        Tuple{Vector{Float64}, Vector{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}
    }

Union of different definitions of the aerodynamic properties of a wing section. See also: [AeroModel](@ref)
  - nothing for INVISCID
  - (`tube_diameter`, camber) for `LEI_AIRFOIL_BREUKELS`
  - (`alpha_range`, `cl_vector`, `cd_vector`, `cm_vector`) for `POLAR_VECTORS`
  - (`alpha_range`, `beta_range`, `cl_matrix`, `cd_matrix`, `cm_matrix`) for `POLAR_MATRICES` 

where `alpha` is the angle of attack [rad], `beta` is trailing edge angle [rad], `cl` the lift coefficient,
`cd` the drag coefficient and `cm` the pitching moment coefficient. The camber of a kite refers to 
the curvature of its airfoil shape. The camber is typically measured as the maximum distance 
between the mean camber line (the line equidistant from the upper and lower surfaces) 
and the chord line of the airfoil.
"""
const AeroData = Union{
        Nothing,
        NTuple{2, Float64},
        Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}},
        Tuple{Vector{Float64}, Vector{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}
    }

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