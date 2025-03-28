# DO NOT DELETE, DO NOT MODIFY!
# This script will not be included in run_tests.jl
# It shall only be executed on one and the same machine for benchmarking the solve functions
using LinearAlgebra
using VortexStepMethod
using BenchmarkTools
using Test

using Pkg

if !("CSV" ∈ keys(Pkg.project().dependencies))
    using TestEnv
    TestEnv.activate()
end

# Step 1: Define wing parameters
n_panels = 20          # Number of panels
span = 20.0            # Wing span [m]
chord = 1.0            # Chord length [m]
v_a = 20.0             # Magnitude of inflow velocity [m/s]
density = 1.225        # Air density [kg/m³]
alpha_deg = 30.0       # Angle of attack [degrees]
alpha = deg2rad(alpha_deg)

# Step 2: Create wing geometry with linear panel distribution
wing = Wing(n_panels, spanwise_distribution=LINEAR)

# Add wing sections - defining only tip sections with inviscid airfoil model
add_section!(wing, 
    [0.0, span/2, 0.0],    # Left tip LE 
    [chord, span/2, 0.0],  # Left tip TE
    INVISCID)
add_section!(wing, 
    [0.0, -span/2, 0.0],   # Right tip LE
    [chord, -span/2, 0.0], # Right tip TE
    INVISCID)

# Step 3: Initialize aerodynamics
wa = BodyAerodynamics([wing])

# Set inflow conditions
vel_app = [cos(alpha), 0.0, sin(alpha)] .* v_a
set_va!(wa, vel_app)

# Step 4: Initialize solvers for both LLT and VSM methods
vsm_solver = Solver(wa; aerodynamic_model_type=VSM)

# Step 5: Solve using both methods
result = @benchmark  solve_base!($vsm_solver, $wa, nothing)  # 34 allocations
println("solve_base():")
println("Allocations: ", result.allocs, ", Mean time: ", round(mean(result.times)/1000), " µs")
# time Python: 32.0 ms  Ryzen 7950x
# time Julia:   0.45 ms Ryzen 7950x
result = @benchmark  sol = solve!($vsm_solver, $wa, nothing) # 68 allocations
println("solve!()")
println("Allocations: ", result.allocs, ", Mean time: ", round(mean(result.times)/1000), " µs")
nothing
