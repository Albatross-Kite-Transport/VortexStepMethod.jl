using LinearAlgebra
using ControlPlots
using VortexStepMethod

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
wing = Wing(n_panels, spanwise_panel_distribution=LINEAR)

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
llt_solver = Solver(aerodynamic_model_type=LLT)
vsm_solver = Solver(aerodynamic_model_type=VSM)

# Step 5: Solve using both methods
results_vsm = solve(vsm_solver, wa)
results_vsm_base = solve_base(vsm_solver, wa)
println("Rectangular wing:")
@time results_vsm_base = solve_base(vsm_solver, wa)
# time Python: 32.0 ms Ryzen 7950x
# time Julia:   0.6 ms Ryzen 7950x
#               0.8 ms laptop, performance mode, battery 

# Create wing geometry
wing = KiteWing("data/ram_air_kite_body.obj", "data/ram_air_kite_foil.dat")
body_aero = BodyAerodynamics([wing])

# Create solvers
vsm_solver = Solver(
    aerodynamic_model_type=VSM,
    is_with_artificial_damping=false
)

# Setting velocity conditions
v_a = 15.0
aoa = 15.0
side_slip = 0.0
yaw_rate = 0.0
aoa_rad = deg2rad(aoa)
vel_app = [
    cos(aoa_rad) * cos(side_slip),
    sin(side_slip),
    sin(aoa_rad)
] * v_a
set_va!(body_aero, vel_app)

# Solving and plotting distributions
results = solve(vsm_solver, body_aero)
results_base = solve_base(vsm_solver, body_aero)
println("RAM-air kite:")
@time results_base = solve_base(vsm_solver, body_aero)

nothing