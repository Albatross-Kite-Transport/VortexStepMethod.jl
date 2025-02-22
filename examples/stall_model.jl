using VortexStepMethod
using LinearAlgebra

using Pkg
if ! ("CSV" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using CSV
using DataFrames

# Find root directory
root_dir = dirname(@__DIR__)
save_folder = joinpath(root_dir, "results", "TUDELFT_V3_LEI_KITE")
mkpath(save_folder)

# Defining discretisation
n_panels = 54
spanwise_panel_distribution = "split_provided"

# Load rib data from CSV
csv_file_path = joinpath(
    root_dir,
    "processed_data",
    "TUDELFT_V3_LEI_KITE",
    "rib_list_from_CAD_LE_TE_and_surfplan_d_tube_camber_19ribs.csv"
)

df = CSV.read(csv_file_path, DataFrame)
rib_list = []
for row in eachrow(df)
    LE = [row.LE_x, row.LE_y, row.LE_z]
    TE = [row.TE_x, row.TE_y, row.TE_z]
    push!(rib_list, (LE, TE, ("lei_airfoil_breukels", [row.d_tube, row.camber])))
end

# Create wing geometry
CAD_wing = Wing(n_panels; spanwise_panel_distribution)
for rib in rib_list
    add_section!(CAD_wing, rib[1], rib[2], rib[3])
end
wing_aero_CAD_19ribs = WingAerodynamics([CAD_wing])

# Create solvers
VSM = Solver(
    aerodynamic_model_type="VSM",
    is_with_artificial_damping=false
)
VSM_with_stall_correction = Solver(
    aerodynamic_model_type="VSM",
    is_with_artificial_damping=true
)

# Setting velocity conditions
v_a = 15.0
aoa = 17.0
side_slip = 0.0
yaw_rate = 0.0
aoa_rad = deg2rad(aoa)
vel_app = [
    cos(aoa_rad) * cos(side_slip),
    sin(side_slip),
    sin(aoa_rad)
] * v_a
wing_aero_CAD_19ribs.va = vel_app

# # Plotting geometry
# plot_geometry(
#     wing_aero_CAD_19ribs,
#     "";
#     data_type=".svg",
#     save_path="",
#     is_save=false,
#     is_show=true,
#     view_elevation=15,
#     view_azimuth=-120
# )

# Solving and plotting distributions
results = solve(VSM, wing_aero_CAD_19ribs)
@time results_with_stall = solve(VSM_with_stall_correction, wing_aero_CAD_19ribs)
@time results_with_stall = solve(VSM_with_stall_correction, wing_aero_CAD_19ribs)

CAD_y_coordinates = [panel.aerodynamic_center[2] for panel in wing_aero_CAD_19ribs.panels]

plot_distribution(
    [CAD_y_coordinates, CAD_y_coordinates],
    [results, results_with_stall],
    ["VSM", "VSM with stall correction"];
    title="CAD_spanwise_distributions_alpha_$(round(aoa, digits=1))_beta_$(round(side_slip, digits=1))_yaw_$(round(yaw_rate, digits=1))_v_a_$(round(v_a, digits=1))",
    data_type=".pdf",
    save_path=joinpath(save_folder, "spanwise_distributions"),
    is_save=false,
    is_show=true
)

# Plotting polar
save_path = joinpath(root_dir, "results", "TUD_V3_LEI_KITE")
path_cfd_lebesque = joinpath(
    root_dir,
    "data",
    "TUDELFT_V3_LEI_KITE",
    "literature_results",
    "V3_CL_CD_RANS_Lebesque_2024_Rey_300e4.csv"
)

plot_polars(
    [VSM, VSM_with_stall_correction],
    [wing_aero_CAD_19ribs, wing_aero_CAD_19ribs],
    [
        "VSM CAD 19ribs",
        "VSM CAD 19ribs , with stall correction",
        "CFD_Lebesque Rey 30e5"
    ];
    literature_path_list=[path_cfd_lebesque],
    angle_range=range(0, 25, length=25),
    angle_type="angle_of_attack",
    angle_of_attack=0,
    side_slip=0,
    v_a=10,
    title="tutorial_testing_stall_model_n_panels_$(n_panels)_distribution_$(spanwise_panel_distribution)",
    data_type=".pdf",
    save_path=joinpath(save_folder, "polars"),
    is_save=true,
    is_show=true
)
nothing