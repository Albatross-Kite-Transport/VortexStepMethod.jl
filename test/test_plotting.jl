using VortexStepMethod
using ControlPlots
using Test

if !@isdefined ram_wing
    body_path = joinpath(tempdir(), "ram_air_kite_body.obj")
    foil_path = joinpath(tempdir(), "ram_air_kite_foil.dat")
    cp("data/ram_air_kite_body.obj", body_path; force=true)
    cp("data/ram_air_kite_foil.dat", foil_path; force=true)
    ram_wing = RamAirWing(body_path, foil_path; alpha_range=deg2rad.(-1:1), delta_range=deg2rad.(-1:1))
end

function create_body_aero()
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
    body_aero = BodyAerodynamics([wing])
    # Set inflow conditions
    vel_app = [cos(alpha), 0.0, sin(alpha)] .* v_a
    set_va!(body_aero, vel_app)
    body_aero
end

plt.ioff()
@testset "Plotting" begin
    fig = plt.figure(figsize=(14, 14))
    res = plt.plot([1,2,3])
    @test fig isa plt.PyPlot.Figure
    @test res isa Vector{plt.PyObject}
    save_plot(fig, "/tmp", "plot")
    @test isfile("/tmp/plot.pdf")
    rm("/tmp/plot.pdf")
    show_plot(fig)
    body_aero = create_body_aero()
    if Sys.islinux()
        fig = plot_geometry(
            body_aero,
            "Rectangular_wing_geometry";
            data_type=".pdf",
            save_path="/tmp",
            is_save=true,
            is_show=false)
        @test fig isa plt.PyPlot.Figure
        @test isfile("/tmp/Rectangular_wing_geometry_angled_view.pdf")
        rm("/tmp/Rectangular_wing_geometry_angled_view.pdf")
        @test isfile("/tmp/Rectangular_wing_geometry_front_view.pdf")
        rm("/tmp/Rectangular_wing_geometry_front_view.pdf")
        @test isfile("/tmp/Rectangular_wing_geometry_side_view.pdf")
        rm("/tmp/Rectangular_wing_geometry_side_view.pdf")
        @test isfile("/tmp/Rectangular_wing_geometry_top_view.pdf")
        rm("/tmp/Rectangular_wing_geometry_top_view.pdf")

        # Step 5: Initialize the solvers
        vsm_solver = Solver(body_aero; aerodynamic_model_type=VSM)
        llt_solver = Solver(body_aero; aerodynamic_model_type=LLT)

        # Step 6: Solve the VSM and LLT
        results_vsm = solve(vsm_solver, body_aero)
        results_llt = solve(llt_solver, body_aero)

        # Step 7: Plot spanwise distributions
        y_coordinates = [panel.aero_center[2] for panel in body_aero.panels]

        fig = plot_distribution(
            [y_coordinates, y_coordinates],
            [results_vsm, results_llt],
            ["VSM", "LLT"],
            title="Spanwise Distributions"
        )
        @test fig isa plt.PyPlot.Figure

        # Step 8: Plot polar curves
        v_a = 20.0            # Magnitude of inflow velocity [m/s]
        angle_range = range(0, 20, 20)
        fig = plot_polars(
            [llt_solver, vsm_solver],
            [body_aero, body_aero],
            ["VSM", "LLT"],
            angle_range=angle_range,
            angle_type="angle_of_attack",
            v_a=v_a,
            title="Rectangular Wing Polars",
            data_type=".pdf",
            save_path="/tmp",
            is_save=true,
            is_show=false
        )
        @test fig isa plt.PyPlot.Figure
        @test isfile("/tmp/Rectangular_Wing_Polars.pdf")
        rm("/tmp/Rectangular_Wing_Polars.pdf")

        # Step 9: Test polar data plotting
        body_aero = BodyAerodynamics([ram_wing])
        fig = plot_polar_data(body_aero; is_show=false)
        @test fig isa plt.PyPlot.Figure
    end
end
plt.ion()
nothing