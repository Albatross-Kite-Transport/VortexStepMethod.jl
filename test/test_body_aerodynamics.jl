using VortexStepMethod
using VortexStepMethod: calculate_cl, calculate_cd_cm, calculate_projected_area, calculate_AIC_matrices!
using LinearAlgebra
using Test
using Logging

include("utils.jl")

@testset "Induction Matrix Creation" begin
    # Setup
    n_panels = 3
    N = n_panels + 1  # number of SECTIONS
    max_chord = 1.0
    span = 2.36
    AR = span^2 / (π * span * max_chord / 4)
    dist = "cos"
    coord = generate_coordinates_el_wing(max_chord, span, N, dist)

    Atot = max_chord / 2 * span / 2 * π
    @debug "N: $N"
    @debug "size(coord): $(size(coord))"

    v_a = 20.0
    aoa = 5.7106 * π / 180
    v_a = [cos(aoa), 0.0, sin(aoa)] .* v_a

    # Create wing geometry
    core_radius_fraction = 1e-20
    coord_left_to_right = flip_created_coord_in_pairs(deepcopy(coord))
    wing = Wing(n_panels; spanwise_panel_distribution=UNCHANGED)
    for idx in 1:2:size(coord_left_to_right, 1)
        add_section!(
            wing,
            coord_left_to_right[idx,:],
            coord_left_to_right[idx+1,:],
            INVISCID
        )
    end
    
    body_aero = BodyAerodynamics([wing])
    set_va!(body_aero, v_a)

    # Calculate reference matrices using thesis functions
    controlpoints, rings, bladepanels, ringvec, coord_L = 
        create_geometry_general(coord, v_a, N, "5fil", LLT)
    
    # Test LLT matrices
    @testset "LLT Matrices" begin
        # Calculate reference matrices
        MatrixU, MatrixV, MatrixW = thesis_induction_matrix_creation(
            deepcopy(ringvec),
            deepcopy(controlpoints),
            deepcopy(rings),
            deepcopy(v_a),
            zeros(N-1),
            nothing,  # data_airf not needed
            nothing,  # conv_crit not needed
            LLT
        )

        # Calculate new matrices
        va_norm_array = fill(norm(v_a), length(coord))
        va_unit_array = repeat(reshape(v_a ./ norm(v_a), 1, 3), length(coord))
        calculate_AIC_matrices!(
            body_aero,
            LLT,
            core_radius_fraction,
            va_norm_array,
            va_unit_array
        )
        AIC_x, AIC_y, AIC_z = @views body_aero.AIC[1, :, :], body_aero.AIC[2, :, :], body_aero.AIC[3, :, :]

        # Compare matrices
        @test isapprox(MatrixU, AIC_x, atol=1e-5)
        @test isapprox(MatrixV, -AIC_y, atol=1e-5)
        @test isapprox(MatrixW, AIC_z, atol=1e-5)
    end

    # Test VSM matrices
    @testset "VSM Matrices" begin
        # Calculate reference matrices for VSM
        controlpoints, rings, bladepanels, ringvec, coord_L = 
            create_geometry_general(coord, v_a, N, "5fil", VSM)
        
        MatrixU, MatrixV, MatrixW = thesis_induction_matrix_creation(
            deepcopy(ringvec),
            deepcopy(controlpoints),
            deepcopy(rings),
            deepcopy(v_a),
            zeros(N-1),
            nothing,
            nothing,
            VSM
        )

        # Calculate new matrices
        va_norm_array = fill(norm(v_a), length(coord))
        va_unit_array = repeat(reshape(v_a ./ norm(v_a), 1, 3), length(coord))
        calculate_AIC_matrices!(
            body_aero,
            VSM,
            core_radius_fraction,
            va_norm_array,
            va_unit_array
        )
        AIC_x, AIC_y, AIC_z = body_aero.AIC[1, :, :], body_aero.AIC[2, :, :], body_aero.AIC[3, :, :]

        # Compare matrices with higher precision for VSM
        @test isapprox(MatrixU, AIC_x, atol=1e-8)
        @test isapprox(MatrixV, -AIC_y, atol=1e-8)
        @test isapprox(MatrixW, AIC_z, atol=1e-8)
    end
end


@testset "Calculate results against output results" begin
    # Setup
    density = 1.225
    N = 40
    max_chord = 1.0
    span = 15.709  # AR = 20
    v_a = 20.0
    AR = span^2 / (π * span * max_chord / 4)
    aoa = deg2rad(5)
    v_a = [cos(aoa), 0.0, sin(aoa)] .* v_a
    model = VSM

    # Setup wing geometry
    dist = "cos"
    core_radius_fraction = 1e-20
    coord = generate_coordinates_el_wing(max_chord, span, N, dist)
    coord_left_to_right = flip_created_coord_in_pairs(deepcopy(coord))
    wing = Wing(N; spanwise_panel_distribution=UNCHANGED)
    for idx in 1:2:length(coord_left_to_right[:, 1])
        @debug "coord_left_to_right[$idx] = $(coord_left_to_right[idx,:])"
        add_section!(
            wing,
            coord_left_to_right[idx,:],
            coord_left_to_right[idx+1,:],
            INVISCID
        )
    end
    
    body_aero = BodyAerodynamics([wing])
    set_va!(body_aero, v_a)

    # Run analysis
    solver_object = Solver(
        aerodynamic_model_type=model, 
        core_radius_fraction=core_radius_fraction
    )
    results_NEW = solve(solver_object, body_aero; reference_point=[0,1,0])

    @test results_NEW isa Dict

    # Calculate forces using uncorrected alpha
    alpha = results_NEW["alpha_uncorrected"]
    dyn_visc = 0.5 * density * norm(v_a)^2
    n_panels = length(body_aero.panels)
    lift = zeros(n_panels)
    drag = zeros(n_panels)
    moment = zeros(n_panels)
    
    for (i, panel) in enumerate(body_aero.panels)
        lift[i] = dyn_visc * calculate_cl(panel, alpha[i]) * panel.chord
        cd_cm = calculate_cd_cm(panel, alpha[i])
        drag[i] = dyn_visc * cd_cm[1] * panel.chord
        moment[i] = dyn_visc * cd_cm[2] * panel.chord^2
        # @info "lift: $lift, drag: $drag, moment: $moment"
    end
    Fmag = hcat(lift, drag, moment)

    # Calculate coefficients using corrected alpha
    alpha = results_NEW["alpha_at_ac"]
    aero_coeffs = hcat(
        [alpha[i] for (i, panel) in enumerate(body_aero.panels)],
        [calculate_cl(panel, alpha[i]) for (i, panel) in enumerate(body_aero.panels)],
        [calculate_cd_cm(panel, alpha[i])[1] for (i, panel) in enumerate(body_aero.panels)],
        [calculate_cd_cm(panel, alpha[i])[2] for (i, panel) in enumerate(body_aero.panels)]
    )
    
    ringvec = [Dict("r0" => panel.width * panel.y_airf) for panel in body_aero.panels]
    controlpoints = [Dict("tangential" => panel.x_airf, "normal" => panel.z_airf) 
                    for panel in body_aero.panels]
    Atot = calculate_projected_area(wing)

    F_rel_ref, F_gl_ref, Ltot_ref, Dtot_ref, CL_ref, CD_ref, CS_ref = 
        output_results(Fmag, aero_coeffs, ringvec, v_a, controlpoints, Atot)

    # Compare results
    @info "Comparing results"
    @info "cl_calculated: $(results_NEW["cl"]), CL_ref: $CL_ref"
    @info "cd_calculated: $(results_NEW["cd"]), CD_ref: $CD_ref"
    @info "cs_calculated: $(results_NEW["cs"]), CS_ref: $CS_ref"
    @info "L_calculated: $(results_NEW["lift"]), Ltot_ref: $Ltot_ref"
    @info "D_calculated: $(results_NEW["drag"]), Dtot_ref: $Dtot_ref"

    # Assert results
    @test isapprox(results_NEW["cl"], CL_ref, rtol=1e-4)
    @test isapprox(results_NEW["cd"], CD_ref, rtol=1e-4)
    @test isapprox(results_NEW["cs"], CS_ref, rtol=1e-4)
    @test isapprox(results_NEW["lift"], Ltot_ref, rtol=1e-4)
    @test isapprox(results_NEW["drag"], Dtot_ref, rtol=1e-4)
    @test isapprox(results_NEW["Fx"], results_NEW["Mz"], rtol=1e-4) # 1 meter arm
    @test isapprox(results_NEW["My"], 0.0, atol=1e-3)
    @test isapprox(results_NEW["Fz"], -results_NEW["Mx"], rtol=1e-4) # 1 meter arm

    # Check array shapes
    @test length(results_NEW["cl_distribution"]) == length(body_aero.panels)
    @test length(results_NEW["cd_distribution"]) == length(body_aero.panels)
end

@testset "Wing Geometry Creation" begin
    function create_geometry(; model=VSM, wing_type=:rectangular, plotting=false, N=40)
        max_chord = 1.0
        span = 17.0
        AR = span^2 / (π * span * max_chord / 4)
        @debug "AR: $AR"
        v_a = 20.0
        aoa = 5.7106 * π / 180
        v_a = [cos(aoa), 0.0, sin(aoa)] .* v_a
    
        coord = if wing_type === :rectangular
            twist = range(-0.5, 0.5, length=N)
            beta = range(-2, 2, length=N)
            generate_coordinates_rect_wing(
                fill(max_chord, N),
                span,
                twist,
                beta,
                N,
                "lin"
            )
        elseif wing_type === :curved
            generate_coordinates_curved_wing(
                max_chord, span, π/4, 5, N, "cos"
            )
        elseif wing_type === :elliptical
            generate_coordinates_el_wing(max_chord, span, N, "cos")
        else
            error("Invalid wing type")
        end
    
        coord_left_to_right = flip_created_coord_in_pairs(deepcopy(coord))
        wing = Wing(N; spanwise_panel_distribution=UNCHANGED)
        for i in 1:2:size(coord_left_to_right, 1)
            add_section!(
                wing,
                coord_left_to_right[i,:],
                coord_left_to_right[i+1,:],
                INVISCID
            )
        end
        body_aero = BodyAerodynamics([wing])
        set_va!(body_aero, v_a)
        
        return body_aero, coord, v_a, model
    end

    for model in [VSM, LLT]
        @debug "model: $model"
        for wing_type in [:rectangular, :curved, :elliptical]
            @debug "wing_type: $wing_type"
            body_aero, coord, v_a, model = create_geometry(
                model=model, wing_type=wing_type
            )
            
            # Generate geometry
            expected_controlpoints, expected_rings, expected_bladepanels, 
                expected_ringvec, expected_coord_L = create_geometry_general(
                coord, v_a, div(size(coord,1), 2), "5fil", model
            )

            for i in 1:length(body_aero.panels)
                @debug "i: $i"
                # Handle control points
                index_reversed = length(body_aero.panels) - i + 1
                panel = body_aero.panels[index_reversed]
                
                evaluation_point = if model === VSM
                    panel.control_point
                else  # LLT
                    panel.aero_center
                end

                @test isapprox(evaluation_point, expected_controlpoints[i]["coordinates"], atol=1e-4)
                @test isapprox(panel.chord, expected_controlpoints[i]["chord"], atol=1e-4)
                @test isapprox(panel.z_airf, expected_controlpoints[i]["normal"], atol=1e-4)
                @test isapprox(panel.x_airf, expected_controlpoints[i]["tangential"], atol=1e-4)
                @test isapprox(
                    hcat(panel.z_airf, panel.x_airf, panel.y_airf),
                    expected_controlpoints[i]["airf_coord"],
                    atol=1e-4
                )
                
                if model === VSM
                    @test isapprox(
                        panel.aero_center,
                        expected_controlpoints[i]["coordinates_aoa"],
                        atol=1e-4
                    )
                end

                # Handle rings
                expected_ring_i = expected_rings[i]
                expected_ring_i_list = [
                    expected_ring_i[1],
                    expected_ring_i[2],
                    expected_ring_i[3],
                    expected_ring_i[4],
                    expected_ring_i[5]
                ]

                filaments = panel.filaments
                filament_list = [
                    filaments[1],
                    filaments[3],
                    filaments[5],
                    filaments[2],
                    filaments[4]
                ]

                for (j, fil) in enumerate(filament_list)
                    if j == 1  # bound filaments
                        @test isapprox(fil.x1, expected_ring_i_list[j]["x1"], atol=1e-4)
                        @test isapprox(fil.x2, expected_ring_i_list[j]["x2"], atol=1e-4)
                    elseif j ∈ (2, 4)  # trailing filaments
                        @test isapprox(fil.x1, expected_ring_i_list[j]["x1"], atol=1e-4)
                        @test isapprox(fil.x2, expected_ring_i_list[j]["x2"], atol=1e-4)
                    else  # semi-infinite filaments
                        @test isapprox(fil.x1, expected_ring_i_list[j]["x1"], atol=1e-4)
                    end
                end

                # Handle bladepanels
                exp_bladepanels = expected_bladepanels[i]
                @test isapprox(panel.LE_point_2, exp_bladepanels["p1"], atol=1e-4)
                @test isapprox(panel.LE_point_1, exp_bladepanels["p2"], atol=1e-4)
                @test isapprox(panel.TE_point_1, exp_bladepanels["p3"], atol=1e-4)
                @test isapprox(panel.TE_point_2, exp_bladepanels["p4"], atol=1e-4)

                # Handle ringvec
                exp_ringvec = expected_ringvec[i]
                r0 = panel.bound_point_1 - panel.bound_point_2
                r3 = evaluation_point - (panel.bound_point_1 + panel.bound_point_2) / 2
                @test isapprox(r0, exp_ringvec["r0"], atol=1e-4)
                @test isapprox(r3, exp_ringvec["r3"], atol=1e-4)

                # Handle coord_L
                @test all(isapprox.(panel.aero_center, expected_coord_L[:, i]))
            end
        end
    end
end
