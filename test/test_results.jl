using VortexStepMethod
using VortexStepMethod: calculate_cl, calculate_cd_cm, calculate_projected_area, calculate_AIC_matrices!, init!
using LinearAlgebra
using Test
using Logging

if !@isdefined ram_wing
    body_path = joinpath(tempdir(), "ram_air_kite_body.obj")
    foil_path = joinpath(tempdir(), "ram_air_kite_foil.dat")
    cp("data/ram_air_kite_body.obj", body_path; force=true)
    cp("data/ram_air_kite_foil.dat", foil_path; force=true)
    ram_wing = RamAirWing(body_path, foil_path; alpha_range=deg2rad.(-1:1), delta_range=deg2rad.(-1:1))
end

@testset "Nonlinear vs Linear - Comprehensive Input Testing" begin
    # Initialize with base parameters
    va = [15.0, 0.0, 0.0]
    theta = zeros(4)
    delta = zeros(4)
    omega = zeros(3)
    
    # Define perturbation magnitudes
    dva_magnitudes = [0.01, 0.01, 0.01]  # Velocity perturbations (m/s)
    dtheta_magnitudes = [deg2rad(0.1), deg2rad(0.5), deg2rad(1.0), deg2rad(2.0)] # Angle perturbations (rad)
    ddelta_magnitudes = [deg2rad(0.1), deg2rad(0.5), deg2rad(1.0), deg2rad(2.0)] # Trailing edge perturbations (rad)
    domega_magnitudes = [deg2rad(0.1), deg2rad(0.5), deg2rad(1.0)]  # Angular rate perturbations (rad/s)
    
    # Create body aerodynamics and solver
    VortexStepMethod.group_deform!(ram_wing, theta, delta; smooth=false)
    body_aero = BodyAerodynamics([ram_wing]; va, omega)
    solver = Solver(body_aero;
        aerodynamic_model_type=VSM,
        is_with_artificial_damping=false,
        atol=1e-8,
        rtol=1e-8,
        solver_type=NONLIN
    )
    
    # Get initial Jacobian and linearization result
    base_inputs = [theta; va; omega; delta]
    jac, lin_res = VortexStepMethod.linearize(
        solver, 
        body_aero, 
        base_inputs; 
        theta_idxs=1:4, 
        va_idxs=5:7, 
        omega_idxs=8:10,
        delta_idxs=11:14,
        moment_frac=0.1
    )
    
    # Verify that linearization results match nonlinear results at operating point
    baseline_res = VortexStepMethod.solve!(solver, body_aero; log=false)
    baseline_res = [solver.sol.force; solver.sol.moment; solver.sol.group_moment_dist]
    @test baseline_res ≈ lin_res
    
    # Define test cases
    test_cases = [
        ("va", 5:7, dva_magnitudes),
        ("theta", 1:4, dtheta_magnitudes),
        ("omega", 8:10, domega_magnitudes),
        ("delta", 11:14, ddelta_magnitudes)
    ]
    
    # For each test case
    for (input_name, indices, magnitudes) in test_cases
        @testset "Input: $input_name" begin
            max_error_ratio = 0.0
            
            # Select one index at a time for focused testing
            for idx in indices
                # Test each magnitude
                for mag in magnitudes
                    # Prepare perturbation vector
                    perturbation = zeros(length(base_inputs))
                    
                    if length(indices) > 1 && input_name != "theta" && input_name != "delta"
                        # For vector quantities (va, omega), perturb all components together
                        perturbation[indices] .= mag
                    else
                        # For control angles, perturb one at a time
                        perturbation[idx] = mag
                    end
                    
                    # Reset to baseline
                    reset_va = copy(va)
                    reset_theta = copy(theta)
                    reset_delta = copy(delta)
                    reset_omega = zeros(3)
                    
                    # Apply the perturbation to the nonlinear model
                    if input_name == "va"
                        reset_va = va + perturbation[5:7]
                    elseif input_name == "theta"
                        reset_theta = theta + perturbation[1:4]
                    elseif input_name == "delta"
                        reset_delta = delta + perturbation[11:14]
                    elseif input_name == "omega"
                        reset_omega = perturbation[8:10]
                    else
                        throw(ArgumentError())
                    end
                    VortexStepMethod.group_deform!(ram_wing, reset_theta, reset_delta; smooth=false)
                    init!(body_aero; init_aero=false, va=reset_va, omega=reset_omega)
                    
                    # Get nonlinear solution
                    nonlin_res = VortexStepMethod.solve!(solver, body_aero, nothing; log=false)
                    nonlin_res = [solver.sol.force; solver.sol.moment; solver.sol.group_moment_dist]
                    @test nonlin_res ≉ baseline_res
                    
                    # Compute linearized prediction
                    lin_prediction = lin_res + jac * perturbation
                    
                    # Calculate error ratio for this case
                    prediction_error = norm(lin_prediction - nonlin_res)
                    baseline_difference = norm(baseline_res - nonlin_res)
                    error_ratio = prediction_error / baseline_difference

                    # Update max error ratio
                    max_error_ratio = max(max_error_ratio, error_ratio)
                    
                    # For small perturbations, test that error ratio is small
                    if idx == first(indices) && mag == first(magnitudes)
                        @test error_ratio < 2e-3
                    end
                end
            end
            
            @info "Max error ratio for $input_name: $max_error_ratio"
        end
    end
    
    # Test combinations of input variations
    @testset "Combined Input Variations" begin
        # For combination testing, we'll create targeted test cases
        # that use only the specific indices we want to test together
        combination_tests = [
            # Name, active indices, perturbation vector, indices mappings
            (
                "Theta + VA", 
                # Use only theta and va indices
                [1:4; 5:7],  
                # Perturbation values for these indices
                [dtheta_magnitudes; dva_magnitudes],
                # Mappings for linearize function
                (theta_idxs=1:4, va_idxs=5:7, omega_idxs=nothing, delta_idxs=nothing)
            ),
            (
                "Theta + Omega", 
                [1:4; 8:10],
                [dtheta_magnitudes; domega_magnitudes],
                (theta_idxs=1:4, va_idxs=nothing, omega_idxs=5:7, delta_idxs=nothing)
            ),
            (
                "VA + Omega", 
                [5:7; 8:10],
                [dva_magnitudes; domega_magnitudes],
                (theta_idxs=nothing, va_idxs=1:3, omega_idxs=4:6, delta_idxs=nothing)
            ),
            (
                "Delta + VA", 
                [5:7; 11:14],
                [dva_magnitudes; ddelta_magnitudes],
                (theta_idxs=nothing, va_idxs=1:3, omega_idxs=nothing, delta_idxs=4:7)
            ),
            (
                "All Inputs", 
                [1:4; 5:7; 8:10; 11:14],
                [dtheta_magnitudes; dva_magnitudes; 
                domega_magnitudes; ddelta_magnitudes],
                (theta_idxs=1:4, va_idxs=5:7, omega_idxs=8:10, delta_idxs=11:14)
            )
        ]
        
        for (combo_name, active_indices, perturbation, idx_mappings) in combination_tests
            @testset "$combo_name" begin
                # Start with a fresh model for each combination test
                VortexStepMethod.group_deform!(ram_wing, theta, delta; smooth=false)
                init!(body_aero; init_aero=false, va, omega)
                
                # Create the appropriate input vector for this combination
                input_vec = Vector{Float64}(undef, length(active_indices))
                
                # Fill with base values
                if !isnothing(idx_mappings.theta_idxs)
                    input_vec[idx_mappings.theta_idxs] .= theta
                end
                if !isnothing(idx_mappings.va_idxs)
                    input_vec[idx_mappings.va_idxs] .= va
                end
                if !isnothing(idx_mappings.omega_idxs)
                    input_vec[idx_mappings.omega_idxs] .= omega
                end
                if !isnothing(idx_mappings.delta_idxs)
                    input_vec[idx_mappings.delta_idxs] .= delta
                end
                
                # Get the Jacobian and linearization result for this specific combination
                jac_combo, lin_res_combo = VortexStepMethod.linearize(
                    solver, 
                    body_aero, 
                    input_vec; 
                    idx_mappings...
                )
                
                # Get baseline results
                baseline_res = VortexStepMethod.solve!(solver, body_aero; log=false)
                baseline_res = [solver.sol.force; solver.sol.moment; solver.sol.group_moment_dist]
                
                # Should match the linearization result
                @test baseline_res ≈ lin_res_combo
                
                # Apply perturbation using the appropriate indices
                perturbed_input = copy(input_vec) + perturbation
                
                # Extract components based on the combination being tested
                perturbed_theta = !isnothing(idx_mappings.theta_idxs) ? 
                    perturbed_input[idx_mappings.theta_idxs] : theta
                
                perturbed_va = !isnothing(idx_mappings.va_idxs) ? 
                    perturbed_input[idx_mappings.va_idxs] : va
                
                perturbed_omega = !isnothing(idx_mappings.omega_idxs) ? 
                    perturbed_input[idx_mappings.omega_idxs] : omega
                
                perturbed_delta = !isnothing(idx_mappings.delta_idxs) ? 
                    perturbed_input[idx_mappings.delta_idxs] : delta
                
                # Apply to nonlinear model
                VortexStepMethod.group_deform!(ram_wing, perturbed_theta, perturbed_delta; smooth=false)
                init!(body_aero; init_aero=false, va=perturbed_va, omega=perturbed_omega)
                
                # Get nonlinear solution with perturbation
                nonlin_res = VortexStepMethod.solve!(solver, body_aero; log=false)
                nonlin_res = [solver.sol.force; solver.sol.moment; solver.sol.group_moment_dist]
                
                # Compute linearized prediction using our specialized Jacobian
                lin_prediction = lin_res_combo + jac_combo * perturbation
                
                # Ensure perturbation had an effect
                @test baseline_res ≉ nonlin_res atol=1e-3
                
                # Calculate error ratio
                prediction_error = norm(lin_prediction - nonlin_res)
                baseline_difference = norm(baseline_res - nonlin_res)
                error_ratio = prediction_error / baseline_difference
                
                @info "$combo_name error metrics" prediction_error baseline_difference error_ratio
                
                # Validate the prediction
                @test lin_prediction ≈ nonlin_res rtol=0.1 atol=1e-3
                @test error_ratio < 0.05
            end
        end
    end
end