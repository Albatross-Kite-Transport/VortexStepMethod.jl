
"""
    Solver

Main solver structure for the Vortex Step Method.See also: [solve](@ref)

# Attributes

## General settings
- `aerodynamic_model_type`::Model = VSM: The model type, see: [Model](@ref)
- density::Float64 = 1.225: Air density [kg/m³] 
- `max_iterations`::Int64 = 1500
- `allowed_error`::Float64 = 1e-5: relative error
- `tol_reference_error`::Float64 = 0.001
- `relaxation_factor`::Float64 = 0.03: Relaxation factor for convergence 

## Damping settings
- `is_with_artificial_damping`::Bool = false: Whether to apply artificial damping
- `artificial_damping`::NamedTuple{(:k2, :k4), Tuple{Float64, Float64}} = (k2=0.1, k4=0.0): Artificial damping parameters

## Additional settings
- `type_initial_gamma_distribution`::InitialGammaDistribution = ELLIPTIC: see: [InitialGammaDistribution](@ref)
- `core_radius_fraction`::Float64 = 1e-20: 
- mu::Float64 = 1.81e-5: Dynamic viscosity [N·s/m²]
- `is_only_f_and_gamma_output`::Bool = false: Whether to only output f and gamma
"""
@with_kw struct Solver
    # General settings
    aerodynamic_model_type::Model = VSM
    density::Float64 = 1.225
    max_iterations::Int64 = 1500
    allowed_error::Float64 = 1e-5
    tol_reference_error::Float64 = 0.001
    relaxation_factor::Float64 = 0.03
    
    # Damping settings
    is_with_artificial_damping::Bool = false
    artificial_damping::NamedTuple{(:k2, :k4), Tuple{Float64, Float64}} =(k2=0.1, k4=0.0)
    
    # Additional settings
    type_initial_gamma_distribution::InitialGammaDistribution = ELLIPTIC
    core_radius_fraction::Float64 = 1e-20
    mu::Float64 = 1.81e-5
    is_only_f_and_gamma_output::Bool = false
end

"""
    solve(solver::Solver, body_aero::BodyAerodynamics, gamma_distribution=nothing; 
          log=false, reference_point=zeros(MVec3))

Main solving routine for the aerodynamic model. Reference point is in the kite body (KB) frame.

# Arguments:
- solver::Solver: The solver to use, could be a VSM or LLT solver. See: [Solver](@ref)
- body_aero::BodyAerodynamics: The aerodynamic body. See: [BodyAerodynamics](@ref)
- gamma_distribution: Initial circulation vector or nothing; Length: Number of segments. [m²/s]

# Keyword Arguments:
- log=false: If true, print the number of iterations and other info.
- reference_point=zeros(MVec3)

# Returns
A dictionary with the results.
"""
function solve(solver::Solver, body_aero::BodyAerodynamics, gamma_distribution=nothing; 
               log=false, reference_point=zeros(MVec3))
    
    # check arguments
    isnothing(body_aero.panels[1].va) && throw(ArgumentError("Inflow conditions are not set, use set_va!(body_aero, va)"))
    
    # Initialize variables
    panels = body_aero.panels
    n_panels = length(panels)
    alpha_array = body_aero.alpha_array
    relaxation_factor = solver.relaxation_factor
    
    # Preallocate arrays
    x_airf_array = zeros(n_panels, 3)
    y_airf_array = zeros(n_panels, 3)
    z_airf_array = zeros(n_panels, 3)
    va_array = zeros(n_panels, 3)
    chord_array = zeros(n_panels)

    # Fill arrays from panels
    for (i, panel) in enumerate(panels)
        x_airf_array[i, :] .= panel.x_airf
        y_airf_array[i, :] .= panel.y_airf
        z_airf_array[i, :] .= panel.z_airf
        va_array[i, :] .= panel.va
        chord_array[i] = panel.chord
    end

    # Calculate unit vectors
    va_norm_array = norm.(eachrow(va_array))
    va_unit_array = va_array ./ va_norm_array

    # Calculate AIC matrices
    calculate_AIC_matrices!(
        body_aero,
        solver.aerodynamic_model_type,
        solver.core_radius_fraction,
        va_norm_array,
        va_unit_array
    )

    # Initialize gamma distribution
    gamma_initial = if isnothing(gamma_distribution)
        if solver.type_initial_gamma_distribution === :elliptic
            calculate_circulation_distribution_elliptical_wing(body_aero)
        else
            zeros(n_panels)
        end
    else
        length(gamma_distribution) == n_panels || 
            throw(ArgumentError("gamma_distribution length must match number of panels"))
        gamma_distribution
    end

    @debug "Initial gamma_new: $gamma_initial"
    # Run main iteration loop
    converged, gamma_new, alpha_array, v_a_array = gamma_loop(
        solver,
        body_aero,
        gamma_initial,
        va_array,
        chord_array,
        x_airf_array,
        y_airf_array,
        z_airf_array,
        panels,
        relaxation_factor;
        log
    )
    # Try again with reduced relaxation factor if not converged
    if !converged && relaxation_factor > 1e-3
        log && @warn "Running again with half the relaxation_factor = $(relaxation_factor/2)"
        converged, gamma_new, alpha_array, v_a_array = gamma_loop(
            solver,
            body_aero,
            gamma_initial,
            va_array,
            chord_array,
            x_airf_array,
            y_airf_array,
            z_airf_array,
            panels,
            relaxation_factor/2;
            log
        )
    end

    # Calculate final results
    results = calculate_results(
        body_aero,
        gamma_new,
        reference_point,
        solver.density,
        solver.aerodynamic_model_type,
        solver.core_radius_fraction,
        solver.mu,
        alpha_array,
        v_a_array,
        chord_array,
        x_airf_array,
        y_airf_array,
        z_airf_array,
        va_array,
        va_norm_array,
        va_unit_array,
        panels,
        solver.is_only_f_and_gamma_output
    )
    return results
end

cross3(x,y) = cross(SVector{3,eltype(x)}(x), SVector{3,eltype(y)}(y))

"""
    gamma_loop(solver::Solver, gamma_new::Vector{Float64}, AIC_x::Matrix{Float64}, 
              AIC_y::Matrix{Float64}, AIC_z::Matrix{Float64}, va_array::Matrix{Float64}, 
              chord_array::Vector{Float64}, x_airf_array::Matrix{Float64}, 
              y_airf_array::Matrix{Float64}, z_airf_array::Matrix{Float64}, 
              panels::Vector{Panel}, relaxation_factor::Float64; log=true)

Main iteration loop for calculating circulation distribution.
"""
function gamma_loop(
    solver::Solver,
    body_aero::BodyAerodynamics,
    gamma_new::Vector{Float64},
    va_array::Matrix{Float64},
    chord_array::Vector{Float64},
    x_airf_array::Matrix{Float64},
    y_airf_array::Matrix{Float64},
    z_airf_array::Matrix{Float64},
    panels::Vector{Panel},
    relaxation_factor::Float64;
    log::Bool = true
)
    converged = false
    n_panels = length(body_aero.panels)
    alpha_array = body_aero.alpha_array
    v_a_array = body_aero.v_a_array
    Umagw_array = similar(v_a_array)

    gamma = copy(gamma_new)
    abs_gamma_new = copy(gamma_new)
    induced_velocity_all = zeros(n_panels, 3)
    relative_velocity_array = similar(va_array)
    relative_velocity_crossz = similar(relative_velocity_array)
    v_acrossz_array = similar(va_array)
    cl_array = zeros(n_panels)
    damp = zeros(length(gamma))
    v_normal_array = zeros(n_panels)
    v_tangential_array = zeros(n_panels)

    AIC_x, AIC_y, AIC_z = body_aero.AIC[1, :, :], body_aero.AIC[2, :, :], body_aero.AIC[3, :, :]

    velocity_view_x = @view induced_velocity_all[:, 1]
    velocity_view_y = @view induced_velocity_all[:, 2]
    velocity_view_z = @view induced_velocity_all[:, 3]

    iters = 0
    for i in 1:solver.max_iterations
        iters += 1
        gamma .= gamma_new
        
        # Calculate induced velocities
        mul!(velocity_view_x, AIC_x, gamma)
        mul!(velocity_view_y, AIC_y, gamma)
        mul!(velocity_view_z, AIC_z, gamma)
        
        relative_velocity_array .= va_array .+ induced_velocity_all
        for i in 1:n_panels
            relative_velocity_crossz[i, :] .= cross3(
                view(relative_velocity_array, i, :),
                view(y_airf_array, i, :)
            )
            v_acrossz_array[i, :] .= cross3(
                view(va_array, i, :),
                view(y_airf_array, i, :)
            )
        end

        for i in 1:n_panels
            v_normal_array[i] = dot(view(z_airf_array, i, :), view(relative_velocity_array, i, :))
            v_tangential_array[i] = dot(view(x_airf_array, i, :), view(relative_velocity_array, i, :))
        end
        alpha_array .= atan.(v_normal_array, v_tangential_array)

        for i in 1:n_panels
            @views v_a_array[i] = norm(relative_velocity_crossz[i, :])
            @views Umagw_array[i] = norm(v_acrossz_array[i, :])
        end
        
        for (i, (panel, alpha)) in enumerate(zip(panels, alpha_array))
            cl_array[i] = calculate_cl(panel, alpha)
        end
        gamma_new .= 0.5 .* v_a_array.^2 ./ Umagw_array .* cl_array .* chord_array

        # Apply damping if needed
        if solver.is_with_artificial_damping
            damp, is_damping_applied = smooth_circulation(gamma, 0.1, 0.5)
            @debug "damp: $damp"
        else
            damp .= 0.0
            is_damping_applied = false
        end
        # Update gamma with relaxation and damping
        gamma_new .= (1 - relaxation_factor) .* gamma .+ 
                    relaxation_factor .* gamma_new .+ damp

        # Check convergence
        abs_gamma_new .= abs.(gamma_new)
        reference_error = maximum(abs_gamma_new)
        reference_error = max(reference_error, solver.tol_reference_error)
        abs_gamma_new .= abs.(gamma_new .- gamma)
        error = maximum(abs_gamma_new)
        normalized_error = error / reference_error

        @debug "Iteration: $i, normalized_error: $normalized_error, is_damping_applied: $is_damping_applied"

        if normalized_error < solver.allowed_error
            converged = true
            break
        end
    end

    if log && converged
        @info "Converged after $iters iterations"
    elseif log
        @warn "NO convergence after $(solver.max_iterations) iterations"
    end

    return converged, gamma_new, alpha_array, v_a_array
end

"""
    smooth_circulation(circulation::Vector{Float64}, 
                      smoothness_factor::Float64, 
                      damping_factor::Float64)

Smooth circulation distribution if needed.

Returns:
- Tuple of smoothed circulation and boolean indicating if smoothing was applied
"""
function smooth_circulation(
    circulation::Vector{Float64},
    smoothness_factor::Float64,
    damping_factor::Float64
)
    # Calculate mean circulation excluding endpoints
    circulation_mean = mean(circulation[2:end-1])
    smoothness_threshold = smoothness_factor * circulation_mean

    # Calculate differences between adjacent points
    differences = diff(circulation[2:end-1])
    @debug "circulation_mean: $circulation_mean, diff: $differences"

    # Check smoothness
    if isempty(differences)
        return zeros(length(circulation)), false
    end

    if maximum(abs.(differences)) <= smoothness_threshold
        return zeros(length(circulation)), false
    end

    # Apply smoothing
    smoothed = copy(circulation)
    for i in 2:length(circulation)-1
        left = circulation[i-1]
        center = circulation[i]
        right = circulation[i+1]
        avg = (left + right) / 2
        smoothed[i] = center + damping_factor * (avg - center)
    end

    # Preserve total circulation
    total_original = sum(circulation)
    total_smoothed = sum(smoothed)
    smoothed .*= total_original / total_smoothed

    damp = smoothed - circulation
    return damp, true
end