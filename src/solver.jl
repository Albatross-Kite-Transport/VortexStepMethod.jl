
"""
    Solver

Main solver structure for the Vortex Step Method.

# Attributes

## General settings
- `aerodynamic_model_type`::Model: The model type, see: [Model](@ref)
- density::Float64: Air density [kg/m³]
- `max_iterations`::Int64
- `reltol`::Float64: relative error
- `tol_reference_error`::Float64
- `relaxation_factor`::Float64

## Damping settings
- `is_with_artificial_damping`::Bool
- `artificial_damping`::NamedTuple{(:k2, :k4), Tuple{Float64, Float64}}

## Additional settings
- `type_initial_gamma_distribution`::InitialGammaDistribution: see: [InitialGammaDistribution](@ref)
- `core_radius_fraction`::Float64
- mu::Float64
- `is_only_f_and_gamma_output`::Bool
"""
struct Solver
    # General settings
    aerodynamic_model_type::Model
    density::Float64
    max_iterations::Int64
    reltol::Float64
    abstol::Float64
    tol_reference_error::Float64
    relaxation_factor::Float64
    
    # Damping settings
    is_with_artificial_damping::Bool
    artificial_damping::NamedTuple{(:k2, :k4), Tuple{Float64, Float64}}
    
    # Additional settings
    type_initial_gamma_distribution::InitialGammaDistribution
    core_radius_fraction::Float64
    mu::Float64
    is_only_f_and_gamma_output::Bool

    function Solver(;
        aerodynamic_model_type::Model    = VSM,
        density::Float64                 = 1.225,
        max_iterations::Int64            = 1500,
        reltol::Float64           = 1e-5, # rel_err
        abstol::Float64           = 1e-5, # abs_err
        tol_reference_error::Float64     = 0.001,
        relaxation_factor::Float64       = 0.03,
        is_with_artificial_damping::Bool = false,
        artificial_damping::NamedTuple{(:k2, :k4), Tuple{Float64, Float64}}=(k2=0.1, k4=0.0),
        type_initial_gamma_distribution::InitialGammaDistribution=ELLIPTIC,
        core_radius_fraction::Float64    = 1e-20,
        mu::Float64                      = 1.81e-5, # TODO do not use magic constants
        is_only_f_and_gamma_output::Bool = false
    )
        new(
            aerodynamic_model_type,
            density,
            max_iterations,
            reltol,
            abstol,
            tol_reference_error,
            relaxation_factor,
            is_with_artificial_damping,
            artificial_damping,
            type_initial_gamma_distribution,
            core_radius_fraction,
            mu,
            is_only_f_and_gamma_output
        )
    end
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
- log=true: If true, print the number of iterations and other info.
- reference_point=zeros(MVec3)

# Returns
A dictionary with the results.
"""
function solve(solver::Solver, body_aero::BodyAerodynamics, gamma_distribution=nothing; 
               log=true, reference_point=zeros(MVec3))
    
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
        if solver.type_initial_gamma_distribution === ELLIPTIC
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
    # if !converged && relaxation_factor > 1e-3
    #     converged, gamma_new, alpha_array, v_a_array = gamma_loop(
    #         solver,
    #         body_aero,
    #         gamma_initial,
    #         va_array,
    #         chord_array,
    #         x_airf_array,
    #         y_airf_array,
    #         z_airf_array,
    #         panels,
    #         relaxation_factor/2;
    #         log
    #     )
    # end

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
    create_smooth_gamma!(gamma::Vector{Float64}, control_points::Vector{Float64}, span_positions::Vector{Float64})

Generate smooth gamma distribution using cubic B-spline interpolation with control points.
Ensures differentiability and bounded derivatives.
"""
function create_smooth_gamma!(gamma::Vector, control_points, span_positions)
    # Create normalized positions for control points
    n_controls = length(control_points)
    control_positions = range(-1.0, 1.0, n_controls)
    
    # Create cubic B-spline interpolation
    itp = CubicSplineInterpolation(control_positions, control_points)
    
    # Evaluate at panel positions
    for (i, pos) in enumerate(span_positions)
        gamma[i] = itp(pos)
    end
end

"""
    fit_control_points(gamma::Vector{Float64}, span_positions::Vector{Float64}, n_controls::Int)

Fit n_controls control points to a gamma distribution using B-spline interpolation.
Returns control point values that best represent the gamma distribution.

# Arguments
- `gamma`: Vector of circulation values
- `span_positions`: Vector of span positions where gamma is evaluated
- `n_controls`: Number of control points to fit

# Returns
- Vector{Float64}: Control point values
"""
function fit_control_points(gamma::Vector, span_positions, n_controls::Int)
    # Create control point positions (evenly spaced)
    control_positions = range(-1.0, 1.0, n_controls)
    
    # Create interpolation of existing gamma distribution
    itp = CubicSplineInterpolation(span_positions, gamma)
    
    # Evaluate at control points
    control_values = zeros(n_controls)
    for (i, pos) in enumerate(control_positions)
        control_values[i] = itp(pos)
    end
    
    return control_values
end

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

    # gamma_new, _ = smooth_circulation(gamma_new, 0.1, 0.03)
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

    n_controls = 13  # Odd number for symmetry
    span_positions = range(-1.0, 1.0, n_panels)

    function f!(res, control_points, p)
        create_smooth_gamma!(gamma, control_points, span_positions)

        # Calculate induced velocities
        mul!(velocity_view_x, AIC_x, gamma)
        mul!(velocity_view_y, AIC_y, gamma)
        mul!(velocity_view_z, AIC_z, gamma)
        
        relative_velocity_array .= va_array .+ induced_velocity_all
        for i in 1:n_panels
            relative_velocity_crossz[i, :] .= cross3(
                view(relative_velocity_array, i, :),
                view(z_airf_array, i, :)
            )
            v_acrossz_array[i, :] .= cross3(
                view(va_array, i, :),
                view(z_airf_array, i, :)
            )
        end

        for i in 1:n_panels
            v_normal_array[i] = dot(view(x_airf_array, i, :), view(relative_velocity_array, i, :))
            v_tangential_array[i] = dot(view(y_airf_array, i, :), view(relative_velocity_array, i, :))
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

        control_new = fit_control_points(gamma_new, span_positions, n_controls)
        res .= abs.(control_new .- control_points)

        # λ = 1e-2
        # for i in 2:length(res)-1
        #     res[i] += 0.5λ * ((res[i] - res[i-1]) - (res[i+1] - res[i]))
        # end

        return nothing
    end

    initial_controls = zeros(n_controls)
    prob = NonlinearProblem(f!, initial_controls, nothing)
    sol = NonlinearSolve.solve(prob, NewtonRaphson(autodiff=AutoFiniteDiff()); reltol=solver.reltol, abstol=solver.abstol, maxiters=solver.max_iterations)
    f!(initial_controls, sol.u, nothing)
    converged = NonlinearSolve.SciMLBase.successful_retcode(sol)
    @show sol.retcode

    if log && converged
        @info "Converged after $(sol.stats.nf) iterations"
    elseif log
        @warn "NO convergence after $(sol.stats.nf) iterations"
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
    # smoothness_threshold = smoothness_factor * circulation_mean

    # Calculate differences between adjacent points
    differences = diff(circulation[2:end-1])
    @debug "circulation_mean: $circulation_mean, diff: $differences"

    # Check smoothness
    if isempty(differences)
        return circulation, false
    end

    # if maximum(abs.(differences)) <= smoothness_threshold
    #     return circulation, false
    # end

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
    return smoothed, true
end