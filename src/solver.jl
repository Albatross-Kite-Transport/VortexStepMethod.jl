
"""
    VSMSolution

Struct for storing the solution of the [solve!](@ref) function. Must contain all info needed by `KiteModels.jl`.

# Attributes
- `panel_width_array`::Vector{Float64}: Width of the panels [m]
- cl_array::Vector{Float64}: Lift coefficients of the panels [-]
- cd_array::Vector{Float64}: Drag coefficients of the panels [-]
- cm_array::Vector{Float64}: Pitching moment coefficients of the panels [-]
- panel_lift::Vector{Float64}: Lift force of the panels [N]
- panel_drag::Vector{Float64}: Drag force of the panels [N]
- panel_moment::Vector{Float64}: Pitching moment around the spanwise vector of the panels [Nm]
- `f_body_3D`::Matrix{Float64}: Matrix of the aerodynamic forces (x, y, z vectors) [N]
- `m_body_3D`::Matrix{Float64}: Matrix of the aerodynamic moments [Nm]
- `gamma_distribution`::Union{Nothing, Vector{Float64}}: Vector containing the panel circulations.
- force::MVec3: Aerodynamic force vector in KB reference frame [N]
- moment::MVec3: Aerodynamic moments [Mx, My, Mz] around the reference point [Nm]
- force_coeffs::MVec3: Aerodynamic force coefficients [CFx, CFy, CFz] [-]
- `moment_coeffs`::MVec3: Aerodynamic moment coefficients [CMx, CMy, CMz] [-]
- `moment_dist`::Vector{Float64}: Pitching moments around the spanwise vector of each panel. [Nm]
- `moment_coeff_dist`::Vector{Float64}: Pitching moment coefficient around the spanwise vector of each panel. [-]
- `solver_status`::SolverStatus: enum, see [SolverStatus](@ref)
"""
@with_kw mutable struct VSMSolution{P,G}
    ### private vectors of solve_base!
    _x_airf_array::Matrix{Float64} = zeros(P, 3)
    _y_airf_array::Matrix{Float64} = zeros(P, 3)
    _z_airf_array::Matrix{Float64} = zeros(P, 3)
    _va_array::Matrix{Float64} = zeros(P, 3)
    _chord_array::Vector{Float64} = zeros(P)
    ### end of private vectors
    panel_width_array::Vector{Float64} = zeros(P)
    cl_array::Vector{Float64} = zeros(P)
    cd_array::Vector{Float64} = zeros(P)
    cm_array::Vector{Float64} = zeros(P)
    panel_lift::Vector{Float64} = zeros(P)
    panel_drag::Vector{Float64} = zeros(P)
    panel_moment::Vector{Float64} = zeros(P)
    f_body_3D::Matrix{Float64} = zeros(3, P)
    m_body_3D::Matrix{Float64} = zeros(3, P)
    gamma_distribution::Union{Nothing, Vector{Float64}} = nothing
    force::MVec3 = zeros(MVec3)          
    moment::MVec3 = zeros(MVec3)       
    force_coeffs::MVec3 = zeros(MVec3)  
    moment_coeffs::MVec3 = zeros(MVec3)  
    moment_dist::MVector{P, Float64} = zeros(MVector{P, Float64})
    moment_coeff_dist::MVector{P, Float64} = zeros(MVector{P, Float64})
    group_moment_dist::MVector{G, Float64} = zeros(MVector{G, Float64})
    group_moment_coeff_dist::MVector{G, Float64} = zeros(MVector{G, Float64})
    solver_status::SolverStatus = FAILURE
end

# Output of the function gamma_loop!
@with_kw mutable struct LoopResult{P}
    converged::Bool              = false
    gamma_new::MVector{P, Float64}   = zeros(MVector{P, Float64})
    alpha_array::MVector{P, Float64} = zeros(MVector{P, Float64}) # TODO: Is this different from BodyAerodynamics.alpha_array ?
    v_a_array::MVector{P, Float64}   = zeros(MVector{P, Float64})
end

@with_kw struct BaseResult{P}
    va_norm_array::MVector{P, Float64} = zeros(MVector{P, Float64})
    va_unit_array::Matrix{Float64} = zeros(P, 3)
end

"""
    Solver

Main solver structure for the Vortex Step Method.See also: [solve](@ref)

# Attributes

## General settings
- `aerodynamic_model_type`::Model = VSM: The model type, see: [Model](@ref)
- density::Float64 = 1.225: Air density [kg/m³] 
- `max_iterations`::Int64 = 1500
- `rtol`::Float64 = 1e-5: relative error
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

## Solution
sol::VSMSolution = VSMSolution(): The result of calling [solve!](@ref) 
"""
@with_kw mutable struct Solver{P,G}
    # General settings
    solver_type::SolverType = LOOP
    aerodynamic_model_type::Model = VSM
    density::Float64 = 1.225
    max_iterations::Int64 = 1500
    rtol::Float64 = 1e-5
    tol_reference_error::Float64 = 0.001
    relaxation_factor::Float64 = 0.03

    # Nonlin solver fields
    prob::Union{NonlinearProblem, Nothing} = nothing
    nonlin_cache::Union{NonlinearSolve.AbstractNonlinearSolveCache, Nothing} = nothing
    atol::Float64 = 1e-5
    
    # Damping settings
    is_with_artificial_damping::Bool = false
    artificial_damping::NamedTuple{(:k2, :k4), Tuple{Float64, Float64}} =(k2=0.1, k4=0.0)
    
    # Additional settings
    type_initial_gamma_distribution::InitialGammaDistribution = ZEROS
    core_radius_fraction::Float64 = 1e-20
    mu::Float64 = 1.81e-5
    is_only_f_and_gamma_output::Bool = false

    # Intermediate results
    lr::LoopResult{P} = LoopResult{P}()
    br::BaseResult{P} = BaseResult{P}()
    cache::Vector{PreallocationTools.LazyBufferCache{typeof(identity), typeof(identity)}} = [LazyBufferCache() for _ in 1:11]
    cache_base::Vector{PreallocationTools.LazyBufferCache{typeof(identity), typeof(identity)}}  = [LazyBufferCache()]
    cache_solve::Vector{PreallocationTools.LazyBufferCache{typeof(identity), typeof(identity)}} = [LazyBufferCache()]
    cache_lin::Vector{PreallocationTools.LazyBufferCache{typeof(identity), typeof(identity)}} = [LazyBufferCache() for _ in 1:4]
    
    # Solution
    sol::VSMSolution{P,G} = VSMSolution{P,G}()
end

function Solver(body_aero; kwargs...)
    P = length(body_aero.panels)
    G = sum([wing.n_groups for wing in body_aero.wings])
    return Solver{P,G}(; kwargs...)
end

"""
    solve!(solver::Solver, body_aero::BodyAerodynamics, gamma_distribution=solver.sol.gamma_distribution; 
          log=false, reference_point=zeros(MVec3), moment_frac=0.1)

Main solving routine for the aerodynamic model. Reference point is in the kite body (KB) frame.
This version is modifying the `solver.sol` struct and is faster than the `solve` function which returns
a dictionary.

# Arguments:
- solver::Solver: The solver to use, could be a VSM or LLT solver. See: [Solver](@ref)
- body_aero::BodyAerodynamics: The aerodynamic body. See: [BodyAerodynamics](@ref)
- gamma_distribution: Initial circulation vector or nothing; Length: Number of segments. [m²/s]

# Keyword Arguments:
- log=false: If true, print the number of iterations and other info.
- reference_point=zeros(MVec3)
- moment_frac=0.1: X-coordinate of normalized panel around which the moment distribution should be calculated.

# Returns
The solution of type [VSMSolution](@ref)
"""
function solve!(solver::Solver, body_aero::BodyAerodynamics, gamma_distribution=solver.sol.gamma_distribution; 
        log=false, reference_point=zeros(MVec3), moment_frac=0.1)

    # calculate intermediate result
    solve_base!(solver, body_aero, gamma_distribution; log, reference_point)
    gamma_new = solver.lr.gamma_new
    if !isnothing(solver.sol.gamma_distribution)
        solver.sol.gamma_distribution .= gamma_new
    else
        solver.sol.gamma_distribution = gamma_new
    end

    # Initialize arrays
    cl_array = solver.sol.cl_array
    cd_array = solver.sol.cd_array
    cm_array = solver.sol.cm_array
    converged = solver.lr.converged
    alpha_array = solver.lr.alpha_array
    alpha_corrected = solver.cache_solve[1][alpha_array]
    v_a_array = solver.lr.v_a_array
    panels = body_aero.panels
   
    panel_width_array = solver.sol.panel_width_array
    solver.sol.moment_dist .= 0
    solver.sol.moment_coeff_dist .= 0
    moment_dist = solver.sol.moment_dist
    moment_coeff_dist = solver.sol.moment_coeff_dist
    density = solver.density
    aerodynamic_model_type = solver.aerodynamic_model_type

    # Calculate coefficients for each panel
    for (i, panel) in enumerate(panels)                                               # zero bytes
        cl_array[i] = calculate_cl(panel, alpha_array[i])
        cd_array[i], cm_array[i] = calculate_cd_cm(panel, alpha_array[i])
        panel_width_array[i] = panel.width
    end

    # create an alias for the three vertical output vectors
    lift = solver.sol.panel_lift
    drag = solver.sol.panel_drag
    panel_moment = solver.sol.panel_moment

    # Compute using fused broadcasting (no intermediate allocations)
    @. lift = cl_array * 0.5 * density * v_a_array^2 * solver.sol._chord_array
    @. drag = cd_array * 0.5 * density * v_a_array^2 * solver.sol._chord_array
    @. panel_moment = cm_array * 0.5 * density * v_a_array^2 * solver.sol._chord_array

    # Calculate alpha corrections based on model type
    if aerodynamic_model_type == VSM                             # 64 bytes
        update_effective_angle_of_attack!(
            alpha_corrected,
            body_aero,
            gamma_new,
            solver.core_radius_fraction,
            solver.sol._z_airf_array,
            solver.sol._x_airf_array,
            solver.sol._va_array,
            solver.br.va_norm_array,
            solver.br.va_unit_array
        )
    elseif aerodynamic_model_type == LLT
        alpha_corrected .= alpha_array
    end

    # Initialize result arrays
    area_all_panels = 0.0

    # Get wing properties
    spanwise_direction = body_aero.wings[1].spanwise_direction
    va_mag = norm(body_aero.va)
    q_inf = 0.5 * density * va_mag^2

    # Calculate wing geometry properties
    projected_area = body_aero.projected_area
    
    for (i, panel) in enumerate(panels)                                               # 8000 bytes

        ### Lift and Drag ###
        # Panel geometry
        panel_area = panel.chord * panel.width
        area_all_panels += panel_area

        # Calculate induced velocity direction
        alpha_corrected_i = alpha_corrected[i]
        dir_induced_va_airfoil = cos(alpha_corrected_i) * panel.x_airf + 
                                 sin(alpha_corrected_i) * panel.z_airf
        normalize!(dir_induced_va_airfoil)

        # Calculate lift and drag directions
        dir_lift_induced_va = dir_induced_va_airfoil × panel.y_airf
        normalize!(dir_lift_induced_va)
        dir_drag_induced_va = spanwise_direction × dir_lift_induced_va
        normalize!(dir_drag_induced_va)

        # Calculate force vectors
        lift_induced_va = lift[i] * dir_lift_induced_va
        drag_induced_va = drag[i] * dir_drag_induced_va
        ftotal_induced_va = lift_induced_va + drag_induced_va

        # Body frame forces
        solver.sol.f_body_3D[:,i] .= ftotal_induced_va .* panel.width

        # Calculate the moments
        # (1) Panel aerodynamic center in body frame:
        panel_ac_body = panel.aero_center  # 3D [x, y, z]
        # (2) Convert local (2D) pitching moment to a 3D vector in body coords.
        #     Use the axis around which the moment is defined,
        #     which is the y-axis pointing "spanwise"
        moment_axis_body = panel.y_airf
        M_local_3D = panel_moment[i] * moment_axis_body * panel.width
        # Vector from panel AC to the chosen reference point:
        r_vector = panel_ac_body - reference_point  # e.g. CG, wing root, etc.
        # Cross product to shift the force from panel AC to ref. point:
        M_shift = r_vector × MVec3(solver.sol.f_body_3D[:,i])
        # Total panel moment about the reference point:
        solver.sol.m_body_3D[:,i] .= M_local_3D + M_shift

        # Calculate the moment distribution (moment on each panel)
        arm = (moment_frac - 0.25) * panel.chord
        moment_dist[i] = ((ftotal_induced_va ⋅ panel.z_airf) * arm + panel_moment[i]) * panel.width
        moment_coeff_dist[i] = moment_dist[i] / (q_inf * projected_area)
    end

    group_moment_dist = solver.sol.group_moment_dist
    group_moment_coeff_dist = solver.sol.group_moment_coeff_dist
    group_moment_dist .= 0.0
    group_moment_coeff_dist .= 0.0
    panel_idx = 1
    group_idx = 1
    for wing in body_aero.wings
        panels_per_group = wing.n_panels ÷ wing.n_groups
        for _ in 1:wing.n_groups
            for _ in 1:panels_per_group
                group_moment_dist[group_idx] += moment_dist[panel_idx]
                group_moment_coeff_dist[group_idx] += moment_coeff_dist[panel_idx]
                panel_idx += 1
            end
            group_idx += 1
        end
    end

    # update the result struct
    solver.sol.force .= MVec3(
        sum(solver.sol.f_body_3D[1,:]),
        sum(solver.sol.f_body_3D[2,:]),
        sum(solver.sol.f_body_3D[3,:])
    )
    solver.sol.moment .= MVec3(
        sum(solver.sol.m_body_3D[1,:]),
        sum(solver.sol.m_body_3D[2,:]),
        sum(solver.sol.m_body_3D[3,:])
    )
    solver.sol.force_coeffs .= solver.sol.force ./ (q_inf * projected_area)
    solver.sol.moment_coeffs .= solver.sol.moment ./ (q_inf * projected_area)
    if converged
        # TODO: Check if the result if feasible if converged
        solver.sol.solver_status = FEASIBLE
    else
        solver.sol.solver_status = FAILURE
    end

    return solver.sol
end

"""
    solve(solver::Solver, body_aero::BodyAerodynamics, gamma_distribution=nothing; 
          log=false, reference_point=zeros(MVec3))

Main solving routine for the aerodynamic model. Reference point is in the kite body (KB) frame.
See also: [solve!](@ref)

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
    # calculate intermediate result
    solve_base!(solver, body_aero, gamma_distribution; log, reference_point)

    # Calculate final results as dictionary
    results = calculate_results(
        body_aero,
        solver.lr.gamma_new,
        reference_point,
        solver.density,
        solver.aerodynamic_model_type,
        solver.core_radius_fraction,
        solver.mu,
        solver.lr.alpha_array,
        solver.lr.v_a_array,
        solver.sol._chord_array,
        solver.sol._x_airf_array,
        solver.sol._y_airf_array,
        solver.sol._z_airf_array,
        solver.sol._va_array,
        solver.br.va_norm_array,
        solver.br.va_unit_array,
        body_aero.panels,
        solver.is_only_f_and_gamma_output
    )
    return results
end

@inline @inbounds function calc_norm_array!(va_norm_array, va_array)
    for i in 1:size(va_array, 1)
        va_norm_array[i] = norm(MVec3(view(va_array, i, :)))
    end
end

function solve_base!(solver::Solver, body_aero::BodyAerodynamics, gamma_distribution=nothing; 
               log=false, reference_point=zeros(MVec3))
    
    # check arguments
    isnothing(body_aero.panels[1].va) && throw(ArgumentError("Inflow conditions are not set, use set_va!(body_aero, va)"))
    
    # Initialize variables
    panels = body_aero.panels
    n_panels = length(panels)
    relaxation_factor = solver.relaxation_factor
    
    # Clear arrays
    solver.sol._x_airf_array .= 0
    solver.sol._y_airf_array .= 0
    solver.sol._z_airf_array .= 0
    solver.sol._va_array .= 0
    solver.sol._chord_array .= 0

    # Fill arrays from panels
    for (i, panel) in enumerate(panels)
        solver.sol._x_airf_array[i, :] .= panel.x_airf
        solver.sol._y_airf_array[i, :] .= panel.y_airf
        solver.sol._z_airf_array[i, :] .= panel.z_airf
        solver.sol._va_array[i, :] .= panel.va
        solver.sol._chord_array[i] = panel.chord
    end

    # Calculate unit vectors
    calc_norm_array!(solver.br.va_norm_array, solver.sol._va_array)
    solver.br.va_unit_array .= solver.sol._va_array ./ solver.br.va_norm_array

    # Calculate AIC matrices
    calculate_AIC_matrices!(body_aero, solver.aerodynamic_model_type, solver.core_radius_fraction, solver.br.va_norm_array,
                            solver.br.va_unit_array)

    # Initialize gamma distribution
    gamma_initial = solver.cache_base[1][solver.sol._chord_array]
    if isnothing(gamma_distribution)
        if solver.type_initial_gamma_distribution == ELLIPTIC
            calculate_circulation_distribution_elliptical_wing(gamma_initial, body_aero)
        else
            gamma_initial .= 0
        end
    else
        length(gamma_distribution) == n_panels || 
            throw(ArgumentError("gamma_distribution length must match number of panels"))
        gamma_initial .= gamma_distribution
    end

    @debug "Initial gamma_new: $gamma_initial"
    solver.lr.gamma_new .= gamma_initial
    # Run main iteration loop
    gamma_loop!(solver, body_aero, panels, relaxation_factor; log)
    # Try again with reduced relaxation factor if not converged
    if solver.solver_type == LOOP && !solver.lr.converged && relaxation_factor > 1e-3
        log && @warn "Running again with half the relaxation_factor = $(relaxation_factor/2)"
        solver.lr.gamma_new .= gamma_initial
        gamma_loop!(solver, body_aero, panels, relaxation_factor/2; log)
    end

    nothing
end

"""
    gamma_loop!(solver::Solver, AIC_x::Matrix{Float64}, 
              AIC_y::Matrix{Float64}, AIC_z::Matrix{Float64},
              panels::Vector{Panel}, relaxation_factor::Float64; log=true)

Main iteration loop for calculating circulation distribution.
"""
function gamma_loop!(
    solver::Solver,
    body_aero::BodyAerodynamics,
    panels::Vector{Panel},
    relaxation_factor::Float64;
    log::Bool = true
)
    va_array = solver.sol._va_array
    chord_array = solver.sol._chord_array
    x_airf_array = solver.sol._x_airf_array
    y_airf_array = solver.sol._y_airf_array
    z_airf_array = solver.sol._z_airf_array
    solver.lr.converged   = false
    n_panels    = length(body_aero.panels)
    solver.lr.alpha_array .= body_aero.alpha_array
    solver.lr.v_a_array   .= body_aero.v_a_array
    
    va_magw_array            = solver.cache[1][solver.lr.v_a_array]
    gamma                    = solver.cache[2][solver.lr.gamma_new]
    abs_gamma_new            = solver.cache[3][solver.lr.gamma_new]
    induced_velocity_all     = solver.cache[4][va_array]
    relative_velocity_array  = solver.cache[5][va_array]
    relative_velocity_crossz = solver.cache[6][va_array]
    v_acrossz_array          = solver.cache[7][va_array]
    cl_array                 = solver.cache[8][solver.lr.gamma_new]
    damp                     = solver.cache[9][solver.lr.gamma_new]
    v_normal_array           = solver.cache[10][solver.lr.gamma_new]
    v_tangential_array       = solver.cache[11][solver.lr.gamma_new]

    AIC_x, AIC_y, AIC_z = body_aero.AIC[1, :, :], body_aero.AIC[2, :, :], body_aero.AIC[3, :, :]

    velocity_view_x = @view induced_velocity_all[:, 1]
    velocity_view_y = @view induced_velocity_all[:, 2]
    velocity_view_z = @view induced_velocity_all[:, 3]

    function calc_gamma_new!(gamma_new, gamma, p)
        # (AIC_x, AIC_y, AIC_z, va_array, induced_velocity_all, 
        #     x_airf_array, y_airf_array, z_airf_array, panels, chord_array) = p
        # Calculate induced velocities
        mul!(velocity_view_x, AIC_x, gamma)
        mul!(velocity_view_y, AIC_y, gamma)
        mul!(velocity_view_z, AIC_z, gamma)
        
        relative_velocity_array .= va_array .+ induced_velocity_all
        for i in 1:n_panels
            relative_velocity_crossz[i, :] .=  MVec3(view(relative_velocity_array, i, :)) ×
                                            MVec3(view(y_airf_array, i, :))
            v_acrossz_array[i, :]          .=  MVec3(view(va_array, i, :)) ×
                                            MVec3(view(y_airf_array, i, :))
        end

        for i in 1:n_panels
            v_normal_array[i] = view(z_airf_array, i, :) ⋅ view(relative_velocity_array, i, :)
            v_tangential_array[i] = view(x_airf_array, i, :) ⋅ view(relative_velocity_array, i, :)
        end
        solver.lr.alpha_array .= atan.(v_normal_array, v_tangential_array)

        for i in 1:n_panels
            @views solver.lr.v_a_array[i] = norm(relative_velocity_crossz[i, :])
            @views va_magw_array[i] = norm(v_acrossz_array[i, :])
        end
        
        for (i, (panel, alpha)) in enumerate(zip(panels, solver.lr.alpha_array))
            cl_array[i] = calculate_cl(panel, alpha)
        end
        gamma_new .= 0.5 .* solver.lr.v_a_array.^2 ./ va_magw_array .* cl_array .* chord_array
        nothing
    end

    if solver.solver_type == NONLIN
        if isnothing(solver.prob)
            function f_nonlin!(d_gamma, gamma, p)
                calc_gamma_new!(solver.lr.gamma_new, gamma, p)
                d_gamma .= solver.lr.gamma_new .- gamma
                nothing
            end
            solver.prob = NonlinearProblem(f_nonlin!, solver.lr.gamma_new, nothing)
            solver.nonlin_cache = init(solver.prob, NewtonRaphson(autodiff=AutoFiniteDiff()); abstol=solver.atol, reltol=solver.rtol)
        end
        
        sol = NonlinearSolve.solve!(solver.nonlin_cache)
        gamma .= sol.u
        solver.lr.gamma_new .= sol.u
        solver.lr.converged = SciMLBase.successful_retcode(sol)
        return nothing
    end

    if solver.solver_type == LOOP
        function f_loop!(gamma_new, gamma, damp)
            gamma .= gamma_new
            calc_gamma_new!(gamma_new, gamma, nothing)
            # Update gamma with relaxation and damping
            @. gamma_new = (1 - relaxation_factor) * gamma + 
                    relaxation_factor * gamma_new + damp
            
            # Apply damping if needed
            if solver.is_with_artificial_damping
                is_damping_applied = smooth_circulation!(damp, gamma, 0.1, 0.5)
                @debug "damp: $damp"
            else
                damp .= 0.0
                is_damping_applied = false
            end
            nothing
        end
        iters = 0
        for i in 1:solver.max_iterations
            iters += 1

            f_loop!(solver.lr.gamma_new, gamma, damp)

            # Check convergence
            abs_gamma_new .= abs.(solver.lr.gamma_new)
            reference_error = maximum(abs_gamma_new)
            reference_error = max(reference_error, solver.tol_reference_error)
            abs_gamma_new .= abs.(solver.lr.gamma_new .- gamma)
            error = maximum(abs_gamma_new)
            normalized_error = error / reference_error

            @debug "Iteration: $i, normalized_error: $normalized_error, is_damping_applied: $is_damping_applied"

            if normalized_error < solver.rtol
                solver.lr.converged = true
                break
            end
        end
    
        if log && solver.lr.converged
            @info "Converged after $iters iterations"
        elseif log
            @warn "NO convergence after $(solver.max_iterations) iterations"
        end
        return nothing
    end
end

"""
    smooth_circulation!(damp, circulation, 
                      smoothness_factor::Float64, 
                      damping_factor::Float64)

Smooth circulation distribution if needed.

Returns:
- Tuple of smoothed circulation and boolean indicating if smoothing was applied
"""
function smooth_circulation!(
    damp,
    circulation,
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

    damp .= smoothed .- circulation
    return true
end

"""
    linearize(solver::Solver, body_aero::BodyAerodynamics, wing::RamAirWing, y::Vector{T}; 
        theta_idxs=1:4, delta_idxs=nothing, va_idxs=nothing, omega_idxs=nothing, kwargs...) where T

Compute the Jacobian matrix for a ram air wing around an operating point using finite differences.

# Arguments
- `solver`: VSM solver instance (must be initialized)
- `body_aero`: Aerodynamic body representation
- `wing`: RamAirWing model to linearize
- `y`: Input vector at operating point, containing a combination of control angles and velocities

# Keyword Arguments
- `theta_idxs`: Indices of twist angles in input vector (default: 1:4)
- `delta_idxs`: Indices of trailing edge deflection angles (default: nothing)
- `va_idxs`: Indices of velocity components `[vx, vy, vz]` (default: nothing)
- `omega_idxs`: Indices of angular velocity components `[ωx, ωy, ωz]` (default: nothing)
- `kwargs...`: Additional arguments passed to the `solve!` function

# Returns
- `jac`: Jacobian matrix (∂outputs/∂inputs)
- `results`: Output vector at the operating point [Fx, Fy, Fz, Mx, My, Mz, group_moments...]

# Example
```julia
# Initialize wing and solver
wing = RamAirWing("path/to/body.obj", "path/to/foil.dat")
body_aero = BodyAerodynamics([wing])
solver = Solver(body_aero)

# Define operating point with 4 control angles, velocity, and angular rates
y_op = [zeros(4);        # 4 twist control angles (rad)
       [15.0, 0.0, 0.0]; # Velocity vector (m/s)
       zeros(3)]         # Angular velocity (rad/s)

# Compute Jacobian
jac, results = linearize(
    solver, body_aero, wing, y_op;
    theta_idxs=1:4,      # Twist angles
    va_idxs=5:7,         # Velocity components
    omega_idxs=8:10      # Angular rates
)
```
"""
function linearize(solver::Solver, body_aero::BodyAerodynamics, y::Vector{T}; 
        theta_idxs=1:4, 
        delta_idxs=nothing,
        va_idxs=nothing,
        omega_idxs=nothing,
        kwargs...) where T

    !(length(body_aero.wings) == 1) && throw(ArgumentError("Linearization only works for a body_aero with one wing"))
    wing = body_aero.wings[1]

    init_va = body_aero.cache[1][body_aero.va]
    init_va .= body_aero.va
    if !isnothing(theta_idxs)
        @views last_theta::Vector{T} = body_aero.cache[2][y[theta_idxs]]
        last_theta .= 0.0
    end
    if !isnothing(delta_idxs)
        @views last_delta::Vector{T} = body_aero.cache[3][y[delta_idxs]]
        last_delta .= 0.0
    end

    # Function to compute forces for given control inputs
    function calc_results!(results, y)
        @views theta_angles = isnothing(theta_idxs) ? nothing : y[theta_idxs]
        @views delta_angles = isnothing(delta_idxs) ? nothing : y[delta_idxs]

        if !isnothing(theta_angles) && isnothing(delta_angles)
            if !all(theta_angles .== last_theta)
                VortexStepMethod.group_deform!(wing, theta_angles, nothing; smooth=false)
                VortexStepMethod.init!(body_aero; init_aero=false)
                last_theta .= theta_angles
            end
        elseif !isnothing(theta_angles) && !isnothing(delta_angles)
            if !all(delta_angles .== last_delta) || !all(theta_angles .== last_theta)
                VortexStepMethod.group_deform!(wing, theta_angles, delta_angles; smooth=false)
                VortexStepMethod.init!(body_aero; init_aero=false)
                last_theta .= theta_angles
                last_delta .= delta_angles
            end
        elseif isnothing(theta_angles) && !isnothing(delta_angles)
            if !all(delta_angles .== last_delta)
                VortexStepMethod.group_deform!(wing, nothing, delta_angles; smooth=false)
                VortexStepMethod.init!(body_aero; init_aero=false)
                last_delta .= delta_angles
            end
        end

        if !isnothing(va_idxs) && isnothing(omega_idxs)
            set_va!(body_aero, y[va_idxs])
        elseif !isnothing(va_idxs) && !isnothing(omega_idxs)
            set_va!(body_aero, y[va_idxs], y[omega_idxs])
        elseif isnothing(va_idxs) && isnothing(omega_idxs)
            set_va!(body_aero, init_va)
        end

        solve!(solver, body_aero; kwargs...)
        results[1:3] .= solver.sol.force
        results[4:6] .= solver.sol.moment
        results[7:end] .= solver.sol.group_moment_dist
        return nothing
    end

    results = zeros(3+3+length(solver.sol.group_moment_dist))
    jac = zeros(length(results), length(y))
    backend = AutoFiniteDiff(absstep=1e2solver.atol, relstep=1e2solver.rtol)
    prep = prepare_jacobian(calc_results!, results, backend, y)
    jacobian!(calc_results!, results, jac, prep, backend, y)

    calc_results!(results, y)
    return jac, results
end

