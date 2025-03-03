"""
    @with_kw mutable struct BodyAerodynamics{P}

Main structure for calculating aerodynamic properties of bodies.

# Fields
- panels::Vector{Panel}: Vector of [Panel](@ref) structs
- wings::Vector{AbstractWing}:   A vector of wings; a body can have multiple wings
- `_va`::MVec3 = zeros(MVec3):   A vector of the apparent wind speed, see: [MVec3](@ref)
- `omega`::MVec3 = zeros(MVec3): A vector of the turn rates around the kite body axes
- `gamma_distribution`::Vector{Float64}=zeros(Float64, P): A vector of the circulation 
                        of the velocity field; Length: Number of segments. [m²/s]
- `alpha_uncorrected`::Vector{Float64}=zeros(Float64, P): unclear, please define
- `alpha_corrected`::Vector{Float64}=zeros(Float64, P):   unclear, please define
- `stall_angle_list`::Vector{Float64}=zeros(Float64, P):  unclear, please define
- alpha_array::Vector{Float64} = zeros(Float64, P)
- v_a_array::Vector{Float64} = zeros(Float64, P)
- work_vectors::NTuple{10,MVec3} = ntuple(_ -> zeros(MVec3), 10)
- AIC::Array{Float64, 3} = zeros(3, P, P)
"""
@with_kw mutable struct BodyAerodynamics{P}
    panels::Vector{Panel}
    wings::Vector{AbstractWing}
    _va::MVec3 = zeros(MVec3)
    omega::MVec3 = zeros(MVec3)
    gamma_distribution::Vector{Float64} = zeros(Float64, P)
    alpha_uncorrected::Vector{Float64} = zeros(Float64, P)
    alpha_corrected::Vector{Float64} = zeros(Float64, P)
    stall_angle_list::Vector{Float64} = zeros(Float64, P)
    alpha_array::Vector{Float64} = zeros(Float64, P)
    v_a_array::Vector{Float64} = zeros(Float64, P)
    work_vectors::NTuple{10,MVec3} = ntuple(_ -> zeros(MVec3), 10)
    AIC::Array{Float64, 3} = zeros(3, P, P)
end

function BodyAerodynamics(
    wings::Vector{T};
    aero_center_location=0.25,
    control_point_location=0.75
) where T <: AbstractWing
    # Initialize panels
    panels = Panel[]
    for wing in wings
        if wing.spanwise_panel_distribution == UNCHANGED
            wing.n_panels = length(wing.sections) - 1
            wing.refined_sections = wing.sections
        else
            wing.refined_sections = Section[Section() for _ in 1:wing.n_panels+1]
        end

        # Create panels
        for _ in 1:wing.n_panels
            panel = Panel()
            push!(panels, panel)
        end
    end

    body_aero = BodyAerodynamics{length(panels)}(; panels, wings)
    init!(body_aero; aero_center_location, control_point_location)
    return body_aero
end

function init!(body_aero::BodyAerodynamics; 
    aero_center_location=0.25,
    control_point_location=0.75)

    idx = 1
    for wing in body_aero.wings
        panel_props = init!(wing; aero_center_location, control_point_location)
        
        # Create panels
        for i in 1:wing.n_panels
            if wing isa KiteWing
                beta = wing.beta_dist[i]
            else
                beta = 0.0
            end
            init!(
                body_aero.panels[idx], 
                wing.refined_sections[i],
                wing.refined_sections[i+1],
                panel_props.aero_centers[i],
                panel_props.control_points[i],
                panel_props.bound_points_1[i],
                panel_props.bound_points_2[i],
                panel_props.x_airf[i],
                panel_props.y_airf[i],
                panel_props.z_airf[i],
                beta
            )
            idx += 1
        end
    end

    @assert all([panel.filaments[1].initialized for panel in body_aero.panels])
    
    # Initialize rest of the struct
    body_aero.stall_angle_list .= calculate_stall_angle_list(body_aero.panels)
    body_aero.alpha_array .= 0.0
    body_aero.v_a_array .= 0.0 
    body_aero.AIC .= 0.0
    return nothing
end

function Base.getproperty(obj::BodyAerodynamics, sym::Symbol)
    if sym === :va
        return getfield(obj, :_va)
    end
    return getfield(obj, sym)
end

function Base.setproperty!(obj::BodyAerodynamics, sym::Symbol, val)
    if sym === :va || sym === :omega
        set_va!(obj, val)
    else
        setfield!(obj, sym, val)
    end
end

"""
    PanelProperties

Structure to hold calculated panel properties.

# Fields
- `aero_centers`::Vector{MVec3}
- `control_points`::Vector{MVec3}
- `bound_points_1`::Vector{MVec3}
- `bound_points_2`::Vector{MVec3}
- `x_airf`::Vector{MVec3}: Vector of unit vectors perpendicular to chord line
- `y_airf`::Vector{MVec3}: Vector of unit vectors parallel to chord line
- `z_airf`::Vector{MVec3}: Vector of unit vectors in spanwise direction
"""
struct PanelProperties
    aero_centers::Vector{MVec3}
    control_points::Vector{MVec3}
    bound_points_1::Vector{MVec3}
    bound_points_2::Vector{MVec3}
    x_airf::Vector{MVec3}
    y_airf::Vector{MVec3}
    z_airf::Vector{MVec3}
end

"""
    calculate_panel_properties(section_list::Vector{Section}, n_panels::Int,
                             aero_center_loc::Float64, control_point_loc::Float64)

Calculate geometric properties for each panel.

Returns:
    PanelProperties containing vectors for each property
"""
function calculate_panel_properties(section_list::Vector{Section}, n_panels::Int,
                                  aero_center_loc::Float64, control_point_loc::Float64)
    # Initialize arrays
    aero_centers = MVec3[]
    control_points = MVec3[]
    bound_points_1 = MVec3[]
    bound_points_2 = MVec3[]
    x_airf = MVec3[]
    y_airf = MVec3[]
    z_airf = MVec3[]
    
    # Define coordinates matrix
    coords = zeros(2 * (n_panels + 1), 3)
    @debug "Shape of coordinates: $(size(coords))"
    
    for i in 1:n_panels
        coords[2i-1, :] .= section_list[i].LE_point
        coords[2i, :]   .= section_list[i].TE_point
        coords[2i+1, :] .= section_list[i+1].LE_point
        coords[2i+2, :] .= section_list[i+1].TE_point
    end
    
    @debug "Coordinates: $coords"
    
    for i in 1:n_panels
        # Define panel points
        section = Dict(
            "p1" => coords[2i-1, :],     # LE_1
            "p2" => coords[2i+1, :],     # LE_2
            "p3" => coords[2i+2, :],     # TE_2
            "p4" => coords[2i, :]        # TE_1
        )
        
        # Calculate control point position
        di = norm(coords[2i-1, :] * 0.75 + coords[2i, :] * 0.25 - 
                 (coords[2i+1, :] * 0.75 + coords[2i+2, :] * 0.25))
        
        ncp = if i == 1
            diplus = norm(coords[2i+1, :] * 0.75 + coords[2i+2, :] * 0.25 - 
                         (coords[2i+3, :] * 0.75 + coords[2i+4, :] * 0.25))
            di / (di + diplus)
        elseif i == n_panels
            dimin = norm(coords[2i-3, :] * 0.75 + coords[2i-2, :] * 0.25 - 
                        (coords[2i-1, :] * 0.75 + coords[2i, :] * 0.25))
            dimin / (dimin + di)
        else
            dimin = norm(coords[2i-3, :] * 0.75 + coords[2i-2, :] * 0.25 - 
                        (coords[2i-1, :] * 0.75 + coords[2i, :] * 0.25))
            diplus = norm(coords[2i+1, :] * 0.75 + coords[2i+2, :] * 0.25 - 
                         (coords[2i+3, :] * 0.75 + coords[2i+4, :] * 0.25))
            0.25 * (dimin / (dimin + di) + di / (di + diplus) + 1)
        end
        
        ncp = 1 - ncp
        
        # Calculate points
        LL_point = (section["p2"] * (1 - ncp) + section["p1"] * ncp) * 0.75 +
                   (section["p3"] * (1 - ncp) + section["p4"] * ncp) * 0.25
        
        VSM_point = (section["p2"] * (1 - ncp) + section["p1"] * ncp) * 0.25 +
                    (section["p3"] * (1 - ncp) + section["p4"] * ncp) * 0.75
        
        bound_1 = section["p1"] * 0.75 + section["p4"] * 0.25
        bound_2 = section["p2"] * 0.75 + section["p3"] * 0.25
        
        # Calculate reference frame vectors
        z_airf_vec = cross(VSM_point - LL_point, section["p1"] - section["p2"])
        z_airf_vec = z_airf_vec / norm(z_airf_vec)
        
        x_airf_vec = VSM_point - LL_point
        x_airf_vec = x_airf_vec / norm(x_airf_vec)
        
        y_airf_vec = bound_1 - bound_2
        y_airf_vec = y_airf_vec / norm(y_airf_vec)
        
        # Store results
        push!(aero_centers, LL_point)
        push!(control_points, VSM_point)
        push!(bound_points_1, bound_1)
        push!(bound_points_2, bound_2)
        push!(x_airf, x_airf_vec)
        push!(y_airf, y_airf_vec)
        push!(z_airf, z_airf_vec)
    end
    
    return PanelProperties(aero_centers, control_points, bound_points_1, 
                          bound_points_2, x_airf, y_airf, z_airf)
end

"""
    calculate_AIC_matrices!(body_aero::BodyAerodynamics, model::Model, 
                         core_radius_fraction::Float64,
                         va_norm_array::Vector{Float64}, 
                         va_unit_array::Matrix{Float64})

Calculate Aerodynamic Influence Coefficient matrices.

See also: [BodyAerodynamics](@ref), [Model](@ref)

Returns:
    Tuple of (`AIC_x`, `AIC_y`, `AIC_z`) matrices
"""
function calculate_AIC_matrices!(body_aero::BodyAerodynamics, model::Model,
                              core_radius_fraction::Float64,
                              va_norm_array::Vector{Float64}, 
                              va_unit_array::Matrix{Float64})
    # Determine evaluation point based on model
    evaluation_point = model === VSM ? :control_point : :aero_center
    evaluation_point_on_bound = model === LLT
    
    # Initialize AIC matrices
    velocity_induced, tempvel, va_unit, U_2D = zeros(MVec3), zeros(MVec3), zeros(MVec3), zeros(MVec3)
    
    # Calculate influence coefficients
    for icp in eachindex(body_aero.panels)
        ep = getproperty(body_aero.panels[icp], evaluation_point)
        for jring in eachindex(body_aero.panels)
            va_unit .= @views va_unit_array[jring, :]
            filaments = body_aero.panels[jring].filaments
            va_norm = va_norm_array[jring]
            calculate_velocity_induced_single_ring_semiinfinite!(
                velocity_induced,
                tempvel,
                filaments,
                ep,
                evaluation_point_on_bound,
                va_norm,
                va_unit,
                1.0,
                core_radius_fraction,
                body_aero.work_vectors
            )
            body_aero.AIC[:, icp, jring] .= velocity_induced
            
            # Subtract 2D induced velocity for VSM
            if icp == jring && model === VSM
                calculate_velocity_induced_bound_2D!(U_2D, body_aero.panels[jring], ep, body_aero.work_vectors)
                body_aero.AIC[:, icp, jring] .-= U_2D
            end
        end
    end
    return nothing
end

"""
    calculate_circulation_distribution_elliptical_wing(body_aero::BodyAerodynamics, gamma_0::Float64=1.0)

Calculate circulation distribution for an elliptical wing.

Returns:
    Vector{Float64}: Circulation distribution along the wing
"""
function calculate_circulation_distribution_elliptical_wing(body_aero::BodyAerodynamics, gamma_0::Float64=1.0)
    length(body_aero.wings) == 1 || throw(ArgumentError("Multiple wings not yet implemented"))
    
    wing_span = body_aero.wings[1].span
    @debug "Wing span: $wing_span"
    
    # Calculate y-coordinates of control points
    y = [panel.control_point[2] for panel in body_aero.panels]
    
    # Calculate elliptical distribution
    gamma_i = gamma_0 * sqrt.(1 .- (2 .* y ./ wing_span).^2)
    
    @debug "Calculated circulation distribution: $gamma_i"
    
    return gamma_i
end

"""
    calculate_stall_angle_list(panels::Vector{Panel};
                             begin_aoa::Float64=9.0,
                             end_aoa::Float64=22.0,
                             step_aoa::Float64=1.0,
                             stall_angle_if_none_detected::Float64=50.0,
                             cl_initial::Float64=-10.0)

Calculate stall angles for each panel.

Returns:
    Vector{Float64}: Stall angles in radians
"""
function calculate_stall_angle_list(panels::Vector{Panel};
                                  begin_aoa::Float64=9.0,
                                  end_aoa::Float64=22.0,
                                  step_aoa::Float64=1.0,
                                  stall_angle_if_none_detected::Float64=50.0,
                                  cl_initial::Float64=-10.0)
    
    aoa_range = deg2rad.(range(begin_aoa, end_aoa, step=step_aoa))
    stall_angles = Float64[]
    
    for panel in panels
        # Default stall angle if none found
        panel_stall = stall_angle_if_none_detected
        
        # Start with minimum cl
        cl_old = cl_initial
        
        # Find stall angle
        for aoa in aoa_range
            cl = calculate_cl(panel, aoa)
            if cl < cl_old
                panel_stall = aoa
                break
            end
            cl_old = cl
        end
        
        push!(stall_angles, panel_stall)
    end
    
    return stall_angles
end

"""
    update_effective_angle_of_attack_if_VSM(body_aero::BodyAerodynamics, gamma::Vector{Float64},
                                          core_radius_fraction::Float64,
                                          z_airf_array::Matrix{Float64},
                                          x_airf_array::Matrix{Float64},
                                          va_array::Matrix{Float64},
                                          va_norm_array::Vector{Float64},
                                          va_unit_array::Matrix{Float64})

Update angle of attack at aerodynamic center for VSM method.

Returns:
    Vector{Float64}: Updated angles of attack
"""
function update_effective_angle_of_attack_if_VSM(body_aero::BodyAerodynamics, 
    gamma::Vector{Float64},
    core_radius_fraction::Float64,
    z_airf_array::Matrix{Float64},
    x_airf_array::Matrix{Float64},
    va_array::Matrix{Float64},
    va_norm_array::Vector{Float64},
    va_unit_array::Matrix{Float64})

    # Calculate AIC matrices at aerodynamic center using LLT method
    calculate_AIC_matrices!(
        body_aero, LLT, core_radius_fraction, va_norm_array, va_unit_array
    )
    AIC_x, AIC_y, AIC_z = @views body_aero.AIC[1, :, :], body_aero.AIC[2, :, :], body_aero.AIC[3, :, :]

    # Calculate induced velocities
    induced_velocity = [
        AIC_x * gamma,
        AIC_y * gamma,
        AIC_z * gamma
    ]
    induced_velocity = hcat(induced_velocity...)
    
    # Calculate relative velocities and angles
    relative_velocity = va_array + induced_velocity
    v_normal = sum(z_airf_array .* relative_velocity, dims=2)
    v_tangential = sum(x_airf_array .* relative_velocity, dims=2)
    alpha_array = atan.(v_normal ./ v_tangential)
    return alpha_array
end

"""
    calculate_results(body_aero::BodyAerodynamics, gamma_new::Vector{Float64}, 
                     density::Float64, aerodynamic_model_type::Model,
                     core_radius_fraction::Float64, mu::Float64,
                     alpha_array::Vector{Float64}, v_a_array::Vector{Float64},
                     chord_array::Vector{Float64}, x_airf_array::Matrix{Float64},
                     y_airf_array::Matrix{Float64}, z_airf_array::Matrix{Float64},
                     va_array::Matrix{Float64}, va_norm_array::Vector{Float64},
                     va_unit_array::Matrix{Float64}, panels::Vector{Panel},
                     is_only_f_and_gamma_output::Bool)

Calculate final aerodynamic results. Reference point is in the kite body (KB) frame.

Returns:
    Dict: Results including forces, coefficients and distributions
"""
function calculate_results(
    body_aero::BodyAerodynamics,
    gamma_new::Vector{Float64},
    reference_point::AbstractVector,
    density::Float64,
    aerodynamic_model_type::Model,
    core_radius_fraction::Float64,
    mu::Float64,
    alpha_array::Vector{Float64},
    v_a_array::Vector{Float64},
    chord_array::Vector{Float64},
    x_airf_array::Matrix{Float64},
    y_airf_array::Matrix{Float64},
    z_airf_array::Matrix{Float64},
    va_array::Matrix{Float64},
    va_norm_array::Vector{Float64},
    va_unit_array::Matrix{Float64},
    panels::Vector{Panel},
    is_only_f_and_gamma_output::Bool,
)

    # Initialize arrays
    n_panels = length(panels)
    cl_array = zeros(n_panels)
    cd_array = zeros(n_panels)
    cm_array = zeros(n_panels)
    panel_width_array = zeros(n_panels)

    # Calculate coefficients for each panel
    for (i, panel) in enumerate(panels)
        cl_array[i] = calculate_cl(panel, alpha_array[i])
        cd_array[i], cm_array[i] = calculate_cd_cm(panel, alpha_array[i])
        panel_width_array[i] = panel.width
    end

    # Calculate forces
    lift = reshape((cl_array .* 0.5 .* density .* v_a_array.^2 .* chord_array), :, 1)
    drag = reshape((cd_array .* 0.5 .* density .* v_a_array.^2 .* chord_array), :, 1)
    moment = reshape((cm_array .* 0.5 .* density .* v_a_array.^2 .* chord_array), :, 1)

    # Calculate alpha corrections based on model type
    alpha_corrected = if aerodynamic_model_type === VSM
        update_effective_angle_of_attack_if_VSM(
            body_aero,
            gamma_new,
            core_radius_fraction,
            z_airf_array,
            x_airf_array,
            va_array,
            va_norm_array,
            va_unit_array
        )
    elseif aerodynamic_model_type === LLT
        alpha_array
    else
        throw(ArgumentError("Unknown aerodynamic model type, should be LLT or VSM"))
    end

    # Verify va is not distributed
    if length(body_aero.va) != 3
        throw(ArgumentError("calculate_results not ready for va_distributed input"))
    end

    # Initialize result arrays
    cl_prescribed_va = Float64[]
    cd_prescribed_va = Float64[]
    cs_prescribed_va = Float64[]
    f_body_3D = zeros(3, n_panels)
    m_body_3D = zeros(3, n_panels)
    area_all_panels = 0.0
    
    # Initialize force sums
    lift_wing_3D_sum = 0.0
    drag_wing_3D_sum = 0.0
    side_wing_3D_sum = 0.0

    # Get wing properties
    spanwise_direction = body_aero.wings[1].spanwise_direction
    va_mag = norm(body_aero.va)
    va = body_aero.va
    va_unit = va / va_mag
    q_inf = 0.5 * density * va_mag^2

    # Main calculation loop
    for (i, panel) in enumerate(panels)
        ### Lift and Drag ###
        # Panel geometry
        panel_area = panel.chord * panel.width
        area_all_panels += panel_area

        # Calculate induced velocity direction
        alpha_corrected_i = alpha_corrected[i]
        induced_va_airfoil = cos(alpha_corrected_i) * panel.x_airf + 
                            sin(alpha_corrected_i) * panel.z_airf
        dir_induced_va_airfoil = induced_va_airfoil / norm(induced_va_airfoil)

        # Calculate lift and drag directions
        dir_lift_induced_va = cross(dir_induced_va_airfoil, panel.y_airf)
        dir_lift_induced_va = dir_lift_induced_va / norm(dir_lift_induced_va)
        dir_drag_induced_va = cross(spanwise_direction, dir_lift_induced_va)
        dir_drag_induced_va = dir_drag_induced_va / norm(dir_drag_induced_va)

        # Calculate force vectors
        lift_induced_va = lift[i] * dir_lift_induced_va
        drag_induced_va = drag[i] * dir_drag_induced_va
        ftotal_induced_va = lift_induced_va + drag_induced_va

        # Calculate forces in prescribed wing frame
        dir_lift_prescribed_va = cross(va, spanwise_direction)
        dir_lift_prescribed_va = dir_lift_prescribed_va / norm(dir_lift_prescribed_va)

        # Calculate force components
        lift_prescribed_va = dot(lift_induced_va, dir_lift_prescribed_va) + 
                           dot(drag_induced_va, dir_lift_prescribed_va)
        drag_prescribed_va = dot(lift_induced_va, va_unit) + 
                           dot(drag_induced_va, va_unit)
        side_prescribed_va = dot(lift_induced_va, spanwise_direction) + 
                           dot(drag_induced_va, spanwise_direction)

        # Body frame forces
        f_body_3D[:,i] .= [
            dot(ftotal_induced_va, [1.0, 0.0, 0.0]),
            dot(ftotal_induced_va, [0.0, 1.0, 0.0]),
            dot(ftotal_induced_va, [0.0, 0.0, 1.0])
        ] .* panel.width

        # Update sums
        lift_wing_3D_sum += lift_prescribed_va * panel.width
        drag_wing_3D_sum += drag_prescribed_va * panel.width  
        side_wing_3D_sum += side_prescribed_va * panel.width

        # Store coefficients
        push!(cl_prescribed_va, lift_prescribed_va / (q_inf * panel.chord))
        push!(cd_prescribed_va, drag_prescribed_va / (q_inf * panel.chord))
        push!(cs_prescribed_va, side_prescribed_va / (q_inf * panel.chord))

        ### Moment ###
        # (1) Panel aerodynamic center in body frame:
        panel_ac_body = panel.aero_center  # 3D [x, y, z]

        # (2) Convert local (2D) pitching moment to a 3D vector in body coords.
        #     Use the axis around which the moment is defined,
        #     which is the y-axis pointing "spanwise"
        moment_axis_body = panel.y_airf

        # Scale by panel width if your 'moment[i]' is 2D moment-per-unit-span:
        M_local_3D = moment[i] * moment_axis_body * panel.width

        # Vector from panel AC to the chosen reference point:
        r_vector = panel_ac_body - reference_point  # e.g. CG, wing root, etc.

        # Cross product to shift the force from panel AC to ref. point:
        M_shift = cross(r_vector, f_body_3D[:,i])

        # Total panel moment about the reference point:
        m_body_3D[:,i] = M_local_3D + M_shift
    end

    if is_only_f_and_gamma_output
        return Dict{String,Any}(
            "F_distribution" => f_body_3D,
            "gamma_distribution" => gamma_new
        )
    end

    # Calculate wing geometry properties
    projected_area = sum(wing -> calculate_projected_area(wing), body_aero.wings)
    wing_span = body_aero.wings[1].span
    aspect_ratio_projected = wing_span^2 / projected_area

    # Calculate geometric angle of attack
    horizontal_direction = [1.0, 0.0, 0.0]
    alpha_geometric = [rad2deg(acos(dot(panel.x_airf, horizontal_direction) /
                     (norm(panel.x_airf) * norm(horizontal_direction))))
                     for panel in panels]

    # Calculate Reynolds number
    max_chord = maximum(panel.chord for panel in panels)
    reynolds_number = density * va_mag * max_chord / mu

    # Create results dictionary
    results = Dict{String,Any}(
        "Fx" => sum(f_body_3D[1,:]),
        "Fy" => sum(f_body_3D[2,:]),
        "Fz" => sum(f_body_3D[3,:]),
        "Mx" => sum(m_body_3D[1,:]),
        "My" => sum(m_body_3D[2,:]),
        "Mz" => sum(m_body_3D[3,:]),
        "lift" => lift_wing_3D_sum,
        "drag" => drag_wing_3D_sum,
        "side" => side_wing_3D_sum,
        "cl" => lift_wing_3D_sum / (q_inf * projected_area),
        "cd" => drag_wing_3D_sum / (q_inf * projected_area),
        "cs" => side_wing_3D_sum / (q_inf * projected_area),
        "cmx" => sum(m_body_3D[1,:]) / (q_inf * projected_area * max_chord),
        "cmy" => sum(m_body_3D[2,:]) / (q_inf * projected_area * max_chord),
        "cmz" => sum(m_body_3D[3,:]) / (q_inf * projected_area * max_chord),
        "cl_distribution" => cl_prescribed_va,
        "cd_distribution" => cd_prescribed_va,
        "cs_distribution" => cs_prescribed_va,
        "F_distribution" => f_body_3D,
        "cfx" => (sum(f_body_3D[1,:]) / (q_inf * projected_area)),
        "cfy" => (sum(f_body_3D[2,:]) / (q_inf * projected_area)),
        "cfz" => (sum(f_body_3D[3,:]) / (q_inf * projected_area)),
        "alpha_at_ac" => alpha_corrected,
        "alpha_uncorrected" => alpha_array,
        "alpha_geometric" => alpha_geometric,
        "gamma_distribution" => gamma_new,
        "area_all_panels" => area_all_panels,
        "projected_area" => projected_area,
        "wing_span" => wing_span,
        "aspect_ratio_projected" => aspect_ratio_projected,
        "Rey" => reynolds_number
    )

    @debug "Results summary:" cl=results["cl"] cd=results["cd"] cs=results["cs"]
    @debug "Forces:" lift=lift_wing_3D_sum drag=drag_wing_3D_sum side=side_wing_3D_sum
    @debug "Areas:" total=area_all_panels projected=projected_area
    @debug "Aspect ratio:" ar=aspect_ratio_projected

    return results
end


"""
    set_va!(body_aero::BodyAerodynamics, va::VelVector, omega=zeros(MVec3))

Set velocity array and update wake filaments.

# Arguments
- body_aero::BodyAerodynamics: The [BodyAerodynamics](@ref) struct to modify
- `va`: Velocity vector of the apparent wind speed           [m/s]
- `omega`: Turn rate vector around x y and z axis            [rad/s]
"""
function set_va!(body_aero::BodyAerodynamics, va::VelVector, omega=zeros(MVec3))
    
    # Calculate va_distribution based on input type
    va_distribution = if all(omega .== 0.0)
        repeat(reshape(va, 1, 3), length(body_aero.panels))
    elseif !all(omega .== 0.0)
        va_dist = zeros(length(body_aero.panels), 3)
        
        for wing in body_aero.wings
            # Get spanwise positions
            spanwise_positions = [panel.control_point for panel in body_aero.panels]
            
            # Calculate velocities for each panel
            for i in 1:wing.n_panels
                omega_va = -omega × spanwise_positions[i]
                va_dist[i, :] .= omega_va .+ va
            end
        end
        va_dist
    end
    
    # Update panel velocities
    for (i, panel) in enumerate(body_aero.panels)
        panel.va .= va_distribution[i,:]
    end
    
    # Update wake elements
    frozen_wake!(body_aero, va_distribution)
    body_aero._va .= va
    return nothing
end

function set_va!(body_aero::BodyAerodynamics, va_distribution::Vector{VelVector}, omega=zeros(MVec3))
    length(va) != length(body_aero.panels) && throw(ArgumentError("Length of va distribution should be equal to number of panels."))
    
    for (i, panel) in enumerate(body_aero.panels)
        panel.va = va_distribution[i]
    end
    
    # Update wake elements
    frozen_wake!(body_aero, va_distribution)
    body_aero._va = va
    return nothing
end