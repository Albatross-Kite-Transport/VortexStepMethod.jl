
"""
    WingAerodynamics

Main structure for calculating aerodynamic properties of wings.
"""
mutable struct WingAerodynamics
    panels::Vector{Panel}
    n_panels::Int
    wings::Vector{AbstractWing}
    _va::Union{Nothing, Vector{Float64}, Tuple{Vector{Float64}, Float64}}
    gamma_distribution::Union{Nothing, Vector{Float64}}
    alpha_uncorrected::Union{Nothing, Vector{Float64}}
    alpha_corrected::Union{Nothing, Vector{Float64}}
    stall_angle_list::Vector{Float64}

    alpha_array::Vector{Float64}
    v_a_array::Vector{Float64}
    work_vectors::NTuple{10,MVec3}
    AIC::Array{Float64, 3}

    function WingAerodynamics(
        wings::Vector{T};
        aerodynamic_center_location::Float64=0.25,
        control_point_location::Float64=0.75
    ) where T <: AbstractWing
        # Initialize panels
        panels = Panel[]
        for wing in wings
            section_list = refine_aerodynamic_mesh(wing)
            n_panels_per_wing = length(section_list) - 1
            
            # Calculate panel properties
            panel_props = calculate_panel_properties(
                section_list,
                n_panels_per_wing,
                aerodynamic_center_location,
                control_point_location
            )
            
            # Create panels
            for i in 1:n_panels_per_wing
                push!(panels, Panel(
                    section_list[i],
                    section_list[i+1],
                    panel_props.aero_centers[i],
                    panel_props.control_points[i],
                    panel_props.bound_points_1[i],
                    panel_props.bound_points_2[i],
                    panel_props.x_airf[i],
                    panel_props.y_airf[i],
                    panel_props.z_airf[i]
                ))
            end
        end
        
        # Initialize rest of the struct
        n_panels = length(panels)
        stall_angle_list = calculate_stall_angle_list(panels)

        alpha_array = zeros(n_panels)
        v_a_array = zeros(n_panels)    
        work_vectors = ntuple(_ -> zeros(MVec3), 10)
        AIC = zeros(3, n_panels, n_panels)

        new(
            panels,
            n_panels,
            wings,
            nothing,  # va
            nothing,  # gamma_distribution
            nothing,  # alpha_uncorrected
            nothing,  # alpha_corrected
            stall_angle_list,
            alpha_array,
            v_a_array,
            work_vectors,
            AIC
        )
    end
end

function Base.getproperty(obj::WingAerodynamics, sym::Symbol)
    if sym === :va
        return getfield(obj, :_va)
    end
    return getfield(obj, sym)
end

function Base.setproperty!(obj::WingAerodynamics, sym::Symbol, val)
    if sym === :va
        set_va!(obj, val)
    else
        setfield!(obj, sym, val)
    end
end

"""
    PanelProperties

Structure to hold calculated panel properties.
"""
struct PanelProperties
    aero_centers::Vector{PosVector}
    control_points::Vector{PosVector}
    bound_points_1::Vector{PosVector}
    bound_points_2::Vector{PosVector}
    x_airf::Vector{Vector{Float64}}
    y_airf::Vector{Vector{Float64}}
    z_airf::Vector{Vector{Float64}}
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
    aero_centers = Vector{Float64}[]
    control_points = Vector{Float64}[]
    bound_points_1 = Vector{Float64}[]
    bound_points_2 = Vector{Float64}[]
    x_airf = Vector{Float64}[]
    y_airf = Vector{Float64}[]
    z_airf = Vector{Float64}[]
    
    # Define coordinates matrix
    coords = zeros(2 * (n_panels + 1), 3)
    @debug "Shape of coordinates: $(size(coords))"
    
    for i in 1:n_panels
        coords[2i-1, :] = section_list[i].LE_point
        coords[2i, :] = section_list[i].TE_point
        coords[2i+1, :] = section_list[i+1].LE_point
        coords[2i+2, :] = section_list[i+1].TE_point
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
        x_airf_vec = cross(VSM_point - LL_point, section["p1"] - section["p2"])
        x_airf_vec = x_airf_vec / norm(x_airf_vec)
        
        y_airf_vec = VSM_point - LL_point
        y_airf_vec = y_airf_vec / norm(y_airf_vec)
        
        z_airf_vec = bound_1 - bound_2
        z_airf_vec = z_airf_vec / norm(z_airf_vec)
        
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
    calculate_AIC_matrices!(wa::WingAerodynamics, model::String, 
                         core_radius_fraction::Float64,
                         va_norm_array::Vector{Float64}, 
                         va_unit_array::Matrix{Float64})

Calculate Aerodynamic Influence Coefficient matrices.

Returns:
    Tuple of (AIC_x, AIC_y, AIC_z) matrices
"""
function calculate_AIC_matrices!(wa::WingAerodynamics, model,
                              core_radius_fraction,
                              va_norm_array, 
                              va_unit_array)
    model in ["VSM", "LLT"] || throw(ArgumentError("Model must be VSM or LLT"))
    # Determine evaluation point based on model
    evaluation_point = model == "VSM" ? :control_point : :aerodynamic_center
    evaluation_point_on_bound = model == "LLT"
    
    # Initialize AIC matrices
    velocity_induced, tempvel, va_unit, U_2D = zeros(MVec3), zeros(MVec3), zeros(MVec3), zeros(MVec3)
    
    # Calculate influence coefficients
    for icp in 1:wa.n_panels
        ep = getproperty(wa.panels[icp], evaluation_point)
        for jring in 1:wa.n_panels
            va_unit .= @views va_unit_array[jring, :]
            filaments = wa.panels[jring].filaments
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
                wa.work_vectors
            )
            wa.AIC[:, icp, jring] .= velocity_induced
            
            # Subtract 2D induced velocity for VSM
            if icp == jring && model == "VSM"
                calculate_velocity_induced_bound_2D!(U_2D, wa.panels[jring], ep, wa.work_vectors)
                wa.AIC[:, icp, jring] .-= U_2D
            end
        end
    end
    return nothing
end

"""
    calculate_circulation_distribution_elliptical_wing(wa::WingAerodynamics, gamma_0::Float64=1.0)

Calculate circulation distribution for an elliptical wing.

Returns:
    Vector{Float64}: Circulation distribution along the wing
"""
function calculate_circulation_distribution_elliptical_wing(wa::WingAerodynamics, gamma_0::Float64=1.0)
    length(wa.wings) == 1 || throw(ArgumentError("Multiple wings not yet implemented"))
    
    wing_span = wa.wings[1].span
    @debug "Wing span: $wing_span"
    
    # Calculate y-coordinates of control points
    y = [panel.control_point[2] for panel in wa.panels]
    
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
    update_effective_angle_of_attack_if_VSM(wa::WingAerodynamics, gamma::Vector{Float64},
                                          core_radius_fraction::Float64,
                                          x_airf_array::Matrix{Float64},
                                          y_airf_array::Matrix{Float64},
                                          va_array::Matrix{Float64},
                                          va_norm_array::Vector{Float64},
                                          va_unit_array::Matrix{Float64})

Update angle of attack at aerodynamic center for VSM method.

Returns:
    Vector{Float64}: Updated angles of attack
"""
function update_effective_angle_of_attack_if_VSM(wa::WingAerodynamics, 
    gamma::Vector{Float64},
    core_radius_fraction::Float64,
    x_airf_array::Matrix{Float64},
    y_airf_array::Matrix{Float64},
    va_array::Matrix{Float64},
    va_norm_array::Vector{Float64},
    va_unit_array::Matrix{Float64})

    # Calculate AIC matrices at aerodynamic center using LLT method
    calculate_AIC_matrices!(
        wa, "LLT", core_radius_fraction, va_norm_array, va_unit_array
    )
    AIC_x, AIC_y, AIC_z = @views wa.AIC[1, :, :], wa.AIC[2, :, :], wa.AIC[3, :, :]

    # Calculate induced velocities
    induced_velocity = [
        AIC_x * gamma,
        AIC_y * gamma,
        AIC_z * gamma
    ]
    induced_velocity = hcat(induced_velocity...)
    
    # Calculate relative velocities and angles
    relative_velocity = va_array + induced_velocity
    v_normal = sum(x_airf_array .* relative_velocity, dims=2)
    v_tangential = sum(y_airf_array .* relative_velocity, dims=2)
    alpha_array = atan.(v_normal ./ v_tangential)

    return alpha_array
end

"""
    calculate_results(wa::WingAerodynamics, gamma_new::Vector{Float64}, 
                     density::Float64, aerodynamic_model_type::String,
                     core_radius_fraction::Float64, mu::Float64,
                     alpha_array::Vector{Float64}, v_a_array::Vector{Float64},
                     chord_array::Vector{Float64}, x_airf_array::Matrix{Float64},
                     y_airf_array::Matrix{Float64}, z_airf_array::Matrix{Float64},
                     va_array::Matrix{Float64}, va_norm_array::Vector{Float64},
                     va_unit_array::Matrix{Float64}, panels::Vector{Panel},
                     is_only_f_and_gamma_output::Bool)

Calculate final aerodynamic results.

Returns:
    Dict: Results including forces, coefficients and distributions
"""
function calculate_results(wa::WingAerodynamics,
    gamma_new::Vector{Float64},
    density::Float64,
    aerodynamic_model_type::String,
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
    is_only_f_and_gamma_output::Bool)

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
    alpha_corrected = if aerodynamic_model_type == "VSM"
        update_effective_angle_of_attack_if_VSM(
            wa,
            gamma_new,
            core_radius_fraction,
            x_airf_array,
            y_airf_array,
            va_array,
            va_norm_array,
            va_unit_array
        )
    elseif aerodynamic_model_type == "LLT"
        alpha_array
    else
        throw(ArgumentError("Unknown aerodynamic model type, should be LLT or VSM"))
    end

    # Verify va is not distributed
    if length(wa.va) != 3
        throw(ArgumentError("calculate_results not ready for va_distributed input"))
    end

    # Initialize result arrays
    cl_prescribed_va = Float64[]
    cd_prescribed_va = Float64[]
    cs_prescribed_va = Float64[]
    f_global_3D = zeros(3, n_panels)
    area_all_panels = 0.0
    
    # Initialize force sums
    lift_wing_3D_sum = 0.0
    drag_wing_3D_sum = 0.0
    side_wing_3D_sum = 0.0

    # Get wing properties
    spanwise_direction = wa.wings[1].spanwise_direction
    va_mag = norm(wa.va)
    va = wa.va
    va_unit = va / va_mag
    q_inf = 0.5 * density * va_mag^2

    # Main calculation loop
    for (i, panel) in enumerate(panels)
        # Panel geometry
        z_airf_span = panel.z_airf
        y_airf_chord = panel.y_airf
        x_airf_normal = panel.x_airf
        panel_area = panel.chord * panel.width
        area_all_panels += panel_area

        # Calculate induced velocity direction
        alpha_corrected_i = alpha_corrected[i]
        induced_va_airfoil = cos(alpha_corrected_i) * y_airf_chord + 
                            sin(alpha_corrected_i) * x_airf_normal
        dir_induced_va_airfoil = induced_va_airfoil / norm(induced_va_airfoil)

        # Calculate lift and drag directions
        dir_lift_induced_va = cross(dir_induced_va_airfoil, z_airf_span)
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

        # Global forces
        f_global_3D[:,i] .= [
            dot(ftotal_induced_va, [1.0, 0.0, 0.0]),
            dot(ftotal_induced_va, [0.0, 1.0, 0.0]),
            dot(ftotal_induced_va, [0.0, 0.0, 1.0])
        ] .* panel.width

        # Update sums
        lift_wing_3D_sum += lift_prescribed_va * panel.width
        drag_wing_3D_sum += drag_prescribed_va * panel.width  
        side_wing_3D_sum += side_prescribed_va * panel.width

        # TODO make this work
        # fx_global_3D_sum += fx_global_3D
        # fy_global_3D_sum += fy_global_3D
        # fz_global_3D_sum += fz_global_3D
        
        # Store coefficients
        push!(cl_prescribed_va, lift_prescribed_va / (q_inf * panel.chord))
        push!(cd_prescribed_va, drag_prescribed_va / (q_inf * panel.chord))
        push!(cs_prescribed_va, side_prescribed_va / (q_inf * panel.chord))

        # TODO translate this
        # fx_global_3D_list.append(fx_global_3D)
        # fy_global_3D_list.append(fy_global_3D)
        # fz_global_3D_list.append(fz_global_3D)
        # f_global_3D_list.append(
        #     np.array([fx_global_3D, fy_global_3D, fz_global_3D])
        # )
    end

    if is_only_f_and_gamma_output
        return Dict{String,Any}(
            "F_distribution" => f_global_3D,
            "gamma_distribution" => gamma_new
        )
    end

    # Calculate wing geometry properties
    projected_area = sum(wing -> calculate_projected_area(wing), wa.wings)
    wing_span = wa.wings[1].span
    aspect_ratio_projected = wing_span^2 / projected_area

    # Calculate geometric angle of attack
    horizontal_direction = [1.0, 0.0, 0.0]
    alpha_geometric = [rad2deg(acos(dot(panel.y_airf, horizontal_direction) /
                     (norm(panel.y_airf) * norm(horizontal_direction))))
                     for panel in panels]

    # Calculate Reynolds number
    max_chord = maximum(panel.chord for panel in panels)
    reynolds_number = density * va_mag * max_chord / mu

    # Create results dictionary
    results = Dict{String,Any}(
        "Fx" => sum(f_global_3D[1,:]),
        "Fy" => sum(f_global_3D[1,:]),
        "Fz" => sum(f_global_3D[1,:]),
        "lift" => lift_wing_3D_sum,
        "drag" => drag_wing_3D_sum,
        "side" => side_wing_3D_sum,
        "cl" => lift_wing_3D_sum / (q_inf * projected_area),
        "cd" => drag_wing_3D_sum / (q_inf * projected_area),
        "cs" => side_wing_3D_sum / (q_inf * projected_area),
        "cl_distribution" => cl_prescribed_va,
        "cd_distribution" => cd_prescribed_va,
        "cs_distribution" => cs_prescribed_va,
        "F_distribution" => f_global_3D,
        "cfx" => (sum(f_global_3D[1,:]) / (q_inf * projected_area)),
        "cfy" => (sum(f_global_3D[2,:]) / (q_inf * projected_area)),
        "cfz" => (sum(f_global_3D[3,:]) / (q_inf * projected_area)),
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
    set_va!(wa::WingAerodynamics, va::Union{Vector{Float64}, Tuple{Vector{Float64}, Float64}})

Set velocity array and update wake filaments.

# Arguments
- `va`: Either a velocity vector or tuple of (velocity vector, yaw_rate)
"""
function set_va!(wa::WingAerodynamics, va)
    # Add length check for va_vec
    if va isa Vector{Float64} && length(va) != 3 && length(va) != wa.n_panels
        throw(ArgumentError("va must be length 3 or match number of panels"))
    end
    # Handle input types
    va_vec, yaw_rate = if va isa Tuple && length(va) == 2
        va
    else
        (va, 0.0)
    end
    
    # Validate input
    va_vec = convert(Vector{Float64}, va_vec)
    
    # Calculate va_distribution based on input type
    va_distribution = if length(va_vec) == 3 && yaw_rate == 0.0
        repeat(reshape(va_vec, 1, 3), wa.n_panels)
    elseif length(va_vec) == wa.n_panels
        va_vec
    elseif yaw_rate != 0.0 && length(va_vec) == 3
        va_dist = Vector{Float64}[]
        
        for wing in wa.wings
            # Get spanwise positions
            spanwise_positions = [panel.control_point[2] for panel in wa.panels]
            
            # Calculate velocities for each panel
            for i in 1:wing.n_panels
                yaw_rate_apparent_velocity = [-yaw_rate * spanwise_positions[i], 0.0, 0.0]
                push!(va_dist, yaw_rate_apparent_velocity + va_vec)
            end
        end
        reduce(vcat, va_dist)
    else
        throw(ArgumentError("Invalid va distribution: length(va)=$(length(va_vec)) ≠ n_panels=$(wa.n_panels)"))
    end
    
    # Update panel velocities
    for (i, panel) in enumerate(wa.panels)
        panel.va = va_distribution[i,:]
    end
    
    # Update wake elements
    wa.panels = frozen_wake(va_distribution, wa.panels)
    wa._va = va_vec
end
