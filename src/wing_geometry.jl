
"""
    Section{T}

Represents a wing section with leading edge, trailing edge, and aerodynamic properties.

# Fields
- `LE_point::MVec3`: Leading edge point coordinates
- `TE_point::MVec3`: Trailing edge point coordinates
- `aero_input`::T: Aerodynamic input, one of three possible types
"""
struct Section{T}
    LE_point::MVec3
    TE_point::MVec3
    aero_input::T
    
    function Section(
            LE_point::PosVector, 
            TE_point::PosVector, 
            aero_input::T
            ) where T
        new{T}(LE_point, TE_point, aero_input)
    end
end

"""
    update_pos!(wing::Wing, LE_point::PosVector, TE_point::PosVector)

Update the section leading edge and trailing edge positions.
"""
function update_pos!(section::Section, LE_point::AbstractVector, TE_point::AbstractVector)
    section.LE_point .= LE_point
    section.TE_point .= TE_point
    return nothing
end

"""
    Wing

Represents a wing composed of multiple sections with aerodynamic properties.

# Fields
- `n_panels::Int64`: Number of panels in aerodynamic mesh
- `spanwise_panel_distribution`::PanelDistribution: [PanelDistribution](@ref)
- `spanwise_direction::Vector{Float64}`: Wing span direction vector
- `sections::Vector{Section}`: List of wing sections, see: [Section](@ref)

"""
mutable struct Wing <: AbstractWing
    n_panels::Int64
    spanwise_panel_distribution::PanelDistribution
    spanwise_direction::PosVector
    sections::Vector{Section}
    
    function Wing(n_panels::Int;
                 spanwise_panel_distribution::PanelDistribution=LINEAR,
                 spanwise_direction::PosVector=MVec3([0.0, 1.0, 0.0]))
        new(n_panels, 
            spanwise_panel_distribution, 
            spanwise_direction, 
            Section[])
    end
end

"""
    add_section!(wing::Wing, LE_point::PosVector, TE_point::PosVector, aero_input)

Add a new section to the wing.

# Arguments:
- wing::Wing: The [Wing](@ref) to which a section shall be added
- LE_point::PosVector: [PosVector](@ref) of the point on the side of the leading edge
- TE_point::PosVector: [PosVector](@ref) of the point on the side of the trailing edge
- aero_input: Can be:
  - :inviscid
  - :`lei_airfoil_breukels`
  - (:`polar_data`, (`alpha_range`, `cl_vector`, `cd_vector`, `cm_vector`))
  - (:`polar_data`, (`alpha_range`, `beta_range`, `cl_matrix`, `cd_matrix`, `cm_matrix`))
"""
function add_section!(wing::Wing, LE_point::Vector{Float64}, 
                     TE_point::Vector{Float64}, aero_input)
    push!(wing.sections, Section(LE_point, TE_point, aero_input))
    return nothing
end

"""
    flip_created_coord_in_pairs_if_needed!(coord::Matrix{Float64})

Ensure coordinates are ordered from positive to negative along y-axis.
"""
function flip_created_coord_in_pairs_if_needed!(coord::Matrix{Float64})
    n_pairs = size(coord, 1) ÷ 2
    reshaped = reshape(coord, (n_pairs, 2, size(coord, 2)))
    
    # Check y-values order
    y_values = reshaped[:, 1, 2]  # y-coordinates of leading edge points
    if !all(y_values[1:end-1] .>= y_values[2:end])
        reshaped = reverse(reshaped, dims=1)
    end
    
    return reshape(reshaped, (size(coord, 1), size(coord, 2)))
end


"""
    refine_aerodynamic_mesh(wing::AbstractWing)

Refine the aerodynamic mesh of the wing based on spanwise panel distribution.

Returns:
    Vector{Section}: List of refined sections
"""
function refine_aerodynamic_mesh(wing::AbstractWing)
    # Sort sections from left to right
    sort!(wing.sections, by=s -> s.LE_point[2], rev=true)
    
    # Calculate number of sections needed
    n_sections = wing.n_panels + 1
    @debug "n_panels: $(wing.n_panels)"
    @debug "n_sections: $n_sections"
    
    # Extract geometry data
    n_current = length(wing.sections)
    LE = zeros(Float64, n_current, 3)
    TE = zeros(Float64, n_current, 3)
    aero_input = Vector{typeof(wing.sections[1].aero_input)}()
    
    for (i, section) in enumerate(wing.sections)
        LE[i,:] = section.LE_point
        TE[i,:] = section.TE_point
        push!(aero_input, section.aero_input)
    end
    
    # Validate input
    if size(LE,1) != size(TE,1) || size(LE,1) != length(aero_input)
        throw(ArgumentError("LE, TE, and aero_input must have the same length"))
    end
    
    # Handle special cases
    if wing.spanwise_panel_distribution == UNCHANGED || length(wing.sections) == n_sections
        return wing.sections
    end

    @info "Refining aerodynamic mesh from $(length(wing.sections)) sections to $n_sections sections."
    
    # Handle two-section case
    if n_sections == 2
        return [
            Section(LE[1,:], TE[1,:], aero_input[1]),
            Section(LE[end,:], TE[end,:], aero_input[end])
        ]
    end
    
    # Handle different distribution types
    if wing.spanwise_panel_distribution == SPLIT_PROVIDED
        return refine_mesh_by_splitting_provided_sections(wing)
    elseif wing.spanwise_panel_distribution in (LINEAR, COSINE, COSINE_VAN_GARREL)
        return refine_mesh_for_linear_cosine_distribution(
            wing.spanwise_panel_distribution,
            n_sections,
            LE,
            TE,
            aero_input
        )
    else
        throw(ArgumentError("Unsupported spanwise panel distribution: $(wing.spanwise_panel_distribution)"))
    end
end

"""
    calculate_new_aero_input(aero_input, 
                            section_index::Int,
                            left_weight::Float64,
                            right_weight::Float64)

Interpolate aerodynamic input between two sections.
"""
function calculate_new_aero_input(aero_input, 
                                section_index::Int,
                                left_weight::Float64,
                                right_weight::Float64)
    
    model_type = isa(aero_input[section_index], Symbol) ? aero_input[section_index] : aero_input[section_index][1]
    model_type_2 = isa(aero_input[section_index + 1], Symbol) ? aero_input[section_index + 1] : aero_input[section_index + 1][1]
    if !(model_type === model_type_2)
        throw(ArgumentError("Different aero models over the span are not supported"))
    end
    
    if model_type === :inviscid
        return :inviscid
        
    elseif model_type === :polar_data
        polar_left = aero_input[section_index][2]
        polar_right = aero_input[section_index + 1][2]
        
        # Unpack polar data
        if length(polar_left) == 4
            alpha_left, CL_left, CD_left, CM_left = polar_left
            alpha_right, CL_right, CD_right, CM_right = polar_right

            # Create common alpha array
            !all(isapprox.(alpha_left, alpha_right)) && @error "Make sure you use the same alpha range for all your interpolations."
            !isa(CL_right, AbstractVector) && @error "Provide polar data in the correct format: (alpha, cl, cd, cm)"
            
            # Weighted interpolation
            CL_data = CL_left .* left_weight .+ CL_right .* right_weight
            CD_data = CD_left .* left_weight .+ CD_right .* right_weight
            CM_data = CM_left .* left_weight .+ CM_right .* right_weight
            
            return (:polar_data, (alpha_left, CL_data, CD_data, CM_data))
            
        elseif length(polar_left) == 5
            alpha_left, beta_left, CL_left, CD_left, CM_left = polar_left
            alpha_right, beta_right, CL_right, CD_right, CM_right = polar_right
            
            # Create common alpha array
            !all(isapprox.(alpha_left, alpha_right)) && @error "Make sure you use the same alpha range for all your interpolations."
            !all(isapprox.(beta_left, beta_right)) && @error "Make sure you use the same alpha range for all your interpolations."
            !isa(CL_right, AbstractMatrix) && @error "Provide polar data in the correct format: (alpha, beta, cl, cd, cm)"

            # Weighted interpolation
            CL_data = CL_left .* left_weight .+ CL_right .* right_weight
            CD_data = CD_left .* left_weight .+ CD_right .* right_weight
            CM_data = CM_left .* left_weight .+ CM_right .* right_weight
            
            return (:polar_data, (alpha_left, beta_left, CL_data, CD_data, CM_data))
        end

    elseif model_type === :lei_airfoil_breukels
        tube_diameter_left = aero_input[section_index][2][1]
        tube_diameter_right = aero_input[section_index + 1][2][1]
        tube_diameter_i = tube_diameter_left * left_weight + tube_diameter_right * right_weight
        
        chamber_height_left = aero_input[section_index][2][2]
        chamber_height_right = aero_input[section_index + 1][2][2]
        chamber_height_i = chamber_height_left * left_weight + chamber_height_right * right_weight
        
        @debug "Interpolation weights" left_weight right_weight
        @debug "Interpolated parameters" tube_diameter_i chamber_height_i
        
        return (:lei_airfoil_breukels, [tube_diameter_i, chamber_height_i])
    else
        throw(ArgumentError("Unsupported aero model: $(model_type)"))
    end
end

"""
    refine_mesh_for_linear_cosine_distribution(
        spanwise_panel_distribution::PanelDistribution,
        n_sections::Int,
        LE::Matrix{Float64},
        TE::Matrix{Float64},
        aero_input)

Refine wing mesh using linear or cosine spacing.

# Arguments
- `spanwise_panel_distribution`: [PanelDistribution](@ref)
- `n_sections`: Number of sections to generate
- `LE`: Matrix of leading edge points
- `TE`: Matrix of trailing edge points
- `aero_input`: Vector of aerodynamic inputs for each section

Returns:
    Vector{Section}: List of refined sections
"""
function refine_mesh_for_linear_cosine_distribution(
    spanwise_panel_distribution::PanelDistribution,
    n_sections::Int,
    LE,
    TE,
    aero_input)

    # 1. Compute quarter chord line
    quarter_chord = LE .+ 0.25 .* (TE .- LE)

    # Calculate segment lengths
    qc_lengths = [norm(quarter_chord[i+1,:] - quarter_chord[i,:]) for i in 1:size(quarter_chord,1)-1]
    qc_total_length = sum(qc_lengths)
    qc_cum_length = vcat(0, cumsum(qc_lengths))

    # 2. Define target lengths
    target_lengths = if spanwise_panel_distribution == LINEAR
        range(0, qc_total_length, n_sections)
    elseif spanwise_panel_distribution in (COSINE, COSINE_VAN_GARREL)
        theta = range(0, π, n_sections)
        qc_total_length .* (1 .- cos.(theta)) ./ 2
    else
        throw(ArgumentError("Unsupported distribution: $spanwise_panel_distribution"))
    end

    # Initialize arrays
    new_quarter_chord = zeros(Float64, n_sections, 3)
    new_LE = zeros(Float64, n_sections, 3)
    new_TE = zeros(Float64, n_sections, 3)
    new_sections = Section[]

    # 3. Calculate new points and interpolate
    for i in 1:n_sections
        target_length = target_lengths[i]

        # Find segment index
        section_index = searchsortedlast(qc_cum_length, target_length) 
        section_index = clamp(section_index, 1, length(qc_cum_length)-1)

        # 4. Calculate weights
        segment_start = qc_cum_length[section_index]
        segment_end = qc_cum_length[section_index+1]
        t = (target_length - segment_start) / (segment_end - segment_start)
        left_weight = 1 - t
        right_weight = t

        # 5. Calculate quarter chord point
        new_quarter_chord[i,:] = quarter_chord[section_index,:] + 
                                t .* (quarter_chord[section_index+1,:] - quarter_chord[section_index,:])

        # 6. Calculate chord vectors
        left_chord = TE[section_index,:] - LE[section_index,:]
        right_chord = TE[section_index+1,:] - LE[section_index+1,:]

        # Normalize chord vectors
        left_chord_norm = left_chord ./ max(norm(left_chord), 1e-12)
        right_chord_norm = right_chord ./ max(norm(right_chord), 1e-12)

        # Interpolate direction
        avg_direction = left_weight .* left_chord_norm .+ right_weight .* right_chord_norm
        avg_direction = avg_direction ./ max(norm(avg_direction), 1e-12)

        # Interpolate length
        left_length = norm(left_chord)
        right_length = norm(right_chord)
        avg_length = left_weight * left_length + right_weight * right_length

        # Final chord vector
        avg_chord = avg_direction .* avg_length

        # Calculate LE and TE points
        new_LE[i,:] = new_quarter_chord[i,:] .- 0.25 .* avg_chord
        new_TE[i,:] = new_quarter_chord[i,:] .+ 0.75 .* avg_chord

        # Interpolate aero properties
        new_aero = calculate_new_aero_input(aero_input, section_index, left_weight, right_weight)

        # Create new section
        push!(new_sections, Section(new_LE[i,:], new_TE[i,:], new_aero))
    end

    # Apply van Garrel distribution if requested
    if spanwise_panel_distribution === :cosine_van_Garrel
        new_sections = calculate_cosine_van_Garrel(new_sections)
    end

    return new_sections
end

"""
    calculate_cosine_van_Garrel(sections::Vector{Section})

Calculate van Garrel cosine distribution of sections.
Reference: http://dx.doi.org/10.13140/RG.2.1.2773.8000

Returns:
    Vector{Section}: Sections with van Garrel distribution
"""
function calculate_cosine_van_Garrel(sections::Vector{Section})
    n = length(sections)
    
    # Calculate chords and quarter chords
    chords = [section.TE_point - section.LE_point for section in sections]
    quarter_chords = [section.LE_point + 0.25 * chord for (section, chord) in zip(sections, chords)]
    
    # Calculate widths
    widths = [norm(quarter_chords[i+1] - quarter_chords[i]) for i in 1:n-1]
    
    # Calculate correction factors
    eta_cp = zeros(n-1)
    
    # First panel
    eta_cp[1] = widths[1] / (widths[1] + widths[2])
    
    # Internal panels
    for j in 2:n-2
        eta_cp[j] = 0.25 * (
            widths[j-1] / (widths[j-1] + widths[j]) +
            widths[j] / (widths[j] + widths[j+1]) + 1
        )
    end
    
    # Last panel
    eta_cp[end] = widths[end-1] / (widths[end-1] + widths[end])
    
    @debug "Correction factors" eta_cp
    
    # Calculate control points
    control_points = [
        quarter_chords[i] + eta * (quarter_chords[i+1] - quarter_chords[i])
        for (i, eta) in enumerate(eta_cp)
    ]
    
    # Generate new sections
    new_sections = [
        Section(
            control_point - 0.25 * chord,  # new LE
            control_point + 0.75 * chord,  # new TE
            sections[i].aero_input         # keep original aero input
        )
        for (i, (control_point, chord)) in enumerate(zip(control_points, chords))
    ]
    
    return new_sections
end


"""
    refine_mesh_by_splitting_provided_sections(wing::AbstractWing)

Refine mesh by splitting provided sections into desired number of panels.

Returns:
    Vector{Section}: Refined sections
"""
function refine_mesh_by_splitting_provided_sections(wing::AbstractWing)
    n_sections_provided = length(wing.sections)
    n_panels_provided = n_sections_provided - 1
    n_panels_desired = wing.n_panels
    
    @debug "Panel counts" n_panels_provided n_panels_desired n_sections_provided
    
    # Check if refinement is needed
    if n_panels_provided == n_panels_desired
        return wing.sections
    end
    
    # Validate panel count relationship
    if n_panels_desired % n_panels_provided != 0
        throw(ArgumentError(
            "Desired panels ($n_panels_desired) must be multiple of provided panels " *
            "($n_panels_provided). Choose: $(n_panels_provided*2), $(n_panels_provided*3), ..."
        ))
    end
    
    # Calculate distribution
    n_new_sections = wing.n_panels + 1 - n_sections_provided
    n_section_pairs = n_sections_provided - 1
    new_sections_per_pair, remaining = divrem(n_new_sections, n_section_pairs)
    
    # Initialize results
    new_sections = Section[]
    
    # Extract geometry data
    LE = [section.LE_point for section in wing.sections]
    TE = [section.TE_point for section in wing.sections]
    aero_input = [section.aero_input for section in wing.sections]
    
    # Process each section pair
    for left_section_index in 1:n_section_pairs
        # Add left section of pair
        push!(new_sections, wing.sections[left_section_index])
        
        # Calculate new sections for this pair
        num_new_sections = new_sections_per_pair + (left_section_index <= remaining ? 1 : 0)
        
        if num_new_sections > 0
            # Prepare pair data
            LE_pair = hcat(LE[left_section_index], LE[left_section_index + 1])'
            TE_pair = hcat(TE[left_section_index], TE[left_section_index + 1])'
            aero_pair = [
                aero_input[left_section_index],
                aero_input[left_section_index + 1]
            ]
            
            # Generate sections for this pair
            new_splitted_sections = refine_mesh_for_linear_cosine_distribution(
                LINEAR,
                num_new_sections + 2,  # +2 for endpoints
                LE_pair,
                TE_pair,
                aero_pair
            )
            
            # Add new sections (excluding endpoints)
            append!(new_sections, new_splitted_sections[2:end-1])
        end
    end
    
    # Add final section
    push!(new_sections, wing.sections[end])
    
    # Validate result
    if length(new_sections) != wing.n_panels + 1
        @warn "Number of panels ($(length(new_sections)-1)) differs from desired ($(wing.n_panels))"
    end
    
    return new_sections
end

"""
    calculate_span(wing::AbstractWing)

Calculate wing span along spanwise direction.

Returns:
    Float64: Wing span
"""
function calculate_span(wing::AbstractWing)
    # Normalize spanwise direction
    vector_axis = wing.spanwise_direction ./ norm(wing.spanwise_direction)
    
    # Get all points
    all_points = reduce(vcat, [[section.LE_point, section.TE_point] 
                              for section in wing.sections])
    
    # Project points and calculate span
    projections = [dot(point, vector_axis) for point in all_points]
    return maximum(projections) - minimum(projections)
end

"""
    calculate_projected_area(wing::AbstractWing, z_plane_vector::Vector{Float64}=[0.0, 0.0, 1.0])

Calculate projected wing area onto plane defined by normal vector.

Returns:
    Float64: Projected area
"""
function calculate_projected_area(wing::AbstractWing, 
                                z_plane_vector::Vector{Float64}=[0.0, 0.0, 1.0])
    # Normalize plane normal vector
    z_plane_vector = z_plane_vector ./ norm(z_plane_vector)
    
    # Project point onto plane
    function project_onto_plane(point::PosVector, normal::Vector{Float64})
        return point .- dot(point, normal) .* normal
    end
    
    # Calculate area by summing trapezoid areas
    projected_area = 0.0
    for i in 1:(length(wing.sections)-1)
        # Get section points
        LE_current = wing.sections[i].LE_point
        TE_current = wing.sections[i].TE_point
        LE_next = wing.sections[i+1].LE_point
        TE_next = wing.sections[i+1].TE_point
        
        # Project points
        LE_current_proj = project_onto_plane(LE_current, z_plane_vector)
        TE_current_proj = project_onto_plane(TE_current, z_plane_vector)
        LE_next_proj = project_onto_plane(LE_next, z_plane_vector)
        TE_next_proj = project_onto_plane(TE_next, z_plane_vector)
        
        # Calculate projected dimensions
        chord_current = norm(TE_current_proj - LE_current_proj)
        chord_next = norm(TE_next_proj - LE_next_proj)
        span = norm(LE_next_proj - LE_current_proj)
        
        # Add trapezoid area
        projected_area += 0.5 * (chord_current + chord_next) * span
    end
    
    return projected_area
end

# Add span property to Wing struct
Base.propertynames(w::AbstractWing) = (fieldnames(typeof(w))..., :span)
function Base.getproperty(w::AbstractWing, s::Symbol)
    if s === :span
        return calculate_span(w)
    else
        return getfield(w, s)
    end
end