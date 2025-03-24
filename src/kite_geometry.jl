"""
    read_faces(filename)

Read vertices and faces from an OBJ file.

# Arguments
- `filename::String`: Path to .obj file

# Returns
- Tuple of (vertices, faces) where:
  - vertices: Vector of 3D coordinates [x,y,z]
  - faces: Vector of triangle vertex indices
"""
function read_faces(filename)
    vertices = []
    faces = []
    
    open(filename) do file
        for line in eachline(file)
            if startswith(line, "v ") && !startswith(line, "vt") && !startswith(line, "vn")
                parts = split(line)
                x = parse(Float64, parts[2])
                y = parse(Float64, parts[3])
                z = parse(Float64, parts[4])
                push!(vertices, [x, y, z])
            elseif startswith(line, "f ")
                parts = split(line)
                # Handle both f v1 v2 v3 and f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3 formats
                indices = map(p -> parse(Int64, split(p, '/')[1]), parts[2:4])
                push!(faces, indices)
            end
        end
    end
    return vertices, faces
end

"""
    find_circle_center_and_radius(vertices)

Find the center and radius of the kite's curvature circle.

# Arguments
- `vertices`: Vector of 3D point coordinates

# Returns
- Tuple of (z_center, radius, gamma_tip) where:
  - z_center: Z-coordinate of circle center
  - radius: Circle radius
  - gamma_tip: Angle of the kite tip from z-axis
"""
function find_circle_center_and_radius(vertices)
    r = zeros(2)
    v_min = zeros(3)
    v_tip = zeros(3)
    v_min .= Inf

    # find the vertex with smallest x in the middle of the kite
    for v in vertices
        if abs(v[2]) ≤ 0.1
            if v[1] < v_min[1]
                v_min .= v
            end
        end
    end

    # Find vertex furthest in -y, -z direction
    max_score = -Inf
    v_tip .= 0.0
    for v in vertices
        # Score each vertex based on -y and -z components
        # lower y and lower z gives higher score
        score = -v[2] - v[3]  # - y - z
        if score > max_score
            max_score = score
            v_tip .= v
        end
    end

    function r_diff!(du, u, p)
        z = u[1]
        r .= Inf
        r[1] = sqrt(v_min[2]^2 + (v_min[3] - z)^2)
        r[2] = sqrt(v_tip[2]^2 + (v_tip[3] - z)^2)
        du[1] = r[1] - r[2]
        return nothing
    end

    prob = NonlinearProblem(r_diff!, [v_min[3]-0.1], nothing)
    result = NonlinearSolve.solve(prob, NewtonRaphson(; autodiff=AutoFiniteDiff(; relstep = 1e-3, absstep = 1e-3)); abstol = 1e-2)
    r_diff!(zeros(1), result, nothing)
    z = result[1]

    gamma_tip = atan(-v_tip[2], (v_tip[3] - z))
    @assert gamma_tip > 0.0

    return z, r[1], gamma_tip
end

"""
    create_interpolations(vertices, circle_center_z, radius, gamma_tip)

Create interpolation functions for leading/trailing edges and area.

# Arguments
- `vertices`: Vector of 3D point coordinates
- `circle_center_z`: Z-coordinate of circle center
- `radius`: Circle radius
- `gamma_tip`: Maximum angular extent

# Returns
- Tuple of (le_interp, te_interp, area_interp) interpolation functions
- Where le_interp and te_interp are tuples themselves, containing the x, y and z interpolations
"""
function create_interpolations(vertices, circle_center_z, radius, gamma_tip, R; interp_steps=40)
    gamma_range = range(-gamma_tip+1e-6, gamma_tip-1e-6, interp_steps)
    stepsize = gamma_range.step.hi
    vz_centered = [v[3] - circle_center_z for v in vertices]
    
    te_gammas = zeros(length(gamma_range))
    le_gammas = zeros(length(gamma_range))
    trailing_edges = zeros(3, length(gamma_range))
    leading_edges = zeros(3, length(gamma_range))
    areas  = zeros(length(gamma_range))
    
    for (j, gamma) in enumerate(gamma_range)
        trailing_edges[1, j] = -Inf
        leading_edges[1, j] = Inf
        for (i, v) in enumerate(vertices)
            # Rotate y coordinate to check box containment
            # rotated_y = v[2] * cos(-gamma) - vz_centered[i] * sin(-gamma)
            gamma_v = atan(-v[2], vz_centered[i])
            if gamma < 0 && gamma - stepsize ≤ gamma_v ≤ gamma
                if v[1] > trailing_edges[1, j]
                    trailing_edges[:, j] .= v
                    te_gammas[j] = gamma_v
                end
                if v[1] < leading_edges[1, j]
                    leading_edges[:, j] .= v
                    le_gammas[j] = gamma_v
                end
            elseif gamma > 0 && gamma ≤ gamma_v ≤ gamma + stepsize
                if v[1] > trailing_edges[1, j]
                    trailing_edges[:, j] .= v
                    te_gammas[j] = gamma_v
                end
                if v[1] < leading_edges[1, j]
                    leading_edges[:, j] .= v
                    le_gammas[j] = gamma_v
                end
            end
        end
        area = norm(leading_edges[:, j] - trailing_edges[:, j]) * stepsize * radius
        last_area = j > 1 ? areas[j-1] : 0.0
        areas[j] = last_area + area
    end

    for j in eachindex(gamma_range)
        leading_edges[:, j] .= R * leading_edges[:, j]
        trailing_edges[:, j] .= R * trailing_edges[:, j]
    end

    le_interp = ntuple(i -> linear_interpolation(te_gammas, leading_edges[i, :],
                                           extrapolation_bc=Line()), 3)
    te_interp = ntuple(i -> linear_interpolation(le_gammas, trailing_edges[i, :],
                                           extrapolation_bc=Line()), 3)
    area_interp = linear_interpolation(gamma_range, areas, extrapolation_bc=Line())
    
    return (le_interp, te_interp, area_interp)
end

"""
    center_to_com!(vertices, faces)

Calculate center of mass of a mesh and translate vertices so that COM is at origin.

# Arguments
- `vertices`: Vector of 3D point coordinates
- `faces`: Vector of vertex indices for each face (can be triangular or non-triangular)

# Returns
- Vector representing the original center of mass before translation

# Notes
- Non-triangular faces are automatically triangulated into triangles
- Assumes uniform surface density
"""
function center_to_com!(vertices, faces)
    area_total = 0.0
    com = zeros(3)
    
    for face in faces
        if length(face) == 3
            # Triangle case
            v1 = vertices[face[1]]
            v2 = vertices[face[2]]
            v3 = vertices[face[3]]
            
            # Calculate triangle area and centroid
            normal = cross(v2 - v1, v3 - v1)
            area = norm(normal) / 2
            centroid = (v1 + v2 + v3) / 3
            
            area_total += area
            com += area * centroid
        else
            throw(ArgumentError("Triangulate faces in a CAD program first"))
        end
    end
    
    com = com / area_total
    !(abs(com[2]) < 0.01) && throw(ArgumentError("Center of mass $com of .obj file has to lie on the xz-plane."))
    @info "Centering vertices of .obj file to the center of mass: $com"
    com[2] = 0.0
    for v in vertices
        v .-= com
    end
    return com
end

"""
    calculate_inertia_tensor(vertices, faces, mass, com)

Calculate the inertia tensor for a triangulated surface mesh, assuming a thin shell with uniform 
surface density.

# Arguments
- `vertices`: Vector of 3D point coordinates representing mesh vertices
- `faces`: Vector of triangle indices, each defining a face of the mesh
- `mass`: Total mass of the shell in kg
- `com`: Center of mass coordinates [x,y,z]

# Method
Uses the thin shell approximation where:
1. Mass is distributed uniformly over the surface area
2. Each triangle contributes to the inertia based on its area and position
3. For each triangle vertex p, contribution to diagonal terms is: area * (sum(p²) - p_i²)
4. For off-diagonal terms: area * (-`p_i` * `p_j`)
5. Final tensor is scaled by mass/(3*total_area) to get correct units

# Returns
- 3×3 matrix representing the inertia tensor in kg⋅m²
"""
function calculate_inertia_tensor(vertices, faces, mass, com)
    # Initialize inertia tensor
    I = zeros(3, 3)
    total_area = 0.0
    
    for face in faces
        v1 = vertices[face[1]] .- com
        v2 = vertices[face[2]] .- com
        v3 = vertices[face[3]] .- com
        
        # Calculate triangle area
        normal = cross(v2 - v1, v3 - v1)
        area = norm(normal) / 2
        total_area += area
        
        # Calculate contribution to inertia tensor
        for i in 1:3
            for j in 1:3
                # Vertices relative to center of mass
                points = [v1, v2, v3]
                
                # Calculate contribution to inertia tensor
                for p in points
                    if i == j
                        # Diagonal terms
                        I[i,i] += area * (sum(p.^2) - p[i]^2)
                    else
                        # Off-diagonal terms
                        I[i,j] -= area * (p[i] * p[j])
                    end
                end
            end
        end
    end
    
    # Scale by mass/total_area to get actual inertia tensor
    return (mass / total_area) * I / 3
end

function calc_inertia_y_rotation(I_b_tensor)
    # Function for nonlinear solver - off-diagonal element should be zero
    function eq!(F, theta, _)
        # Rotation matrix around y-axis
        R_y = [
            cos(theta[1])  0  sin(theta[1]);
            0              1  0;
            -sin(theta[1]) 0  cos(theta[1])
        ]
        # Transform inertia tensor
        I_rotated = R_y * I_b_tensor * R_y'
        # We want the off-diagonal xz elements to be zero
        F[1] = I_rotated[1,3] 
    end
    
    theta0 = [0.0]
    prob = NonlinearProblem(eq!, theta0, nothing)
    sol = NonlinearSolve.solve(prob, NewtonRaphson())
    theta_opt = sol.u[1]
    
    R_b_p = [
        cos(theta_opt)  0  sin(theta_opt);
        0               1  0;
        -sin(theta_opt) 0  cos(theta_opt)
    ]
    # Calculate diagonalized inertia tensor
    I_diag = R_b_p * I_b_tensor * R_b_p'
    @assert isapprox(I_diag[1,3], 0.0, atol=1e-5)
    return I_diag, R_b_p
end

"""
    interpolate_matrix_nans!(matrix::Matrix{Float64})

Replace NaN values in a matrix by interpolating from nearest non-NaN neighbors.
Uses an expanding search radius until valid neighbors are found.

# Arguments
- `matrix`: Matrix containing NaN values to be interpolated
"""
function interpolate_matrix_nans!(matrix::Matrix{Float64})
    rows, cols = size(matrix)
    nans_found = 0
    while any(isnan, matrix)
        for i in 1:rows, j in 1:cols
            if isnan(matrix[i,j])
                # Search in expanding radius until we find valid neighbors
                radius = 1
                values = Float64[]
                weights = Float64[]
                
                while isempty(values) && radius < max(rows, cols)
                    # Check all points at current Manhattan distance
                    for di in -radius:radius, dj in -radius:radius
                        if abs(di) + abs(dj) == radius  # Points exactly at distance 'radius'
                            ni, nj = i + di, j + dj
                            if 1 ≤ ni ≤ rows && 1 ≤ nj ≤ cols && !isnan(matrix[ni,nj])
                                # Weight by inverse distance
                                dist = sqrt(di^2 + dj^2)
                                push!(values, matrix[ni,nj])
                                push!(weights, 1/dist)
                            end
                        end
                    end
                    radius += 1
                end
                
                if !isempty(values)
                    # Calculate weighted average of found values
                    matrix[i,j] = sum(values .* weights) / sum(weights)
                    nans_found += 1
                else
                    throw(ArgumentError("Could not remove NaN"))
                end
            end
        end
    end
    @info "Removed $nans_found NaNs from the matrix."
    return matrix
end


"""
    RamAirWing

Represents a curved wing that inherits from Wing with additional geometric properties.

# Fields
- All fields from Wing:
  - `n_panels::Int64`: Number of panels in aerodynamic mesh
  - `spanwise_distribution`::PanelDistribution: see: [PanelDistribution](@ref)
  - `spanwise_direction::MVec3`: Wing span direction vector
  - `sections::Vector{Section}`: List of wing sections, see: [Section](@ref)
  -  refined_sections::Vector{Section}
  - `remove_nan::Bool`: Wether to remove the NaNs from interpolations or not
- Additional fields:
  - `circle_center_z::Vector{Float64}`: Center of circle coordinates
  - gamma_tip::Float64: Angle between the body frame z axis and the vector going from the kite circular shape center to the wing tip.
  - `inertia_tensor`::Matrix{Float64}: see: [`calculate_inertia_tensor`](@ref)
  - radius::Float64: Radius of curvature
  - le_interp::NTuple{3, Extrapolation}: see: [Extrapolation](https://juliamath.github.io/Interpolations.jl/stable/extrapolation/)
  - te_interp::NTuple{3, Extrapolation}
  - area_interp::Extrapolation

"""
mutable struct RamAirWing <: AbstractWing
    n_panels::Int64
    spanwise_distribution::PanelDistribution
    panel_props::PanelProperties
    spanwise_direction::MVec3
    sections::Vector{Section}
    refined_sections::Vector{Section}
    remove_nan::Bool
    
    # Additional fields for RamAirWing
    non_deformed_sections::Vector{Section}
    mass::Float64
    gamma_tip::Float64
    inertia_tensor::Matrix{Float64}
    center_of_mass::Vector{Float64}
    radius::Float64
    le_interp::NTuple{3, Extrapolation}
    te_interp::NTuple{3, Extrapolation}
    area_interp::Extrapolation
    theta_dist::Vector{Float64}
    delta_dist::Vector{Float64}
end

"""
    RamAirWing(obj_path, dat_path; alpha=0.0, crease_frac=0.75, wind_vel=10., mass=1.0, 
             n_panels=54, n_sections=n_panels+1, spanwise_distribution=UNCHANGED, 
             spanwise_direction=[0.0, 1.0, 0.0], remove_nan::Bool=true)

Constructor for a [RamAirWing](@ref) that allows to use an `.obj` and a `.dat` file as input.

# Parameters
- obj_path: Path to the `.obj` file used for creating the geometry
- dat_path: Path to the `.dat` file, a standard format for 2d foil geometry

# Keyword Parameters
- alpha=0.0: Angle of attack of each segment relative to the x axis [rad]
- crease_frac=0.75: The x coordinate around which the trailing edge rotates on a normalized 2d foil, 
                    used in the xfoil polar generation
- wind_vel=10.0: Apparent wind speed in m/s, used in the xfoil polar generation
- mass=1.0: Mass of the wing in kg, used for the inertia calculations 
- `n_panels`=56: Number of panels.
- `n_sections`=n_panels+1: Number of sections (there is a section on each side of each panel.)
- `spanwise_distribution`=UNCHANGED: see: [PanelDistribution](@ref)
- `spanwise_direction`=[0.0, 1.0, 0.0]
- `remove_nan::Bool`: Wether to remove the NaNs from interpolations or not
"""
function RamAirWing(obj_path, dat_path; alpha=0.0, crease_frac=0.75, wind_vel=10., mass=1.0, 
                  n_panels=56, n_sections=n_panels+1, spanwise_distribution=UNCHANGED, 
                  spanwise_direction=[0.0, 1.0, 0.0], remove_nan=true, align_to_principal=false,
                  alpha_range=deg2rad.(-5:1:20), delta_range=deg2rad.(-5:1:20), interp_steps=n_sections
                  )

    !isapprox(spanwise_direction, [0.0, 1.0, 0.0]) && throw(ArgumentError("Spanwise direction has to be [0.0, 1.0, 0.0], not $spanwise_direction"))

    # Load or create polars
    (!endswith(dat_path, ".dat")) && (dat_path *= ".dat")
    (!isfile(dat_path)) && error("DAT file not found: $dat_path")
    polar_path = dat_path[1:end-4] * "_polar.bin"

    (!endswith(obj_path, ".obj")) && (obj_path *= ".obj")
    (!isfile(obj_path)) && error("OBJ file not found: $obj_path")
    info_path = obj_path[1:end-4] * "_info.bin"

    if !ispath(info_path)
        @info "Reading $obj_path"
        vertices, faces = read_faces(obj_path)
        center_of_mass = center_to_com!(vertices, faces)
        inertia_tensor = calculate_inertia_tensor(vertices, faces, mass, zeros(3))

        if align_to_principal
            inertia_tensor, R_b_p = calc_inertia_y_rotation(inertia_tensor)
            circle_center_z, radius, gamma_tip = find_circle_center_and_radius(vertices)
            le_interp, te_interp, area_interp = create_interpolations(vertices, circle_center_z, radius, gamma_tip, R_b_p; interp_steps)
        else
            circle_center_z, radius, gamma_tip = find_circle_center_and_radius(vertices)
            le_interp, te_interp, area_interp = create_interpolations(vertices, circle_center_z, radius, gamma_tip, I(3); interp_steps)
        end

        @info "Writing $info_path"
        serialize(info_path, (inertia_tensor, center_of_mass, radius, gamma_tip, 
            le_interp, te_interp, area_interp))
    end

    @info "Loading kite info from $info_path and polars from $polar_path"
    try
        (inertia_tensor::Matrix, center_of_mass::Vector,
            radius::Real, gamma_tip::Real, le_interp, te_interp, area_interp) = deserialize(info_path)

        if !ispath(polar_path)
            width = 2gamma_tip * radius
            area = area_interp(gamma_tip)
            @eval Main begin
                foil_path, polar_path, v_wind, area, width, x_turn, alpha_range, delta_range =
                    $dat_path, $polar_path, $wind_vel, $gamma_tip, $width, $crease_frac, $alpha_range, $delta_range
                include(joinpath(dirname(@__FILE__), "../scripts/polars.jl"))
            end
        end

        (alpha_range, delta_range, cl_matrix::Matrix, cd_matrix::Matrix, cm_matrix::Matrix) = deserialize(polar_path)
        if remove_nan
            interpolate_matrix_nans!(cl_matrix)
            interpolate_matrix_nans!(cd_matrix)
            interpolate_matrix_nans!(cm_matrix)
        end

        # Create sections
        sections = Section[]
        refined_sections = Section[]
        non_deformed_sections = Section[]
        for gamma in range(-gamma_tip, gamma_tip, n_sections)
            aero_data = (collect(alpha_range), collect(delta_range), cl_matrix, cd_matrix, cm_matrix)
            LE_point = [le_interp[i](gamma) for i in 1:3]
            TE_point = [te_interp[i](gamma) for i in 1:3]
            push!(sections, Section(LE_point, TE_point, POLAR_MATRICES, aero_data))
            push!(refined_sections, Section(LE_point, TE_point, POLAR_MATRICES, aero_data))
            push!(non_deformed_sections, Section(LE_point, TE_point, POLAR_MATRICES, aero_data))
        end

        panel_props = PanelProperties{n_panels}()

        RamAirWing(n_panels, spanwise_distribution, panel_props, spanwise_direction, sections, 
            refined_sections, remove_nan, non_deformed_sections,
            mass, gamma_tip, inertia_tensor, center_of_mass, radius,
            le_interp, te_interp, area_interp, zeros(n_panels), zeros(n_panels))

    catch e
        if e isa BoundsError
            @error "Delete $info_path and $polar_path and try again."
        end
        rethrow(e)
    end
end

function smooth_deform!(wing::RamAirWing, theta_angles::AbstractVector, delta_angles::AbstractVector)
    !(wing.n_panels % length(theta_angles) == 0) && throw(ArgumentError("Number of angles has to be a multiple of number of panels"))
    panels = wing.panels
    theta_dist = wing.theta_dist
    delta_dist = wing.delta_dist

    dist_idx = 1
    for angle_idx in eachindex(theta_angles)
        for _ in 1:(wing.n_panels ÷ length(theta_angles))
            theta_dist[dist_idx] = theta_angles[angle_idx]
            delta_dist[dist_idx] = delta_angles[angle_idx]
            dist_idx += 1
        end
    end
    @assert (dist_idx == wing.n_panels)

    window_size = wing.n_panels ÷ length(theta_angles)
    if length(panels) > window_size
        smoothed = copy(theta_dist)
        for i in (window_size÷2 + 1):(length(panels) - window_size÷2)
            smoothed[i] = mean(theta_dist[(i - window_size÷2):(i + window_size÷2)])
        end
        theta_dist .= smoothed
        for i in (window_size÷2 + 1):(length(panels) - window_size÷2)
            smoothed[i] = mean(delta_dist[(i - window_size÷2):(i + window_size÷2)])
        end
        delta_dist .= smoothed
    end
    deform!(wing)
    return nothing
end

"""
    deform!(wing::RamAirWing, theta_dist::AbstractVector, delta_dist::AbstractVector; width)

Deform wing by applying theta and delta distributions.

# Arguments
- `wing`: RamAirWing to deform
- `theta_dist`: the angle distribution between of the kite and the body x-axis in radians of each panel
- `delta_dist`: the deformation of the trailing edges of each panel

# Effects
Updates wing.sections with deformed geometry
"""
function deform!(wing::RamAirWing, theta_dist::AbstractVector, delta_dist::AbstractVector)
    !(length(theta_dist) == wing.n_panels) && throw(ArgumentError("theta_dist and panels are of different lengths"))
    !(length(delta_dist) == wing.n_panels) && throw(ArgumentError("delta_dist and panels are of different lengths"))
    wing.theta_dist .= theta_dist
    wing.delta_dist .= delta_dist

    deform!(wing)
end

function deform!(wing::RamAirWing)
    local_y = zeros(MVec3)
    chord = zeros(MVec3)
    normal = zeros(MVec3)

    for i in 1:wing.n_panels
        section1 = wing.non_deformed_sections[i]
        section2 = wing.non_deformed_sections[i+1]
        local_y .= normalize(section1.LE_point - section2.LE_point)
        chord .= section1.TE_point .- section1.LE_point
        normal .= chord × local_y
        @. wing.sections[i].TE_point = section1.LE_point + cos(wing.theta_dist[i]) * chord - sin(wing.theta_dist[i]) * normal
    end
    return nothing
end
