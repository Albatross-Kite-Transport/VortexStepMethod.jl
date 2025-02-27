
"""
    set_plot_style(titel_size=16)

Set the default style for plots using LaTeX.

# Arguments:
- `titel_size`: size of the plot title in points (default: 16)
"""
function set_plot_style(titel_size=16)
    rcParams = plt.PyDict(plt.matplotlib."rcParams")
    rcParams["text.usetex"] = true
    rcParams["font.family"] = "serif"
    rcParams["font.serif"] = ["Computer Modern Roman"]
    rcParams["axes.titlesize"] = titel_size
    rcParams["axes.labelsize"] = 12
    rcParams["axes.linewidth"] = 1
    rcParams["lines.linewidth"] = 1
    rcParams["lines.markersize"] = 6
    rcParams["xtick.labelsize"] = 10
    rcParams["ytick.labelsize"] = 10
    rcParams["legend.fontsize"] = 10
    rcParams["figure.titlesize"] = 16
    rcParams["pgf.texsystem"] = "pdflatex"  # Use pdflatex
    rcParams["pgf.rcfonts"] = false
    rcParams["figure.figsize"] = (10, 6)  # Default figure size
end


"""
    save_plot(fig, save_path, title; data_type=".pdf")

Save a plot to a file.

# Arguments
- `fig`: Plots figure object
- `save_path`: Path to save the plot
- `title`: Title of the plot

# Keyword arguments
- `data_type`: File extension (default: ".pdf")
"""
function VortexStepMethod.save_plot(fig, save_path, title; data_type=".pdf")
    isnothing(save_path) && throw(ArgumentError("save_path should be provided"))

    !isdir(save_path) && mkpath(save_path)
    full_path = joinpath(save_path, title * data_type)

    @debug "Attempting to save figure to: $full_path"
    @debug "Current working directory: $(pwd())"

    try
        fig.savefig(full_path)
        @debug "Figure saved as $data_type"

        if isfile(full_path)
            @debug "File successfully saved to $full_path"
            @debug "File size: $(filesize(full_path)) bytes"
        else
            @info "File does not exist after save attempt: $full_path"
        end
    catch e
        @error "Error saving figure: $e"
        @error "Error type: $(typeof(e))"
        rethrow(e)
    end
end

"""
    show_plot(fig; dpi=130)

Display a plot at specified DPI.

# Arguments
- `fig`: Plots figure object

# Keyword arguments
- `dpi`: Dots per inch for the figure (default: 130)
"""
function VortexStepMethod.show_plot(fig; dpi=130)
    plt.display(fig)
end

"""
    plot_line_segment!(ax, segment, color, label; width=3)

Plot a line segment in 3D with arrow.

# Arguments
- `ax`: Plot axis
- `segment`: Array of two points defining the segment
- `color`: Color of the segment
- `label`: Label for the legend

# Keyword Arguments
- `width`: Line width (default: 3)
"""
function plot_line_segment!(ax, segment, color, label; width=3)
    ax.plot(
        [segment[1][1], segment[2][1]],
        [segment[1][2], segment[2][2]],
        [segment[1][3], segment[2][3]],
        color=color, label=label, linewidth=width
    )

    dir = segment[2] - segment[1]
    ax.quiver(
        [segment[1][1]], [segment[1][2]], [segment[1][3]],
        [dir[1]], [dir[2]], [dir[3]],
        color=color
    )
end

"""
    set_axes_equal!(ax; zoom=1.8)

Set 3D plot axes to equal scale.

# Arguments
- ax: 3D plot axis

# Keyword arguments
zoom: zoom factor (default: 1.8)
"""
function set_axes_equal!(ax; zoom=1.8)
    x_lims = ax.get_xlim3d() ./ zoom
    y_lims = ax.get_ylim3d() ./ zoom
    z_lims = ax.get_zlim3d() ./ zoom

    x_range = abs(x_lims[2] - x_lims[1])
    y_range = abs(y_lims[2] - y_lims[1])
    z_range = abs(z_lims[2] - z_lims[1])

    max_range = max(x_range, y_range, z_range)

    x_mid = mean(x_lims)
    y_mid = mean(y_lims)
    z_mid = mean(z_lims)

    ax.set_xlim3d((x_mid - max_range / 2, x_mid + max_range / 2))
    ax.set_ylim3d((y_mid - max_range / 2, y_mid + max_range / 2))
    ax.set_zlim3d((z_mid - max_range / 2, z_mid + max_range / 2))
end

"""
    create_geometry_plot(body_aero::BodyAerodynamics, title, view_elevation, view_azimuth; zoom=1.8)

Create a 3D plot of wing geometry including panels and filaments.

# Arguments
- body_aero: struct of type BodyAerodynamics
- title: plot title
- view_elevation: initial view elevation angle [°]
- view_azimuth: initial view azimuth angle [°]

# Keyword arguments
- zoom: zoom factor (default: 1.8)
"""
function create_geometry_plot(body_aero::BodyAerodynamics, title, view_elevation, view_azimuth; 
                              zoom=1.8)
    set_plot_style(28)

    panels = body_aero.panels
    va = isa(body_aero.va, Tuple) ? body_aero.va[1] : body_aero.va

    # Extract geometric data
    corner_points = [panel.corner_points for panel in panels]
    control_points = [panel.control_point for panel in panels]
    aero_centers = [panel.aero_center for panel in panels]

    # Create plot
    fig = plt.figure(figsize=(14, 14))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_title(title)

    # Plot panels
    legend_used = Dict{String,Bool}()
    for (i, panel) in enumerate(panels)
        # Plot panel edges and surfaces
        corners = panel.corner_points
        x_corners = corners[1, :]
        y_corners = corners[2, :]
        z_corners = corners[3, :]

        push!(x_corners, x_corners[1])
        push!(y_corners, y_corners[1])
        push!(z_corners, z_corners[1])

        ax.plot(x_corners,
            y_corners,
            z_corners,
            color=:grey,
            linewidth=1,
            label=i == 1 ? "Panel Edges" : "")

        # Plot control points and aerodynamic centers
        ax.scatter([control_points[i][1]], [control_points[i][2]], [control_points[i][3]],
            color=:green, label=i == 1 ? "Control Points" : "")
        ax.scatter([aero_centers[i][1]], [aero_centers[i][2]], [aero_centers[i][3]],
            color=:blue, label=i == 1 ? "Aerodynamic Centers" : "")

        # Plot filaments
        filaments = calculate_filaments_for_plotting(panel)
        legends = ["Bound Vortex", "side1", "side2", "wake_1", "wake_2"]

        for (filament, legend) in zip(filaments, legends)
            x1, x2, color = filament
            @debug "Legend: $legend"
            show_legend = !get(legend_used, legend, false)
            plot_line_segment!(ax, [x1, x2], color, show_legend ? legend : "")
            legend_used[legend] = true
        end
    end

    # Plot velocity vector
    max_chord = maximum(panel.chord for panel in panels)
    va_mag = norm(va)
    va_vector_begin = -2 * max_chord * va / va_mag
    va_vector_end = va_vector_begin + 1.5 * va / va_mag
    plot_line_segment!(ax, [va_vector_begin, va_vector_end], :lightblue, "va")

    # Add legends for the first occurrence of each label
    handles, labels = ax.get_legend_handles_labels()
    by_label = Dict(zip(labels, handles))
    ax.legend(values(by_label), keys(by_label), bbox_to_anchor=(0, 0, 1.1, 1))

    # Set labels and make axes equal
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    set_axes_equal!(ax; zoom)

    # Set the initial view
    ax.view_init(elev=view_elevation, azim=view_azimuth)

    # Ensure the figure is fully rendered
    # fig.canvas.draw()
    plt.tight_layout(rect=(0, 0, 1, 0.97))

    return fig
end

"""
    plot_geometry(body_aero::BodyAerodynamics, title; 
                  data_type=".pdf", save_path=nothing, 
                  is_save=false, is_show=false, 
                  view_elevation=15, view_azimuth=-120)

Plot wing geometry from different viewpoints and optionally save/show plots.

# Arguments:
- `body_aero`: the [BodyAerodynamics](@ref) to plot
- title: plot title

# Keyword arguments:
- `data_type``: string with the file type postfix (default: ".pdf")
- `save_path`: path for saving the graphic (default: `nothing`)
- `is_save`: boolean value, indicates if the graphic shall be saved (default: `false`)
- `is_show`: boolean value, indicates if the graphic shall be displayed (default: `false`)
- `view_elevation`: initial view elevation angle (default: 15) [°]
- `view_azimuth`: initial view azimuth angle (default: -120) [°]

"""
function VortexStepMethod.plot_geometry(body_aero::BodyAerodynamics, title;
    data_type=".pdf",
    save_path=nothing,
    is_save=false,
    is_show=false,
    view_elevation=15,
    view_azimuth=-120)

    if is_save
        plt.ioff()
        # Angled view
        fig = create_geometry_plot(body_aero, "$(title)_angled_view", 15, -120)
        save_plot(fig, save_path, "$(title)_angled_view", data_type=data_type)

        # Top view
        fig = create_geometry_plot(body_aero, "$(title)_top_view", 90, 0)
        save_plot(fig, save_path, "$(title)_top_view", data_type=data_type)

        # Front view
        fig = create_geometry_plot(body_aero, "$(title)_front_view", 0, 0)
        save_plot(fig, save_path, "$(title)_front_view", data_type=data_type)

        # Side view
        fig = create_geometry_plot(body_aero, "$(title)_side_view", 0, -90)
        save_plot(fig, save_path, "$(title)_side_view", data_type=data_type)
    end

    if is_show
        plt.ion()
        fig = create_geometry_plot(body_aero, title, view_elevation, view_azimuth)
        plt.display(fig)
    else
        fig = create_geometry_plot(body_aero, title, view_elevation, view_azimuth)
    end
    fig
end

"""
    plot_distribution(y_coordinates_list, results_list, label_list;  
                      title="spanwise_distribution", data_type=".pdf",
                      save_path=nothing, is_save=false, is_show=true)

Plot spanwise distributions of aerodynamic properties.

# Arguments
- `y_coordinates_list`: List of spanwise coordinates
- `results_list`: List of result dictionaries
- `label_list`: List of labels for different results

# Keyword arguments
- `title`: Plot title (default: "spanwise_distribution")
- `data_type`: File extension for saving (default: ".pdf")
- `save_path`: Path to save plots (default: nothing)
- `is_save`: Whether to save plots (default: false)
- `is_show`: Whether to display plots (default: true)
"""
function VortexStepMethod.plot_distribution(y_coordinates_list, results_list, label_list;
    title="spanwise_distribution",
    data_type=".pdf",
    save_path=nothing,
    is_save=false,
    is_show=true)

    length(results_list) == length(label_list) || throw(ArgumentError(
        "Number of results ($(length(results_list))) must match number of labels ($(length(label_list)))"
    ))

    # Set the plot style
    set_plot_style()

    # Initializing plot
    fig, axs = plt.subplots(3, 3, figsize=(16, 10))
    fig.suptitle(title, fontsize=16)

    # CL plot
    for (y_coordinates_i, result_i, label_i) in zip(y_coordinates_list, results_list, label_list)
        value = "$(round(result_i["cl"], digits=2))"
        if label_i == "LLT"
            label = label_i * L" $~C_L$: " * value
        else
            label = label_i * L" $C_L$: " * value
        end
        axs[1, 1].plot(
            y_coordinates_i,
            result_i["cl_distribution"],
            label=label
        )
    end
    axs[1, 1].set_title(L"$C_L$ Distribution", size=16)
    axs[1, 1].set_xlabel(L"Spanwise Position $y/b$")
    axs[1, 1].set_ylabel(L"Lift Coefficient $C_L$")
    axs[1, 1].legend()

    # CD plot
    for (y_coordinates_i, result_i, label_i) in zip(y_coordinates_list, results_list, label_list)
        value = "$(round(result_i["cl"], digits=2))"
        if label_i == "LLT"
            label = label_i * L" $~C_D$: " * value
        else
            label = label_i * L" $C_D$: " * value
        end
        axs[1, 2].plot(
            y_coordinates_i,
            result_i["cd_distribution"],
            label=label
        )
    end
    axs[1, 2].set_title(L"$C_D$ Distribution", size=16)
    axs[1, 2].set_xlabel(L"Spanwise Position $y/b$")
    axs[1, 2].set_ylabel(L"Drag Coefficient $C_D$")
    axs[1, 2].legend()

    # Gamma Distribution
    for (y_coordinates_i, result_i, label_i) in zip(y_coordinates_list, results_list, label_list)
        axs[1, 3].plot(
            y_coordinates_i,
            result_i["gamma_distribution"],
            label=label_i
        )
    end
    axs[1, 3].set_title(L"\Gamma~Distribution", size=16)
    axs[1, 3].set_xlabel(L"Spanwise Position $y/b$")
    axs[1, 3].set_ylabel(L"Circulation~\Gamma")
    axs[1, 3].legend()

    # Geometric Alpha
    for (y_coordinates_i, result_i, label_i) in zip(y_coordinates_list, results_list, label_list)
        axs[2, 1].plot(
            y_coordinates_i,
            result_i["alpha_geometric"],
            label=label_i
        )
    end
    axs[2, 1].set_title(L"$\alpha$ Geometric", size=16)
    axs[2, 1].set_xlabel(L"Spanwise Position $y/b$")
    axs[2, 1].set_ylabel(L"Angle of Attack $\alpha$ (deg)")
    axs[2, 1].legend()

    # Calculated/ Corrected Alpha
    for (y_coordinates_i, result_i, label_i) in zip(y_coordinates_list, results_list, label_list)
        axs[2, 2].plot(
            y_coordinates_i,
            result_i["alpha_at_ac"],
            label=label_i
        )
    end
    axs[2, 2].set_title(L"$\alpha$ result (corrected to aerodynamic center)", size=16)
    axs[2, 2].set_xlabel(L"Spanwise Position $y/b$")
    axs[2, 2].set_ylabel(L"Angle of Attack $\alpha$ (deg)")
    axs[2, 2].legend()

    # Uncorrected Alpha plot
    for (y_coordinates_i, result_i, label_i) in zip(y_coordinates_list, results_list, label_list)
        axs[2, 3].plot(
            y_coordinates_i,
            result_i["alpha_uncorrected"],
            label=label_i
        )
    end
    axs[2, 3].set_title(L"$\alpha$ Uncorrected (if VSM, at the control point)", size=16)
    axs[2, 3].set_xlabel(L"Spanwise Position $y/b$")
    axs[2, 3].set_ylabel(L"Angle of Attack $\alpha$ (deg)")
    axs[2, 3].legend()

    # Force Components
    for (idx, component) in enumerate(["x", "y", "z"])
        axs[3, idx].set_title("Force in $component direction", size=16)
        axs[3, idx].set_xlabel(L"Spanwise Position $y/b$")
        axs[3, idx].set_ylabel(raw"$F_\mathrm" * "{$component}" * raw"$")
        for (y_coords, results, label) in zip(y_coordinates_list, results_list, label_list)
            # Extract force components for the current direction (idx)
            forces = results["F_distribution"][idx, :]
            # Verify dimensions match
            if length(y_coords) != length(forces)
                @warn "Dimension mismatch in force plotting" length(y_coords) length(forces) component
                continue  # Skip this component instead of throwing error
            end
            space = ""
            if label == "LLT"
                space = "~"
            end
            axs[3, idx].plot(
                y_coords,
                forces,
                label="$label" * space * raw"$~\Sigma~F_\mathrm" * "{$component}:~" *
                      raw"$" * "$(round(results["F$component"], digits=2)) N"
            )
            axs[3, idx].legend()
        end
    end

    fig.tight_layout()

    # Save and show plot
    if is_save
        save_plot(fig, save_path, title, data_type=data_type)
    end

    if is_show
        show_plot(fig)
    end

    return fig
end

"""
    generate_polar_data(solver, body_aero::BodyAerodynamics, angle_range;     
                        angle_type="angle_of_attack", angle_of_attack=0.0,
                        side_slip=0.0, v_a=10.0)

Generate polar data for aerodynamic analysis over a range of angles.

# Arguments
- `solver`: Aerodynamic solver object
- `body_aero`: Wing aerodynamics struct
- `angle_range`: Range of angles to analyze

# Keyword arguments
- `angle_type`: Type of angle variation ("angle_of_attack" or "side_slip")
- `angle_of_attack`: Initial angle of attack [rad]
- `side_slip`: Initial side slip angle in [rad]
- `v_a`: norm of apparent wind speed [m/s]

# Returns
- Tuple of polar data array and Reynolds number
"""
function generate_polar_data(
    solver,
    body_aero::BodyAerodynamics,
    angle_range;
    angle_type="angle_of_attack",
    angle_of_attack=0.0,
    side_slip=0.0,
    v_a=10.0
)
    n_panels = length(body_aero.panels)
    n_angles = length(angle_range)

    # Initialize arrays
    cl = zeros(n_angles)
    cd = zeros(n_angles)
    cs = zeros(n_angles)
    gamma_distribution = zeros(n_angles, n_panels)
    cl_distribution = zeros(n_angles, n_panels)
    cd_distribution = zeros(n_angles, n_panels)
    cs_distribution = zeros(n_angles, n_panels)
    reynolds_number = zeros(n_angles)

    # Previous gamma for initialization
    gamma = nothing

    set_plot_style()

    for (i, angle_i) in enumerate(angle_range)
        # Set angle based on type
        if angle_type == "angle_of_attack"
            α = deg2rad(angle_i)
            β = side_slip
        elseif angle_type == raw"side_slip"
            α = angle_of_attack
            β = deg2rad(angle_i)
        else
            throw(ArgumentError("angle_type must be 'angle_of_attack' or 'side_slip'"))
        end

        # Update inflow conditions
        set_va!(
            body_aero,
            [
                cos(α) * cos(β),
                sin(β),
                sin(α)
            ] * v_a
        )

        # Solve and store results
        results = solve(solver, body_aero, gamma_distribution[i, :])

        cl[i] = results["cl"]
        cd[i] = results["cd"]
        cs[i] = results["cs"]
        gamma_distribution[i, :] = results["gamma_distribution"]
        cl_distribution[i, :] = results["cl_distribution"]
        cd_distribution[i, :] = results["cd_distribution"]
        cs_distribution[i, :] = results["cs_distribution"]
        reynolds_number[i] = results["Rey"]

        # Store gamma for next iteration
        gamma = gamma_distribution[i, :]
    end

    polar_data = [
        angle_range,
        cl,
        cd,
        cs,
        gamma_distribution,
        cl_distribution,
        cd_distribution,
        cs_distribution,
        reynolds_number
    ]

    return polar_data, reynolds_number[1]
end

"""
    plot_polars(solver_list, body_aero_list, label_list;
    literature_path_list=String[], angle_range=range(0, 20, 2), angle_type="angle_of_attack", 
    angle_of_attack=0.0, side_slip=0.0, v_a=10.0, 
    title="polar", data_type=".pdf", save_path=nothing, 
    is_save=true, is_show=true)

Plot polar data comparing different solvers and configurations.

# Arguments
- `solver_list`: List of aerodynamic solvers
- `body_aero_list`: List of wing aerodynamics objects
- `label_list`: List of labels for each configuration

# Keyword arguments
- `literature_path_list`: Optional paths to literature data files
- `angle_range`: Range of angles to analyze [°]
- `angle_type`: "`angle_of_attack`" or "`side_slip`"; (default: `angle_of_attack`) 
- `angle_of_attack:` AoA to be used for plotting the polars (default: 0.0) [rad]
- `side_slip`: side slip angle (default: 0.0) [rad]
- v_a: norm of apparent wind speed (default: 10.0) [m/s]
- title: plot title
- `data_type`: File extension for saving (default: ".pdf")
- `save_path`: Path to save plots (default: nothing)
- `is_save`: Whether to save plots (default: true)
- `is_show`: Whether to display plots (default: true)
"""
function VortexStepMethod.plot_polars(
    solver_list,
    body_aero_list,
    label_list;
    literature_path_list=String[],
    angle_range=range(0, 20, 2),
    angle_type="angle_of_attack",
    angle_of_attack=0.0,
    side_slip=0.0,
    v_a=10.0,
    title="polar",
    data_type=".pdf",
    save_path=nothing,
    is_save=true,
    is_show=true
)
    # Validate inputs
    total_cases = length(body_aero_list) + length(literature_path_list)
    if total_cases != length(label_list) || length(solver_list) != length(body_aero_list)
        throw(ArgumentError("Mismatch in number of solvers ($(length(solver_list))), " *
                            "cases ($total_cases), and labels ($(length(label_list)))"))
    end
    main_title = replace(title, " " => "_")

    # Generate polar data
    polar_data_list = []
    for (i, (solver, body_aero)) in enumerate(zip(solver_list, body_aero_list))
        polar_data, rey = generate_polar_data(
            solver, body_aero, angle_range;
            angle_type,
            angle_of_attack,
            side_slip,
            v_a
        )
        push!(polar_data_list, polar_data)
        # Update label with Reynolds number
        label_list[i] = "$(label_list[i]) Re = $(round(Int64, rey*1e-5))e5"
    end

    # Load literature data if provided
    if !isempty(literature_path_list)
        for path in literature_path_list
            # Read all data first
            data = readdlm(path, ',')
            # Skip the header row by taking data from row 2 onwards
            data = data[2:end, :]
            push!(polar_data_list, [data[:, 3], data[:, 1], data[:, 2], zeros(length(data[:, 1]))])
        end
    end

    # Initializing plot
    fig, axs = plt.subplots(2, 2, figsize=(14, 14))

    # Number of computational results (excluding literature)
    n_solvers = length(solver_list)
    for (i, (polar_data, label)) in enumerate(zip(polar_data_list, label_list))
        if i < n_solvers
            linestyle = "-"
            marker = "*"
            markersize = 7
        else
            linestyle = "-"
            marker = "."
            markersize = 5
        end
        if contains(label, "LLT")
            label = replace(label, "e5" => raw"\cdot10^5")
            label = replace(label, " " => raw"~")
            label = replace(label, "LLT" => raw"\mathrm{LLT}{~\,}")
            label = raw"$" * label * raw"$"
        else
            label = replace(label, "e5" => raw"\cdot10^5")
            label = replace(label, " " => "~")
            label = replace(label, "VSM" => raw"\mathrm{VSM}")
            label = raw"$" * label * raw"$"
        end
        axs[1, 1].plot(
            polar_data[1],
            polar_data[2],
            label=label,
            linestyle=linestyle,
            marker=marker,
            markersize=markersize,
        )
        # Limit y-range if CL > 10
        if maximum(polar_data[2]) > 10
            axs[1, 1].set_ylim([-0.5, 2])
        end
        title = raw"$C_L" * raw"$" * " vs $angle_type [°]"
        axs[1, 1].set_title(title)
        axs[1, 1].set_xlabel("$angle_type [°]")
        axs[1, 1].set_ylabel(L"$C_L$")
        axs[1, 1].legend()
    end

    for (i, (polar_data, label)) in enumerate(zip(polar_data_list, label_list))
        if i < n_solvers
            linestyle = "-"
            marker = "*"
            markersize = 7
        else
            linestyle = "-"
            marker = "."
            markersize = 5
        end
        if contains(label, "LLT")
            label = replace(label, "e5" => raw"\cdot10^5")
            label = replace(label, " " => raw"~")
            label = replace(label, "LLT" => raw"\mathrm{LLT}{~\,}")
            label = raw"$" * label * raw"$"
        else
            label = replace(label, "e5" => raw"\cdot10^5")
            label = replace(label, " " => "~")
            label = replace(label, "VSM" => raw"\mathrm{VSM}")
            label = raw"$" * label * raw"$"
        end
        axs[1, 2].plot(
            polar_data[1],
            polar_data[3],
            label=label,
            linestyle=linestyle,
            marker=marker,
            markersize=markersize,
        )
        # Limit y-range if CL > 10
        if maximum(polar_data[2]) > 10
            axs[1, 2].set_ylim([-0.5, 2])
        end
        title = raw"$C_D" * raw"$" * " vs $angle_type [°]"
        axs[1, 2].set_title(title)
        axs[1, 2].set_xlabel("$angle_type [°]")
        axs[1, 2].set_ylabel(L"$C_D$")
        axs[1, 2].legend()
    end


    for (i, (polar_data, label)) in enumerate(zip(polar_data_list, label_list))
        if i < n_solvers
            linestyle = "-"
            marker = "*"
            markersize = 7
        else
            linestyle = "-"
            marker = "."
            markersize = 5
        end
        if contains(label, "LLT")
            label = replace(label, "e5" => raw"\cdot10^5")
            label = replace(label, " " => raw"~")
            label = replace(label, "LLT" => raw"\mathrm{LLT}{~\,}")
            label = raw"$" * label * raw"$"
        else
            label = replace(label, "e5" => raw"\cdot10^5")
            label = replace(label, " " => "~")
            label = replace(label, "VSM" => raw"\mathrm{VSM}")
            label = raw"$" * label * raw"$"
        end
        axs[2, 1].plot(
            polar_data[1],
            polar_data[4],
            label=label,
            linestyle=linestyle,
            marker=marker,
            markersize=markersize,
        )
        # Limit y-range if CL > 10
        if maximum(polar_data[2]) > 10
            axs[2, 1].set_ylim([-0.5, 2])
        end
        title = raw"$C_S" * raw"$" * " vs $angle_type [°]"
        axs[2, 1].set_title(title)
        axs[2, 1].set_xlabel("$angle_type [°]")
        axs[2, 1].set_ylabel(L"$C_S$")
        axs[2, 1].legend()
    end

    for (i, (polar_data, label)) in enumerate(zip(polar_data_list, label_list))
        if i < n_solvers
            linestyle = "-"
            marker = "*"
            markersize = 7
        else
            linestyle = "-"
            marker = "."
            markersize = 5
        end
        if contains(label, "LLT")
            label = replace(label, "e5" => raw"\cdot10^5")
            label = replace(label, " " => raw"~")
            label = replace(label, "LLT" => raw"\mathrm{LLT}{~\,}")
            label = raw"$" * label * raw"$"
        else
            label = replace(label, "e5" => raw"\cdot10^5")
            label = replace(label, " " => "~")
            label = replace(label, "VSM" => raw"\mathrm{VSM}")
            label = raw"$" * label * raw"$"
        end
        axs[2, 2].plot(
            polar_data[3],
            polar_data[2],
            label=label,
            linestyle=linestyle,
            marker=marker,
            markersize=markersize,
        )
        # Limit y-range if CL > 10
        if maximum(polar_data[2]) > 10 || maximum(polar_data[3]) > 10
            axs[2, 2].set_ylim([-0.5, 2])
            axs[2, 2].set_xlim([-0.5, 2])
        end
        title = raw"$C_L" * raw"$" * " vs " * raw"$C_D" * raw"$"
        axs[2, 2].set_title(title)
        axs[2, 2].set_xlabel(L"$C_D$")
        axs[2, 2].set_ylabel(L"$C_L$")
        axs[2, 2].legend()
    end

    fig.tight_layout(h_pad=3.5, rect=(0.01, 0.01, 0.99, 0.99))

    # Save and show plot
    if is_save && !isnothing(save_path)
        save_plot(fig, save_path, main_title; data_type)
    end

    if is_show
        show_plot(fig)
    end

    return fig
end