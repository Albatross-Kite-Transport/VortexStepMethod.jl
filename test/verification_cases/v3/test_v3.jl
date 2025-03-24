function get_CAD_matching_uri()
    CAD_matching_pos_values =
            [0.000000e00 0.000000e00 0.000000e00
            1.538773e03 4.113307e03 5.530496e03
            -1.762300e01 3.967978e03 6.471622e03
            -2.376830e02 3.134335e03 7.476759e03
            -3.837330e02 1.959729e03 8.078914e03
            -4.565620e02 6.642520e02 8.339101e03
            -4.565620e02 -6.642520e02 8.339101e03
            -3.837330e02 -1.959729e03 8.078914e03
            -2.376830e02 -3.134335e03 7.476759e03
            -1.762300e01 -3.967978e03 6.471622e03
            1.538773e03 -4.113307e03 5.530496e03
            1.703467e03 -3.955506e03 6.467819e03
            2.002516e03 -3.116753e03 7.454254e03
            2.110145e03 -1.946587e03 8.041739e03
            2.158559e03 -6.607210e02 8.294064e03
            2.158559e03 6.607210e02 8.294064e03
            2.110145e03 1.946587e03 8.041739e03
            2.002516e03 3.116753e03 7.454254e03
            1.703467e03 3.955506e03 6.467819e03
            8.595800e02 4.139660e03 5.654227e03
            8.595800e02 -4.139660e03 5.654227e03
            3.138480e02 0.000000e00 1.252129e03
            3.875460e02 6.909170e02 1.479174e03
            1.053321e03 1.772499e03 4.343864e03
            1.441205e03 2.708913e03 5.601785e03
            1.528031e03 1.338349e03 5.966178e03
            3.875460e02 -6.909170e02 1.479174e03
            1.053321e03 -1.772499e03 4.343864e03
            1.441205e03 -2.708913e03 5.601785e03
            1.528031e03 -1.338349e03 5.966178e03
            -6.552600e01 1.321471e03 4.213046e03
            -4.248900e01 2.046976e03 5.525725e03
            -9.173800e01 1.262274e03 5.961848e03
            -6.552600e01 -1.321471e03 4.213046e03
            -4.248900e01 -2.046976e03 5.525725e03
            -9.173800e01 -1.262274e03 5.961848e03]
    return CAD_matching_pos_values
end
function struct2aero_geometry(coord_struc)
    coord = zeros(20, 3)
    coord[1, :] = coord_struc[21, :]
    coord[2, :] = coord_struc[11, :]
    coord[3, :] = coord_struc[10, :]
    coord[4, :] = coord_struc[12, :]
    coord[5, :] = coord_struc[9, :]
    coord[6, :] = coord_struc[13, :]
    coord[7, :] = coord_struc[8, :]
    coord[8, :] = coord_struc[14, :]
    coord[9, :] = coord_struc[7, :]
    coord[10, :] = coord_struc[15, :]
    coord[11, :] = coord_struc[6, :]
    coord[12, :] = coord_struc[16, :]
    coord[13, :] = coord_struc[5, :]
    coord[14, :] = coord_struc[17, :]
    coord[15, :] = coord_struc[4, :]
    coord[16, :] = coord_struc[18, :]
    coord[17, :] = coord_struc[3, :]
    coord[18, :] = coord_struc[19, :]
    coord[19, :] = coord_struc[20, :]
    coord[20, :] = coord_struc[2, :]
    return coord
end


function get_v3_case_params()
    wing_type = "LEI_kite"
    dist = "lin"
    N_split = 4
    aoas = collect(-4:2:22)
    Umag = 22
    # convergence criteria
    max_iterations = 1500
    allowed_error = 1e-5
    relaxation_factor = 0.03
    core_radius_fraction = 1e-20

    # Wing geometry
    coord_struc = get_CAD_matching_uri()
    coord = struct2aero_geometry(coord_struc) / 1000

    N = length(coord) ÷ 2

    # LE thickness at each section [m]
    # 10 sections
    LE_thicc = 0.1

    # Camber for each section (ct in my case)
    camber = 0.095

    # # Refine structural mesh into more panels
    # coord = thesis_functions.refine_LEI_mesh(coord, N - 1, N_split)
    # N = int(len(coord) / 2)  # Number of section after refining the mesh

    # # Definition of airfoil coefficients
    # # Based on Breukels (2011) correlation model
    # aoas_for_polar = np.arange(-80, 80, 0.1)
    # data_airf = np.empty((len(aoas_for_polar), 4))
    # for j in range(len(aoas_for_polar))
    #     alpha = aoas_for_polar[j]
    #     Cl, Cd, Cm = thesis_functions.LEI_airf_coeff(LE_thicc, camber, alpha)
    #     data_airf[j, 0] = alpha
    #     data_airf[j, 1] = Cl
    #     data_airf[j, 2] = Cd
    #     data_airf[j, 3] = Cm
    # end

    # Atot = test_utils.calculate_projected_area(coord)
    # coord_input_params = [coord, LE_thicc, camber]
    # case_parameters = [
    #     coord_input_params,
    #     aoas,
    #     wing_type,
    #     Umag,
    #     0,
    #     Atot,
    #     max_iterations,
    #     allowed_error,
    #     relaxation_factor,
    #     core_radius_fraction,
    #     data_airf,
    # ]
    #
    # return case_parameters
end