
using Test
using VortexStepMethod
using VortexStepMethod: create_interpolations, find_circle_center_and_radius, calculate_inertia_tensor, calculate_com, read_faces
using LinearAlgebra
using Interpolations
using Serialization

@testset "Kite Geometry Tests" begin
    # Test data
    test_obj_path = joinpath(tempdir(), "test.obj")
    test_dat_path = joinpath(tempdir(), "test.dat")
    
    @testset "OBJ File Reading" begin
        # Create minimal test OBJ file
        test_vertices = """
        v 0.0 0.0 0.0
        v 1.0 0.0 0.0
        v 0.0 1.0 0.0
        f 1 2 3
        """
        write(test_obj_path, test_vertices)
        
        vertices, faces = read_faces(test_obj_path)
        
        @test length(vertices) == 3
        @test length(faces) == 1
        @test vertices[1] ≈ [0.0, 0.0, 0.0]
        @test vertices[2] ≈ [1.0, 0.0, 0.0]
        @test vertices[3] ≈ [0.0, 1.0, 0.0]
        @test faces[1] == [1, 2, 3]
    end
    
    @testset "Center of Mass Calculation" begin
        vertices = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
        faces = [[1, 2, 3]]
        
        com = calculate_com(vertices, faces)
        expected_com = [1/3, 1/3, 0.0]
        
        @test isapprox(com, expected_com, rtol=1e-5)
    end
    
    @testset "Inertia Tensor Calculation" begin
        vertices = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
        faces = [[1, 2, 3]]
        mass = 1.0
        com = [1/3, 1/3, 0.0]
        
        I = calculate_inertia_tensor(vertices, faces, mass, com)
        
        # Test properties of inertia tensor
        @test size(I) == (3,3)
        @test isapprox(I, I', rtol=1e-10)  # Symmetric
        @test all(diag(I) .≥ 0)  # Non-negative diagonal
    end
    
    @testset "Circle Fitting" begin
        # Create simple curved wing vertices
        r = 5.0
        z_center = 2.0
        vertices = []
        for θ in range(-π/4, π/4, length=100)
            push!(vertices, [0.0, r*sin(θ), z_center + r*cos(θ)])
        end
        
        z, radius, gamma_tip = find_circle_center_and_radius(vertices)
        
        @test isapprox(z, z_center, rtol=1e-2)
        @test isapprox(radius, r, rtol=1e-2)
        @test gamma_tip ≈ π/4 rtol=1e-2
    end
    
    r = 5.0
    @testset "Interpolation Creation" begin
        vertices = []
        z_center = 2.0
        for θ in range(-π/4, π/4, length=10)
            push!(vertices, [0.0, r*sin(θ), z_center + r*cos(θ)])
            push!(vertices, [1.0, r*sin(θ), z_center + r*cos(θ)])
        end
        
        # Create test airfoil data file
        test_dat_path = joinpath(tempdir(), "test.dat")
        write(test_dat_path, "1.0 0.0\n0.0 0.0\n-1.0 0.0\n")
        
        # Create polar data
        alphas = -1.0:1.0:1.0
        d_trailing_edge_angles = -1.0:1.0:1.0
        cl_matrix = zeros(length(alphas), length(d_trailing_edge_angles))
        cd_matrix = zeros(length(alphas), length(d_trailing_edge_angles))
        cm_matrix = zeros(length(alphas), length(d_trailing_edge_angles))
        
        # Fill matrices with sample data
        for i in eachindex(alphas)
            for j in eachindex(d_trailing_edge_angles)
                cl_matrix[i,j] = sin(deg2rad(alphas[i]))
                cd_matrix[i,j] = 0.01 + 0.1*sin(deg2rad(alphas[i]))^2
                cm_matrix[i,j] = -0.1*sin(deg2rad(alphas[i]))
            end
        end
        
        # Serialize polar data
        polar_path = test_dat_path[1:end-4] * "_polar.bin"
        serialize(polar_path, (alphas, d_trailing_edge_angles, cl_matrix, cd_matrix, cm_matrix))
        
        # Create and serialize obj file
        faces = [[i, i+1, i+2] for i in 1:2:length(vertices)-2]
        open(test_obj_path, "w") do io
            for v in vertices
                println(io, "v $(v[1]) $(v[2]) $(v[3])")
            end
            for f in faces
                println(io, "f $(f[1]) $(f[2]) $(f[3])")
            end
        end
        
        # Create info file
        info_path = test_obj_path[1:end-4] * "_info.bin"
        le_interp, te_interp, area_interp = create_interpolations(vertices, z_center, r, π/4)
        center_of_mass = calculate_com(vertices, faces)
        inertia_tensor = calculate_inertia_tensor(vertices, faces, 1.0, center_of_mass)
        
        serialize(info_path, (
            center_of_mass, 
            inertia_tensor, 
            z_center, 
            r, 
            π/4, 
            le_interp, 
            te_interp, 
            area_interp
        ))
        
        # Test interpolation at middle point
        @test le_interp(0.0) ≈ 0.0 rtol=1e-2
        @test te_interp(0.0) ≈ 1.0 rtol=1e-2
    end
    
    @testset "KiteWing Construction" begin
        wing = KiteWing(test_obj_path, test_dat_path)
        
        @test wing.n_panels == 54  # Default value
        @test wing.spanwise_panel_distribution == UNCHANGED
        @test wing.spanwise_direction ≈ [0.0, 1.0, 0.0]
        @test length(wing.sections) > 0  # Should have sections now
        @test wing.mass ≈ 1.0
        @test length(wing.center_of_mass) == 3
        @test wing.radius ≈ r rtol=1e-2
        @test wing.gamma_tip ≈ π/4 rtol=1e-2
    end
    
    rm(test_obj_path)
    rm(test_dat_path)
end