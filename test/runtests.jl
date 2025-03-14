using Test

cd("..")
println("Running tests...")
@testset verbose = true "Testing VortexStepMethod..." begin
    include("bench.jl")
    include("test_bound_filament.jl")
    include("test_panel.jl")
    include("test_semi_infinite_filament.jl")
    include("test_body_aerodynamics.jl")
    include("test_kite_geometry.jl")
    include("test_wing_geometry.jl")
    include("test_plotting.jl")
    include("aqua.jl")
end