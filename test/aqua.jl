using Aqua
@testset "Aqua.jl" begin
    Aqua.test_all(
      VortexStepMethod;
      stale_deps=(ignore=[:Xfoil, :Timers, :PyCall],),
      deps_compat=(ignore=[:PyCall],)
    )
  end
