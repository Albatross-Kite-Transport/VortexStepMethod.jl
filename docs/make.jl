using Pkg
if ("TestEnv" ∈ keys(Pkg.project().dependencies))
    if ! ("Documents" ∈ keys(Pkg.project().dependencies))
        using TestEnv; TestEnv.activate()
    end
end
using ControlPlots
using VortexStepMethod
using Documenter

DocMeta.setdocmeta!(VortexStepMethod, :DocTestSetup, :(using VortexStepMethod); recursive=true)

makedocs(;
    modules=[VortexStepMethod,
             isdefined(Base, :get_extension) ? 
             Base.get_extension(VortexStepMethod, :VortexStepMethodControlPlotsExt) :
             VortexStepMethod.VortexStepMethodControlPlotsExt],
    authors="Uwe Fechner <uwe.fechner.msc@gmail.com>, Bart van de Lint <bart@vandelint.net> and contributors",
    sitename="VortexStepMethod.jl",
    checkdocs=:none,
    format = Documenter.HTML(prettyurls = haskey(ENV, "CI")),
    pages=[
        "Home" => "index.md",
        "How it works" => "explanation.md",
        "Exported Functions" => "functions.md",
        "Exported Types" => "types.md",
        "Private Functions" => "private_functions.md",
        "Private Types" => "private_types.md",
        "Reference Frames" => "reference_frames.md"
    ],
)

deploydocs(;
    repo="github.com/Albatross-Kite-Transport/VortexStepMethod.jl",
    devbranch="main",
)
