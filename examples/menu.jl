using ControlPlots
using VortexStepMethod
using REPL.TerminalMenus

url = "https://albatross-kite-transport.github.io/VortexStepMethod.jl/dev"

options = ["rectangular_wing = include(\"rectangular_wing.jl\")",
           "ram_air_kite = include(\"ram_air_kite.jl\")",
           "lei_kite = include(\"lei_kite.jl\")",
           "bench = include(\"bench.jl\")",
           "cleanup = include(\"cleanup.jl\")",
           "help = run(`firefox $url`)",
           "quit"]

function example_menu()
    active = true
    while active
        menu = RadioMenu(options, pagesize=8)
        choice = request("\nChoose function to execute or `q` to quit: ", menu)

        if choice != -1 && choice != length(options)
            eval(Meta.parse(options[choice]))
        else
            println("Left menu. Press <ctrl><d> to quit Julia!")
            active = false
        end
    end
end

example_menu()