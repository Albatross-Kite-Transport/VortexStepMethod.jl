[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://Albatross-Kite-Transport.github.io/VortexStepMethod.jl/dev) [![Build Status](https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Coverage](https://codecov.io/gh/Albatross-Kite-Transport/VortexStepMethod.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Albatross-Kite-Transport/VortexStepMethod.jl)


# Aerodynamic models of 3D wings using the Vortex Step Method 

The Vortex Step Method (VSM) is an enhanced lifting line method that improves upon the classic approach by solving the circulation system at the three-quarter chord position, among the most important details. This adjustment allows for more accurate calculations of lift and drag forces, particularly addressing the shortcomings in induced drag prediction. 
VSM is further refined by coupling it with 2D viscous airfoil polars, making it well-suited for complex geometries, 
including low aspect ratio wings, as well as configurations with sweep, dihedral, and anhedral angles.

The software presented here includes a couple of examples: a rectangular wing, a leading-edge inflatable kite and a ram-air kite.

This package was translated from the Python code version 1.0.0 available at https://github.com/ocayon/Vortex-Step-Method with some extensions as documented in [News.md](https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl/blob/main/NEWS.md).

## Installation
Install [Julia 1.10](https://ufechner7.github.io/2024/08/09/installing-julia-with-juliaup.html) or later, 
if you haven't already. On Linux, make sure that Python3, Matplotlib and LaTeX are installed:
```
sudo apt install python3-matplotlib
sudo apt install texlive-full texlive-fonts-extra cm-super
```
Furthermore, the packages `TestEnv` and `ControlPlots` must be installed globally:
```
julia -e 'using Pkg; Pkg.add("TestEnv"); Pkg.add(ControlPlots)'
```)

Before installing this software it is suggested to create a new project, for example like this:
```bash
mkdir test
cd test
julia --project=.
```
Then add VortexStepMethod from  Julia's package manager, by typing:
```julia
using Pkg
pkg"add https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl"
``` 
at the Julia prompt. You can run the unit tests with the command:
```julia
pkg"test VortexStepMethod"
```

## Running the examples
If you have git installed, check out this repo because it makes it easier to understand the code:
```bash
mkdir repos
cd repos
git clone https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl
cd VortexStepMethod.jl
```
You can launch Julia with:
```bash
julia --project
```
or with:
```bash
./bin/run_julia
```
In Julia, first update the packages:
```julia
using Pkg
Pkg.update()
```
and then you can display a menu with the available examples:
```julia
using ControlPlots
using VortexStepMethod
menu()
```
To browse the code, it is suggested to use [VSCode](https://code.visualstudio.com/) with the Julia plugin.

## Input
Three kinds of input data is needed:

- The wing geometry, defined by section:
  - for the rectangualar wing two sections, two points in CAD reference frame + polars  
    (three different options to provide them) per section
  - kite wing: model of polars included, n sections to define

- The airflow and turn rate:
  - `v_app` vector and `omega` (turn rate) vector in Kite Body (KB) reference frame

- The configuration:
  - how many panels  
    --> two sections make a panel.

Apart from the wing geometry there is no input file yet, the input has to be defined in the code.

### Example for defining the required input:
```julia

# Step 1: Define wing parameters
n_panels = 20          # Number of panels
span = 20.0            # Wing span [m]
chord = 1.0            # Chord length [m]
v_a = 20.0             # Magnitude of inflow velocity [m/s]
density = 1.225        # Air density [kg/m³]
alpha_deg = 30.0       # Angle of attack [degrees]
alpha = deg2rad(alpha_deg)

# Step 2: Create wing geometry with linear panel distribution
wing = Wing(n_panels, spanwise_panel_distribution=LINEAR)

# Add wing sections - defining only tip sections with inviscid airfoil model
add_section!(wing, 
    [0.0, span/2, 0.0],    # Left tip LE 
    [chord, span/2, 0.0],  # Left tip TE
    INVISCID)
add_section!(wing, 
    [0.0, -span/2, 0.0],   # Right tip LE
    [chord, -span/2, 0.0], # Right tip TE
    INVISCID)

# Step 3: Initialize aerodynamics
wa = BodyAerodynamics([wing])

# Set inflow conditions
vel_app = [cos(alpha), 0.0, sin(alpha)] .* v_a
set_va!(wa, vel_app)
```
It is possible to import the wing geometry using an `.obj` file as shown in the example `ram_air_kite.jl`.

Surfplan files can be converted to an input for `VortexStepMethod.jl` using the [SurfplanAdapter](https://github.com/jellepoland/SurfplanAdapter).

## Output
- CL, CD, CS (side force coefficient)
- the spanwise distribution of forces  
  --> moment coefficients (will be implemented in release 1.1) 

## Citation
If you use this project in your research, please consider citing it. 
Citation details can be found in the [CITATION.cff](https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl/blob/main/CITATION.cff) file included in this repository.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## WAIVER
Technische Universiteit Delft hereby disclaims all copyright interest in the package “VortexStepMethod.jl” written by the Author(s).

Prof.dr. H.G.C. (Henri) Werij, Dean of Aerospace Engineering

### Copyright
Copyright (c) 2022 Oriol Cayon

Copyright (c) 2024 Oriol Cayon, Jelle Poland, TU Delft

Copyright (c) 2025 Oriol Cayon, Jelle Poland, Bart van de Lint, Uwe Fechner
