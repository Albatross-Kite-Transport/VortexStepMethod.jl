## Unreleased
### Added
- The function `install_examples()` which allows to easily install the examples without using `git`
- The function `solve!` returns a struct now. The function `solve`that returns a dict is still available.
- The moment coefficients distribution in `solve!`
- The script `install` to the `bin` folder for users who checked out this git repository
- The script `bench2.jl` was added for allocation testing of the `solve!` function
### Changed
- `deform!` by a distribution instead of just a left and right angle
- Read the y-coordinates in the correct direction from the `ram_air_kite_body.obj` file
- In the `menu.jl`, changed `help` to `help_me`. It works better now, no more warnings on Linux, it should also work on MacOS now
- The coordinate frames of the panels now use the same convention as the kite body frame
- The page "Glossary" of the documentation is quite complete now
- The center of mass field of the `RamAirWing` is removed, and the geometry is created such that `[0, 0, 0]` is the center of mass
- The enumeration `WingType` was added and replaces the symbols, used before
- The allocations of the function `solve!` where reduced by a factor of three
### Fixed
- The function `calculate_circulation_distribution_elliptical_wing()` was never called
- Fix the calculation of force coefficients in `solve!`
- The continues integration scripts (CI.yml) use now separate runs for the test coverage and for the allocation tests.

## VortexStepMethod v1.1.0 2025-03-04
### Added
- Dynamically deform the RamAirWing by twisting the left side and right side, and deforming the trailing edges using deform! #19
- Set turn rate `omega = [omega_x, omega_y, omega_z]` in kite body frame using set_va! #49
- Add moment coefficient calculations around specified point to `solve` #17
- Add moment distribution of the moment around the local panel y-axes around user-defined points on the panels to `solve!` #90
- Add function `solve!()` which returns a `VSMSolution` struct #87
- Add the option to remove the NaNs in `aero_data` vectors or matrices using the `remove_nan` keyword in the `Wing` and `RamAirWing` constructors #98
### Changed
- Add origin argument to `BodyAerodynamics` constructor #66
- Improve documentation

## Initial Release
This project is based on version 1.0 of the Python project [Vortex-Step-Method](https://github.com/ocayon/Vortex-Step-Method)

## Noteworthy Differences of v1.0.1 to the Python version
- implemented in Julia, therefore about 50 times faster
- an importer for `.obj` wing geometry files was added (see: [#10](https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl/issues/10))
- automatic creation of polars using Xfoil.jl was added (see: [#43](https://github.com/Albatross-Kite-Transport/VortexStepMethod.jl/pull/43))
- a ram-air kite example was added
- `Umag` was replaced with `v_a` as variable name for the norm of the apparent wind speed
- memory allocations were significantly reduced
- a menu (examples/menu.jl) for running the examples was added
- plotting was moved to an extension #55
- added improved online documentation
