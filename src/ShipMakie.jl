"""
    ShipMakie

Makie recipes for the foam ship-CFD stack. Provides reusable plot
types for free-surface elevation, velocity slices, vorticity volumes,
hull silhouettes, and rotor / rudder geometries. The package depends
only on `Makie` so it can be loaded alongside any Makie backend
(CairoMakie, GLMakie, WGLMakie) without pulling in the heavy plotting
stack.

The recipes here consume raw simulation outputs from the foam stack:
- WaterLily's `flow.u`, `flow.p`
- VoF's `vof.α`
- ShipShapes' `Wigley`/`Containership` SDFs
- LiftingSurfaces' `BladedRotor`, `Rudder`

Each recipe is a thin shim that converts the simulation's
cell-units convention (face-staggered grids, 0.5-cell ghost offset)
into Makie-friendly inputs.

# Available recipes

| Recipe              | What it draws                                  |
|---------------------|------------------------------------------------|
| `etaheatmap`        | 2D free-surface elevation η(x, y) from α       |
| `velocityslice`     | 2D velocity-component heatmap on a coord slice |
| `vorticityvolume`   | 3D vorticity-magnitude volume (iso or MIP)     |
| `hullsilhouette`    | 2D body silhouette derived from an SDF         |
| `probeline`         | 1D profile along an x/y/z line                 |
| `bladedrotor3d`     | 3D wireframe of a BladedRotor                  |
"""
module ShipMakie

using Makie
using StaticArrays

# Helpers (private) — loaded first so recipes can use them
include("utils.jl")

# Public recipes — each lives in its own file
include("etaheatmap.jl")
include("velocityslice.jl")
include("hullsilhouette.jl")
include("probeline.jl")
include("vorticityvolume.jl")
include("bladedrotor3d.jl")

# Curated colormaps for ship-CFD work
include("colormaps.jl")

export etaheatmap, etaheatmap!,
       velocityslice, velocityslice!,
       hullsilhouette, hullsilhouette!,
       probeline, probeline!,
       vorticityvolume, vorticityvolume!,
       bladedrotor3d, bladedrotor3d!

end # module
