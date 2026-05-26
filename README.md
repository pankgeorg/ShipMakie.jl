# ShipMakie.jl

Makie recipes for the foam ship-CFD stack. Reusable plot types for
WaterLily / VoF / Turbulence / LiftingSurfaces / ShipShapes outputs.

## Recipes

| Recipe              | What it draws                                  |
|---------------------|------------------------------------------------|
| `etaheatmap`        | 2D free-surface elevation η(x, y) from α       |
| `velocityslice`     | 2D velocity-component heatmap on a coord slice |
| `vorticityvolume`   | 3D ω_mag or λ₂ isosurface / MIP                |
| `hullsilhouette`    | 2D body silhouette derived from an SDF         |
| `probeline`         | 1D profile along x/y/z                         |
| `bladedrotor3d`     | 3D wireframe of a BladedRotor                  |

All recipes follow Makie conventions: `etaheatmap(α; …)` creates a new
figure/axis, `etaheatmap!(ax, α; …)` plots into an existing axis. The
input observables are reactive — passing an `Observable{Array}` makes
the plot update when the array changes (use this for animations).

## Quick start

```julia
using CairoMakie, ShipMakie

# 1. free-surface elevation
fig, ax, p = etaheatmap(vof.α; waterline_z = NZ/2)

# 2. centreline velocity slice
ax2 = Axis(fig[1, 2])
velocityslice!(ax2, sim.flow.u; slice_axis = :y, component = 1)

# 3. hull silhouette overlay
hullsilhouette!(ax2, p -> wigley_sdf(p, L, B, T);
    grid_size = size(vof.α), slice_axis = :y)

# 4. vortex volume (needs GLMakie for 3D)
using GLMakie
fig3 = Figure()
ax3 = Axis3(fig3[1, 1])
vorticityvolume!(ax3, sim.flow.u; field = :lambda2, algorithm = :iso)

# 5. rotor wireframe overlay
bladedrotor3d!(ax3, 3, 2.4, 0.5; center = (prop_xc, prop_yc, prop_zc))
```

## Design notes

- **Only depends on `Makie`** — backend choice (CairoMakie / GLMakie /
  WGLMakie) is left to the downstream user.
- **No simulation logic**. The recipes are thin shims over WaterLily's
  face-staggered velocity convention. Use them as the visualisation
  layer for any cell-centred CFD output that matches the layout.
- **Observable-friendly**. All input arrays go through `lift(...)` so
  the plots respond to time evolution without per-frame rebuild cost.
- **Cell-units convention**. Cell `I` has centre at world `I - 1.5`
  (WaterLily ghost convention). Cell-centre coords are emitted on the
  output axes so heatmaps overlay cleanly.

## Status

Tested with CairoMakie. 19 tests across 8 testsets pass. See
`test/runtests.jl` for the canonical examples.

3D recipes (`vorticityvolume`, `bladedrotor3d`) render with CairoMakie
but the 3D volume rendering quality is much better with `GLMakie`.
