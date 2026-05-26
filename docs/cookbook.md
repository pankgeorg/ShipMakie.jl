# ShipMakie cookbook

A collection of worked examples for each recipe. Each example assumes
you have an active simulation with the shared bindings shown in the
preamble below.

## Preamble

```julia
using GLMakie, ShipMakie    # or CairoMakie for 2D-only
using ShipShapes: StaticArrays
const SVector = StaticArrays.SVector

# Simulation state from your sim — adapt to your variable names
NX, NY, NZ = size(vof.α)
H_w_c = NZ / 2                                  # still-water cell index
hull_xc, hull_yc, hull_zc = NX/5, NY/2, H_w_c   # hull origin in cells
L_c, B_c, T_c = 36, 8, 5                        # hull dimensions

# Per-cell hull SDF closure in world coords
hull_sdf_world(p) = ShipShapes.wigley_sdf(
    SVector(p[1] - hull_xc, p[2] - hull_yc, p[3] - hull_zc),
    Float64(L_c), Float64(B_c), Float64(T_c))

# Mask that hides the hull's plan-view footprint on η heatmaps
mask_inside(i, j) =
    abs(Float32(i - 1.5) - hull_xc) ≤ L_c / 2 &&
    abs(Float32(j - 1.5) - hull_yc) ≤ B_c / 2
```

## Quickstart — the 4-panel default scene

The fastest way to look at a simulation:

```julia
fig = ShipMakie.default_scene(sim.flow.u, vof.α;
    hull_sdf    = hull_sdf_world,
    waterline_z = H_w_c,
    mask        = mask_inside,
    hull_box    = (hull_xc, hull_yc, L_c, B_c),
    rotor_center = (prop_xc, prop_yc, prop_zc),
    rotor_R     = R_prop, rotor_R_hub = 0.5,
)
save("scene.png", fig)
```

## 2D recipes

### `etaheatmap` — free-surface elevation

Reads the VoF colour function `α`, finds the cell-by-cell water-line
crossing, and plots η(x, y).

```julia
fig = Figure(size = (800, 400))
ax  = Axis(fig[1, 1]; aspect = DataAspect())
etaheatmap!(ax, vof.α;
    waterline_z = H_w_c,           # midplane is the still-water surface
    mask        = mask_inside,     # hide the hull plan-view footprint
    hull_box    = (hull_xc, hull_yc, L_c, B_c),
    colorrange  = (-0.4, 0.4))     # symmetric around 0
```

Pro tip: for animations, fix `colorrange` (don't let it auto-scale)
so frame-to-frame η differences are visible.

### `velocityslice` — velocity component on a slice

```julia
ax = Axis(fig[1, 1]; aspect = DataAspect())
velocityslice!(ax, sim.flow.u;
    slice_axis = :y,           # slice perpendicular to y
    index      = NY ÷ 2,       # midplane
    component  = 1,            # u_x
    colorrange = (-1.5, 3.0))
```

### `streamslice` — 2D streamlines

```julia
streamslice!(ax, sim.flow.u;
    slice_axis = :y, index = NY ÷ 2,
    density    = 0.8,
    linewidth  = 1.0)
```

### `vectorslice` — 2D velocity arrows

```julia
vectorslice!(ax, sim.flow.u;
    slice_axis = :y, stride = 4,
    lengthscale = 1.0, arrowsize = 8)
```

Use `stride` to avoid clutter on dense grids.

### `hullsilhouette` — body outline on a 2D slice

Pairs with `velocityslice` or `streamslice` for "where's the hull
in this slice?" overlay.

```julia
hullsilhouette!(ax, hull_sdf_world;
    grid_size = (NX, NY, NZ),
    slice_axis = :y, slice_index = NY ÷ 2,
    color = :black, linewidth = 2)
```

### `probeline` — 1D field profile

```julia
# Centreline z-profile of α at midship
probeline!(ax, vof.α;
    along = :z,
    fixed_indices = (round(Int, hull_xc), round(Int, hull_yc)))
```

`along ∈ {:x, :y, :z}`. For 3D fields, supply two `fixed_indices`;
for 2D, supply one.

## 3D recipes (need GLMakie or WGLMakie)

### `hullmesh3d` — hull surface from SDF

```julia
ax = Axis3(fig[1, 1]; aspect = :data)
hullmesh3d!(ax, hull_sdf_world;
    grid_size = (NX, NY, NZ),
    color = :navy)
```

Use `transparency = true` for see-through hulls when stacking with
other 3D plots.

### `freesurface3d` — water surface mesh

```julia
freesurface3d!(ax, vof.α;
    waterline_z = H_w_c,
    mask        = mask_inside,
    colorrange  = (-0.4, 0.4),
    alpha       = 0.55)
```

### `vorticityvolume` — vorticity volume rendering

```julia
vorticityvolume!(ax, sim.flow.u;
    field     = :omega_mag,    # or :lambda2 for Jeong–Hussain
    algorithm = :mip,          # or :iso (with isovalue=...)
    colormap  = :algae,
    absorption = 4.0)
```

For sharp vortex tubes use `field = :lambda2, algorithm = :iso`.
For a smooth volumetric look use `:omega_mag, :mip`.

### `bladedrotor3d` — rotor wireframe

```julia
bladedrotor3d!(ax, 3, R_prop, 0.5 * R_prop;
    center     = (prop_xc, prop_yc, prop_zc),
    rotor_axis = (1.0, 0.0, 0.0),    # axis = +x (standard)
    color      = :orange,
    linewidth  = 3)
```

### `rudder3d` — rudder plate with δ rotation

```julia
rudder3d!(ax, rudder_chord, rudder_span;
    center     = (rud_xc, rud_yc, rud_zc),
    rudder_axis = (0.0, 0.0, 1.0),    # span = +z
    δ          = deg2rad(8.0),
    color      = :lime,
    linewidth  = 3)
```

### `pressureisosurface` — pressure iso surfaces

```julia
pressureisosurface!(ax, sim.flow.p;
    level = [-0.5, 0.5],           # explicit levels; or automatic
    color = [:tomato, :steelblue],
    transparency = true)
```

### `streamlines3d` — 3D streamlines from seeds

```julia
# Seed in a transverse plane upstream of the bow
seeds = NTuple{3, Float64}[]
for _ in 1:80
    y = hull_yc - B_c + 2*B_c*rand()
    z = hull_zc - T_c + 2*T_c*rand()
    push!(seeds, (Float64(hull_xc - L_c/2 - 3), y, z))
end
streamlines3d!(ax, sim.flow.u, seeds;
    nsteps = 250, dt = 0.4,
    colormap = :plasma, linewidth = 2.0)
```

`ShipMakie.random_seeds(grid_size, n)` is a convenience for uniform
scatter; for purposeful seeding (like the bow plane above), construct
the vector manually.

## Layering recipes

The standard "everything in one Axis3" recipe stack:

```julia
fig = Figure(size = (1400, 900))
ax = Axis3(fig[1, 1]; aspect = :data,
    azimuth = 0.40π, elevation = 0.20π)

hullmesh3d!(ax, hull_sdf_world; grid_size = (NX, NY, NZ), color = :navy)
freesurface3d!(ax, vof.α; alpha = 0.55, colorrange = (-0.4, 0.4))
vorticityvolume!(ax, sim.flow.u; field = :omega_mag)
bladedrotor3d!(ax, 3, R_prop, 0.5 * R_prop; center = (prop_xc, prop_yc, prop_zc))
rudder3d!(ax, 4.0, 5.0; center = (rud_xc, rud_yc, rud_zc), δ = deg2rad(8))
```

Order matters for transparency: render opaque things (`hullmesh3d`)
first, then translucent (`freesurface3d`), then volumes
(`vorticityvolume`).

## Headless 3D rendering (GLMakie + xvfb-run)

On a Linux box without a display:

```bash
xvfb-run -a julia my_script.jl
```

This wraps the script in a virtual X server. GLMakie loads, renders
to an off-screen framebuffer, and `save("foo.png", fig)` produces a
correct image.

## Animation pattern

```julia
GLMakie.activate!()
fig = Figure(); ax = Axis3(fig[1, 1])

# Set up Observables — Makie watches these and updates the plot
α_obs = Observable(vof.α)
u_obs = Observable(sim.flow.u)

etaheatmap!(Axis(fig[2, 1]), α_obs; waterline_z = H_w_c)
vorticityvolume!(ax, u_obs; field = :omega_mag)

record(fig, "out.gif", 1:60; framerate = 24) do frame
    WaterLily.mom_step!(sim.flow, sim.pois)
    step_vof_mules!(vof, sim; dt = sim.flow.Δt[end-1])
    α_obs[] = vof.α     # trigger reactive update
    u_obs[] = sim.flow.u
end
```

Or skip Observables and save per-frame PNGs (faster for non-interactive
renders):

```julia
for frame in 1:60
    ...
    save("frames/$(lpad(frame, 5, '0')).png", fig)
end
run(`convert -delay 4 frames/*.png out.gif`)   # ImageMagick
```

The example scripts in `examples/` use the latter pattern.

## Performance notes

- `vorticityvolume` recomputes ω/λ₂ on every frame. For a static
  snapshot this is fine; for animation, the cost is ~`O(NX·NY·NZ)`
  per frame which is small relative to a WaterLily step.
- `hullmesh3d` samples the SDF on the grid each frame. If the hull
  is fixed, compute the SDF array once and pass it as the input
  instead of the closure.
- `streamlines3d` integrates `nsteps` per seed at `dt`. For 80 seeds
  × 250 steps that's 20 000 RK4 evaluations; on the order of
  10–100 ms per frame, modest compared to plotting.

## Where to go from here

- `examples/` — every recipe demonstrated in isolation and combined.
- `test/runtests.jl` — minimal "does this recipe build a Figure?"
  tests; copy these as starter templates.
- `src/utils.jl` — `_axis_idx`, `cell_centres`, `eta_from_alpha`
  helpers reusable in custom recipes.
