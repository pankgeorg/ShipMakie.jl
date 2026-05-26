"""
    default_scene(u, α; hull_sdf = nothing, …)

Return a 4-panel `Figure` with the standard ship-CFD view:
- Top-left:  η free-surface heatmap from `α` (plan view)
- Top-right: u_x slice (side view) with optional hull silhouette
- Bottom-left: streamslice on the centreline plane
- Bottom-right: 3D scene with hullmesh + freesurface + vortex volume

This is the "I just want to look at my simulation" preset. All the
attributes are exposed as keyword arguments. Returns the `Figure` so
the caller can save / display it.

# Keyword arguments

- `hull_sdf` — optional callable `p -> Real` returning the body SDF
  in world coords. If supplied, hull silhouettes / mesh are drawn.
- `waterline_z` — z cell index of still-water surface; default
  `size(α, 3) / 2`.
- `hull_box` — optional `(xc, yc, L, B)` plan-view bounding box.
- `mask` — optional `(i, j) -> Bool` to NaN-out columns inside the
  hull on the η heatmap.
- `slice_index` — y cell index for the side slice; default mid.
- `rotor_center`, `rotor_R`, `rotor_R_hub`, `rotor_N` — if supplied,
  draws the rotor wireframe on top of the 3D scene.
- `figure_size` — passed to `Figure`. Default `(1400, 900)`.

# Example

```julia
using GLMakie, ShipMakie
fig = ShipMakie.default_scene(sim.flow.u, vof.α;
    hull_sdf = p -> wigley_sdf(p, L, B, T),
    waterline_z = NZ/2,
    hull_box = (hull_xc, hull_yc, L_c, B_c),
    mask = (i, j) -> ... ,
    rotor_center = (prop_xc, prop_yc, prop_zc),
    rotor_R = R_prop, rotor_R_hub = 0.5,
)
save("scene.png", fig)
```
"""
function default_scene(u::AbstractArray{T, 4}, α::AbstractArray{S, 3};
        hull_sdf      = nothing,
        waterline_z   = nothing,
        hull_box      = nothing,
        mask          = nothing,
        slice_index   = nothing,
        rotor_center  = nothing,
        rotor_R       = 2.0,
        rotor_R_hub   = 0.5,
        rotor_N       = 3,
        figure_size   = (1400, 900),
        omega_field   = :omega_mag,
        omega_alg     = :mip,
        omega_cmap    = :algae,
    ) where {T, S}
    nx, ny, nz = size(α)
    wl_z = waterline_z === nothing ? nz / 2 : waterline_z
    sl_y = slice_index === nothing ? ny ÷ 2 : slice_index

    fig = Makie.Figure(size = figure_size, backgroundcolor = :white)

    # Top-left: η heatmap
    ax1 = Makie.Axis(fig[1, 1]; aspect = Makie.DataAspect(),
        title = "η free surface (plan view)",
        xlabel = "x (cells)", ylabel = "y (cells)")
    etaheatmap!(ax1, α;
        waterline_z = wl_z, mask = mask, hull_box = hull_box,
        colorrange = (-0.4, 0.4))

    # Top-right: u_x side slice
    ax2 = Makie.Axis(fig[1, 2]; aspect = Makie.DataAspect(),
        title = "u_x centreline (side)",
        xlabel = "x (cells)", ylabel = "z (cells)")
    velocityslice!(ax2, u;
        slice_axis = :y, index = sl_y, component = 1,
        colorrange = (-1.5, 3.0))
    if hull_sdf !== nothing
        hullsilhouette!(ax2, hull_sdf;
            grid_size = (nx, ny, nz), slice_axis = :y, slice_index = sl_y)
    end

    # Bottom-left: streamslice
    ax3 = Makie.Axis(fig[2, 1]; aspect = Makie.DataAspect(),
        title = "streamslice (centreline)",
        xlabel = "x (cells)", ylabel = "z (cells)")
    streamslice!(ax3, u;
        slice_axis = :y, index = sl_y, density = 0.7)
    if hull_sdf !== nothing
        hullsilhouette!(ax3, hull_sdf;
            grid_size = (nx, ny, nz), slice_axis = :y, slice_index = sl_y)
    end

    # Bottom-right: 3D scene
    ax4 = Makie.Axis3(fig[2, 2];
        title = "3D scene",
        aspect = :data,
        azimuth = 0.40π, elevation = 0.20π)
    if hull_sdf !== nothing
        hullmesh3d!(ax4, hull_sdf;
            grid_size = (nx, ny, nz), color = :navy)
    end
    freesurface3d!(ax4, α;
        waterline_z = wl_z, mask = mask,
        colorrange = (-0.4, 0.4), alpha = 0.55)
    vorticityvolume!(ax4, u;
        field = omega_field, algorithm = omega_alg,
        colormap = omega_cmap, absorption = 4.0)
    if rotor_center !== nothing
        bladedrotor3d!(ax4, rotor_N, rotor_R, rotor_R_hub;
            center = rotor_center, color = :orange, linewidth = 3)
    end

    return fig
end

"""
    record_default_scene(filename::String, u_obs, α_obs, step!;
                          nframes = 60, framerate = 20,
                          hull_sdf = nothing, …)

Convenience wrapper that records an animation of the
`default_scene` view to `filename` (mp4 / gif / webm — backend
decides). Uses Makie's `record` with Observable inputs so per-frame
plot updates are reactive.

# Arguments

- `filename` — output file path (extension determines format).
- `u_obs::Observable{Array}` — wraps the live velocity field.
- `α_obs::Observable{Array}` — wraps the live colour function.
- `step!()` — a no-arg function that advances the simulation one
  step (does the WaterLily / VoF update and any other bookkeeping).
  After it returns, the caller is responsible for updating `u_obs[]`
  and `α_obs[]` to the new arrays.

# Keywords

- `nframes`, `framerate` — animation parameters.
- All other kwargs are forwarded to `default_scene`.

# Example

```julia
u_obs = Observable(sim.flow.u)
α_obs = Observable(vof.α)
function step!()
    WaterLily.mom_step!(sim.flow, sim.pois; udf=udf, pois_tol=1f-6)
    step_vof_mules!(vof, sim; dt=sim.flow.Δt[end-1])
    u_obs[] = sim.flow.u
    α_obs[] = vof.α
end
ShipMakie.record_default_scene("scene.mp4", u_obs, α_obs, step!;
    nframes = 60, framerate = 24,
    hull_sdf = my_hull_sdf,
    rotor_center = (px, py, pz))
```
"""
function record_default_scene(filename::AbstractString,
        u_obs::Makie.Observable, α_obs::Makie.Observable, step_fn::Function;
        nframes::Int = 60, framerate::Int = 20,
        kwargs...)
    fig = default_scene(u_obs[], α_obs[]; kwargs...)
    Makie.record(fig, filename, 1:nframes; framerate = framerate) do _frame
        step_fn()
        # Trigger reactive updates (caller mutated the wrapped arrays)
        Makie.notify(u_obs)
        Makie.notify(α_obs)
    end
    return filename
end

