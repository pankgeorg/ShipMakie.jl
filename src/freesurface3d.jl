"""
    freesurface3d(α; waterline_z = size(α, 3) / 2)

Plot the free surface in 3D as a height-mapped surface, using the
η(x, y) elevation extracted from the VoF colour function `α`. The
height channel is set to `η + waterline_z` (world z), so the surface
appears in absolute coordinates and lines up with hull / rotor
features.

# Attributes

- `waterline_z` — z cell index of the still-water surface (default
  `size(α, 3) / 2`).
- `mask` — optional `(i, j) -> Bool` to NaN-out columns inside the
  hull silhouette.
- `colormap` — default `:RdBu`.
- `colorrange` — default `automatic`; symmetric scaling around 0.
- `alpha` — surface transparency (0–1), default 0.85.

# Example

```julia
fig = Figure(); ax = Axis3(fig[1, 1])
freesurface3d!(ax, vof.α; waterline_z = NZ/2)
```
"""
Makie.@recipe(FreeSurface3d, α) do scene
    Makie.Attributes(
        waterline_z = Makie.automatic,
        mask        = nothing,
        colormap    = :RdBu,
        colorrange  = Makie.automatic,
        alpha       = 0.85,
    )
end

function Makie.plot!(p::FreeSurface3d)
    α_obs = p[1]
    # η(x, y) elevation + absolute world-frame surface z
    η_obs = Makie.lift(α_obs, p.waterline_z, p.mask) do α, wl_z, mask
        nx, ny, nz = size(α)
        wl = wl_z === Makie.automatic ? nz / 2 : wl_z
        eta_from_alpha(α, wl;
            mask_inside = mask isa Function ? mask : nothing)
    end
    z_obs = Makie.lift(η_obs, α_obs, p.waterline_z) do η, α, wl_z
        nz = size(α, 3)
        wl = wl_z === Makie.automatic ? nz / 2 : wl_z
        η .+ wl
    end
    xs_obs = Makie.lift(α_obs) do α; cell_centres(size(α, 1)); end
    ys_obs = Makie.lift(α_obs) do α; cell_centres(size(α, 2)); end
    Makie.surface!(p, xs_obs, ys_obs, z_obs;
        color      = η_obs,    # use η as the colour channel
        colormap   = p.colormap,
        colorrange = p.colorrange,
        alpha      = p.alpha,
    )
    return p
end
