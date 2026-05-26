"""
    etaheatmap(α; waterline_z = size(α, 3) / 2,
                 mask = nothing,
                 hull_box = nothing)

Plot the free-surface elevation `η(x, y)` extracted from a VoF colour
function `α::Array{T, 3}` as a 2D heatmap. The elevation is computed
relative to `waterline_z` (the z index where the still-water surface
sits; defaults to mid-domain).

# Attributes

- `waterline_z` — the z cell index of the still-water surface. Cells
  with `α ≥ 0.5` below this and `α < 0.5` above contribute to the
  surface estimate.
- `mask` — optional function `(i, j) -> Bool` returning `true` for
  cells that should be masked (rendered as NaN). Useful to hide the
  hull's plan-view footprint.
- `hull_box` — optional `(xc, yc, L, B)` to draw a rectangular hull
  outline on top of the heatmap. Set to `nothing` to skip.
- `colormap` — default `:RdBu` (symmetric around 0).
- `colorrange` — default `automatic`; set to `(-η_max, η_max)` to
  enforce symmetric scaling across frames.

# Example

```julia
using GLMakie, ShipMakie
fig, ax, p = etaheatmap(sim_α; waterline_z = 16)
```
"""
Makie.@recipe(EtaHeatmap, α) do scene
    Makie.Attributes(
        waterline_z = Makie.automatic,
        mask        = nothing,
        hull_box    = nothing,
        colormap    = :RdBu,
        colorrange  = Makie.automatic,
        nan_color   = :grey80,
        highclip    = :darkblue,
        lowclip     = :darkred,
    )
end

function Makie.plot!(p::EtaHeatmap)
    α_obs = p[1]   # the input α field as an Observable
    # We delay η extraction to a `lift` so it tracks updates to α.
    η_obs = Makie.lift(α_obs, p.waterline_z, p.mask) do α, wl_z, mask
        nx, ny, nz = size(α)
        wl = wl_z === Makie.automatic ? nz / 2 : wl_z
        eta_from_alpha(α, wl;
            mask_inside = mask isa Function ? mask : nothing)
    end
    xs_obs = Makie.lift(α_obs) do α
        cell_centres(size(α, 1))
    end
    ys_obs = Makie.lift(α_obs) do α
        cell_centres(size(α, 2))
    end
    Makie.heatmap!(p, xs_obs, ys_obs, η_obs;
        colormap   = p.colormap,
        colorrange = p.colorrange,
        nan_color  = p.nan_color,
        highclip   = p.highclip,
        lowclip    = p.lowclip,
    )
    # Hull box overlay (optional)
    if p.hull_box[] !== nothing
        box_obs = Makie.lift(p.hull_box) do hb
            xc, yc, L, B = hb
            [(xc - L/2, yc - B/2),
             (xc + L/2, yc - B/2),
             (xc + L/2, yc + B/2),
             (xc - L/2, yc + B/2),
             (xc - L/2, yc - B/2)]
        end
        Makie.lines!(p, box_obs; color = :black, linewidth = 1.5)
    end
    return p
end
