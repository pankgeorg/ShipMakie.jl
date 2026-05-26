"""
    vectorslice(u; slice_axis = :y, index = automatic, stride = 4)

Plot a 2D vector arrow field on a coord-aligned slice of a WaterLily
face-staggered velocity array. Useful alongside `velocityslice` for
showing direction in addition to magnitude.

# Attributes

- `slice_axis` — `:x`, `:y`, or `:z` (default `:y`).
- `index` — slice index (default midplane).
- `stride` — show one arrow every `stride` cells (default 4, to avoid
  clutter).
- `lengthscale` — arrow length multiplier, default 1.
- `colormap` — arrow colour map (by magnitude), default `:viridis`.
- `arrowsize` — head size, default 8.

# Example

```julia
vectorslice!(ax, sim.flow.u; slice_axis = :y, stride = 6)
```
"""
Makie.@recipe(VectorSlice, u) do scene
    Makie.Attributes(
        slice_axis  = :y,
        index       = Makie.automatic,
        stride      = 4,
        lengthscale = 1.0,
        colormap    = :viridis,
        arrowsize   = 8.0,
    )
end

function Makie.plot!(p::VectorSlice)
    u_obs = p[1]
    arrows_obs = Makie.lift(u_obs, p.slice_axis, p.index, p.stride) do u, ax, idx, st
        ax_i = _axis_idx(ax)
        sz = size(u)
        k = idx === Makie.automatic ? sz[ax_i] ÷ 2 : idx
        # Pick the two non-slice component arrays
        if ax_i == 1
            (Na, Nb) = (sz[2], sz[3])
            ua = @view u[k, :, :, 2]
            ub = @view u[k, :, :, 3]
        elseif ax_i == 2
            (Na, Nb) = (sz[1], sz[3])
            ua = @view u[:, k, :, 1]
            ub = @view u[:, k, :, 3]
        else
            (Na, Nb) = (sz[1], sz[2])
            ua = @view u[:, :, k, 1]
            ub = @view u[:, :, k, 2]
        end
        # Sample every `st` cells, skipping the ghost rim
        xs   = Float32[];   ys   = Float32[]
        us   = Float32[];   vs   = Float32[]
        mags = Float32[]
        for j in 2:st:Nb-1, i in 2:st:Na-1
            push!(xs, Float32(i - 1.5))
            push!(ys, Float32(j - 1.5))
            push!(us, ua[i, j])
            push!(vs, ub[i, j])
            push!(mags, sqrt(ua[i, j]^2 + ub[i, j]^2))
        end
        return (xs, ys, us, vs, mags)
    end
    Makie.arrows!(p,
        Makie.lift(a -> a[1], arrows_obs),  # xs
        Makie.lift(a -> a[2], arrows_obs),  # ys
        Makie.lift(a -> a[3], arrows_obs),  # us
        Makie.lift(a -> a[4], arrows_obs);  # vs
        lengthscale = p.lengthscale,
        arrowsize   = p.arrowsize,
        color       = Makie.lift(a -> a[5], arrows_obs),
        colormap    = p.colormap,
    )
    return p
end
