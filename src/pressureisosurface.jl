"""
    pressureisosurface(p; level = automatic)

3D isosurface(s) of a scalar pressure field `p::Array{T, 3}`. Useful
to visualise bow-wave high-pressure regions and stern low-pressure
regions in the same scene.

# Attributes

- `level` — isovalue (or vector of levels) for the iso surfaces.
  Default chooses ±0.7·max(|p|) automatically.
- `color` — surface colour (or vector matching `level`).
- `transparency` — default `true` so other recipes can be layered.

# Example

```julia
pressureisosurface!(ax, sim.flow.p)
```
"""
Makie.@recipe(PressureIsosurface, p) do scene
    Makie.Attributes(
        level        = Makie.automatic,
        color        = [:tomato, :steelblue],
        transparency = true,
    )
end

function Makie.plot!(p::PressureIsosurface)
    field_obs = p[1]
    level_obs = Makie.lift(field_obs, p.level) do f, l
        if l === Makie.automatic
            # ±0.7 × the symmetric 95th percentile
            v = sort(abs.(vec(f)))
            q = v[max(1, Int(round(length(v) * 0.95)))]
            return [-0.7q, +0.7q]
        else
            return l isa AbstractVector ? Float64.(l) : [Float64(l)]
        end
    end
    xs_obs = Makie.lift(field_obs) do f
        c = cell_centres(size(f, 1)); (Float64(c[1]), Float64(c[end]))
    end
    ys_obs = Makie.lift(field_obs) do f
        c = cell_centres(size(f, 2)); (Float64(c[1]), Float64(c[end]))
    end
    zs_obs = Makie.lift(field_obs) do f
        c = cell_centres(size(f, 3)); (Float64(c[1]), Float64(c[end]))
    end
    Makie.contour!(p, xs_obs, ys_obs, zs_obs, field_obs;
        levels       = level_obs,
        color        = p.color,
        transparency = p.transparency,
    )
    return p
end
