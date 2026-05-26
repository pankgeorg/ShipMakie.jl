"""
    hullmesh3d(sdf; grid_size = (32, 32, 32), level = 0.0)

Plot the hull's 3D surface as an isosurface of the SDF at `level = 0`.
Uses Makie's 3D `contour!` (Marching Cubes under the hood).

# Inputs

- `sdf` — a callable `sdf(p::SVector{3}) → Real` (or `sdf(p, t) → Real`,
  we sample at `t = 0`).
- `grid_size` — `(Nx, Ny, Nz)` cells to sample the SDF on.

# Attributes

- `level` — SDF level to render (default 0; the hull surface).
- `color` — surface colour, default `:slategray3`.
- `transparency` — `true` for see-through, default `false`.
- `colormap` — only used when `color` is a vector; default `:slategray3`.

# Example

```julia
fig = Figure(); ax = Axis3(fig[1, 1])
hullmesh3d!(ax, p -> wigley_sdf(p, L, B, T); grid_size = (NX, NY, NZ))
```
"""
Makie.@recipe(HullMesh3d, sdf) do scene
    Makie.Attributes(
        grid_size    = (32, 32, 32),
        level        = 0.0,
        color        = :slategray3,
        transparency = false,
        colormap     = :greys,
    )
end

function MakieEvalSDF(sdf, x, y, z)
    p = SVector{3, Float64}(x, y, z)
    try
        return sdf(p)
    catch
        return sdf(p, 0.0)
    end
end

function Makie.plot!(p::HullMesh3d)
    sdf = p[1][]
    field_obs = Makie.lift(p.grid_size) do gs
        Nx, Ny, Nz = gs
        f = zeros(Float64, Nx, Ny, Nz)
        @inbounds for k in 1:Nz, j in 1:Ny, i in 1:Nx
            f[i, j, k] = MakieEvalSDF(sdf,
                Float64(i - 1.5), Float64(j - 1.5), Float64(k - 1.5))
        end
        return f
    end
    xs_obs = Makie.lift(p.grid_size) do gs
        c = cell_centres(gs[1]); (Float64(c[1]), Float64(c[end]))
    end
    ys_obs = Makie.lift(p.grid_size) do gs
        c = cell_centres(gs[2]); (Float64(c[1]), Float64(c[end]))
    end
    zs_obs = Makie.lift(p.grid_size) do gs
        c = cell_centres(gs[3]); (Float64(c[1]), Float64(c[end]))
    end
    Makie.contour!(p, xs_obs, ys_obs, zs_obs, field_obs;
        levels       = Makie.lift(p.level) do l; [l]; end,
        color        = p.color,
        colormap     = p.colormap,
        transparency = p.transparency,
    )
    return p
end
