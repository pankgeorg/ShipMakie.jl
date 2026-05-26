"""
    hullsilhouette(sdf, grid_size::NTuple{2, Int}, slice_axis::Symbol = :z,
                   slice_index::Int = grid_size[3] ÷ 2;
                   level = 0.0)

Plot the body silhouette derived from a signed-distance function
`sdf(p) → Real` evaluated on a 2D plane through the simulation domain.
The contour at `sdf = level` is the body surface.

# Inputs

- `sdf` — a callable `sdf(p::SVector{3}) → Real`, or `sdf(p, t) → Real`
  if your AutoBody convention requires the time arg. We probe at
  `t = 0` for the time-dependent variant.
- `grid_size` — `(Nx, Ny, Nz)` to set the world-coord range.
- `slice_axis` — which axis is fixed (`:x`, `:y`, or `:z`).
- `slice_index` — the cell index along the fixed axis.
- `level` — the SDF contour level (default 0; body surface).

# Attributes

- `color` — silhouette stroke colour, default `:black`.
- `linewidth` — default 1.5.

# Example

```julia
hullsilhouette(p -> wigley_sdf(p, L, B, T),
               (NX, NY, NZ), :y, NY ÷ 2)
```
"""
Makie.@recipe(HullSilhouette, sdf) do scene
    Makie.Attributes(
        grid_size   = (32, 32, 32),
        slice_axis  = :z,
        slice_index = Makie.automatic,
        level       = 0.0,
        color       = :black,
        linewidth   = 1.5,
    )
end

# Evaluate the user's SDF at world point p. Handles both
# `sdf(p)` and `sdf(p, t)` signatures.
@inline function _eval_sdf(sdf, p)
    try
        return sdf(p)
    catch
        return sdf(p, 0.0)
    end
end

function Makie.plot!(p::HullSilhouette)
    sdf = p[1][]
    field_obs = Makie.lift(p.grid_size, p.slice_axis, p.slice_index) do gs, ax, idx
        ax_i = _axis_idx(ax)
        Nx, Ny, Nz = gs
        Ni, Nj = ax_i == 1 ? (Ny, Nz) :
                 ax_i == 2 ? (Nx, Nz) :
                             (Nx, Ny)
        k = idx === Makie.automatic ? gs[ax_i] ÷ 2 : idx
        f = zeros(Float64, Ni, Nj)
        @inbounds for j in 1:Nj, i in 1:Ni
            xx = if ax_i == 1
                SVector{3, Float64}(k - 1.5, i - 1.5, j - 1.5)
            elseif ax_i == 2
                SVector{3, Float64}(i - 1.5, k - 1.5, j - 1.5)
            else
                SVector{3, Float64}(i - 1.5, j - 1.5, k - 1.5)
            end
            f[i, j] = _eval_sdf(sdf, xx)
        end
        f
    end
    xs_obs = Makie.lift(p.grid_size, p.slice_axis) do gs, ax
        ax_i = _axis_idx(ax)
        cell_centres(gs[ax_i == 1 ? 2 : 1])
    end
    ys_obs = Makie.lift(p.grid_size, p.slice_axis) do gs, ax
        ax_i = _axis_idx(ax)
        cell_centres(gs[ax_i == 3 ? 2 : 3])
    end
    Makie.contour!(p, xs_obs, ys_obs, field_obs;
        levels    = Makie.lift(p.level) do l; [l]; end,
        color     = p.color,
        linewidth = p.linewidth,
    )
    return p
end
