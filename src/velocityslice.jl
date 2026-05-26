"""
    velocityslice(u, slice_axis::Symbol = :y, index::Int = size(u, slice_axis_idx(slice_axis)) ÷ 2;
                  component::Int = 1)

Plot a 2D heatmap of one velocity component on a coordinate-aligned
slice of a WaterLily face-staggered velocity array `u::Array{T, 4}`
(size `(NX, NY, NZ, 3)`).

# Attributes

- `slice_axis` — `:x`, `:y`, or `:z`; which slice_axis to slice perpendicular to.
- `index` — integer cell index along `slice_axis`. Defaults to the midplane.
- `component` — which u-component to plot (1=u_x, 2=u_y, 3=u_z).
- `colormap` — default `:viridis`.
- `colorrange` — default `automatic`.

The plot uses cell-centre world coordinates on the two remaining axes.

# Example

```julia
# centreline u_x slice (side view)
velocityslice(sim.flow.u, :y, size(sim.flow.u, 2) ÷ 2; component = 1)
```
"""
Makie.@recipe(VelocitySlice, u) do scene
    Makie.Attributes(
        slice_axis       = :y,
        index      = Makie.automatic,
        component  = 1,
        colormap   = :viridis,
        colorrange = Makie.automatic,
    )
end

# Helper: slice_axis symbol → integer (1=x, 2=y, 3=z)
@inline _axis_idx(s::Symbol) = s === :x ? 1 : s === :y ? 2 : s === :z ? 3 :
    throw(ArgumentError("slice_axis must be :x, :y, or :z; got $(repr(s))"))

function Makie.plot!(p::VelocitySlice)
    u_obs = p[1]
    plane_obs = Makie.lift(u_obs, p.slice_axis, p.index, p.component) do u, ax, idx, c
        ax_i  = _axis_idx(ax)
        sz    = size(u)
        k     = idx === Makie.automatic ? sz[ax_i] ÷ 2 : idx
        # Build a 2D slice. Other-slice_axis ranges, then component.
        # u has shape (NX, NY, NZ, 3); slice along ax_i at index k.
        if ax_i == 1
            return collect(@view u[k, :, :, c])
        elseif ax_i == 2
            return collect(@view u[:, k, :, c])
        else
            return collect(@view u[:, :, k, c])
        end
    end
    xs_obs = Makie.lift(u_obs, p.slice_axis) do u, ax
        ax_i = _axis_idx(ax)
        # First "other" slice_axis after ax_i (the slice_axis Makie shows on x of plot)
        other = ax_i == 1 ? 2 : 1
        cell_centres(size(u, other))
    end
    ys_obs = Makie.lift(u_obs, p.slice_axis) do u, ax
        ax_i = _axis_idx(ax)
        # Second "other" slice_axis (the y of plot)
        other = ax_i == 3 ? 2 : 3
        cell_centres(size(u, other))
    end
    Makie.heatmap!(p, xs_obs, ys_obs, plane_obs;
        colormap   = p.colormap,
        colorrange = p.colorrange,
    )
    return p
end
