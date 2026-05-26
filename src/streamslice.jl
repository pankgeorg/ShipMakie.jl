"""
    streamslice(u; slice_axis = :y, index = automatic, ...)

Plot 2D streamlines on a coordinate-aligned slice of a WaterLily
face-staggered velocity array `u::Array{T, 4}` (size `(NX, NY, NZ, 3)`).

The streamline plot uses Makie's `streamplot` recipe under the hood;
this wrapper handles the staggered-cell → cell-centre conversion and
takes the slice along the chosen axis.

# Attributes

- `slice_axis` — `:x`, `:y`, or `:z` (default `:y`, centreline side view).
- `index` — slice index along `slice_axis` (default midplane).
- `density` — streamline density (Makie default 1.0).
- `colormap` — streamline colour, default `:viridis`.
- `linewidth` — default `1`.
- `arrow_size` — default `8`.

# Example

```julia
# centreline streamlines (side view), velocity in the x-z plane
streamslice(sim.flow.u; slice_axis = :y, index = NY ÷ 2)
```
"""
Makie.@recipe(StreamSlice, u) do scene
    Makie.Attributes(
        slice_axis = :y,
        index      = Makie.automatic,
        density    = 1.0,
        colormap   = :viridis,
        linewidth  = 1.0,
        arrow_size = 8.0,
    )
end

function Makie.plot!(p::StreamSlice)
    u_obs = p[1]
    # Extract the two velocity components in the slice plane.
    # For slice_axis = :y at fixed j, the plane is (x, z) with
    # components u_x and u_z. Etc.
    plane_obs = Makie.lift(u_obs, p.slice_axis, p.index) do u, ax, idx
        ax_i = _axis_idx(ax)
        sz = size(u)
        k = idx === Makie.automatic ? sz[ax_i] ÷ 2 : idx
        # Pick the two non-slice component indices for plotting
        if ax_i == 1
            return (collect(@view u[k, :, :, 2]),  # u_y
                    collect(@view u[k, :, :, 3]))  # u_z
        elseif ax_i == 2
            return (collect(@view u[:, k, :, 1]),  # u_x
                    collect(@view u[:, k, :, 3]))  # u_z
        else
            return (collect(@view u[:, :, k, 1]),  # u_x
                    collect(@view u[:, :, k, 2]))  # u_y
        end
    end
    # Build a closure (x, y) -> Point2f for streamplot.
    # `streamplot` requires a function form; we sample the (u_a, u_b)
    # arrays at the nearest cell.
    bounds_obs = Makie.lift(u_obs, p.slice_axis) do u, ax
        ax_i = _axis_idx(ax)
        # Two non-slice axes determine x and y of the plot.
        if ax_i == 1
            return (size(u, 2), size(u, 3))
        elseif ax_i == 2
            return (size(u, 1), size(u, 3))
        else
            return (size(u, 1), size(u, 2))
        end
    end
    fn_obs = Makie.lift(plane_obs, bounds_obs) do (ua, ub), (Na, Nb)
        function f(pt)
            # pt[1], pt[2] in cell-centre world coords; convert to cell index
            ia = clamp(round(Int, pt[1] + 1.5), 1, Na)
            ib = clamp(round(Int, pt[2] + 1.5), 1, Nb)
            return Point2f(ua[ia, ib], ub[ia, ib])
        end
        return f
    end
    xrange_obs = Makie.lift(bounds_obs) do (Na, _Nb)
        (cell_centres(Na)[1], cell_centres(Na)[end])
    end
    yrange_obs = Makie.lift(bounds_obs) do (_Na, Nb)
        (cell_centres(Nb)[1], cell_centres(Nb)[end])
    end
    Makie.streamplot!(p, fn_obs, xrange_obs, yrange_obs;
        density    = p.density,
        colormap   = p.colormap,
        linewidth  = p.linewidth,
        arrow_size = p.arrow_size,
    )
    return p
end
