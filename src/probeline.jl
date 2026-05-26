"""
    probeline(field::AbstractArray{T, D}, along::Symbol, fixed_indices...)

Plot a 1D profile of a scalar `field` along a coord-aligned line.
For 3D fields, supply two `fixed_indices` (the other two axes' cell
indices); the plot's x-along is the cell-centre coord along `along`.

# Examples

```julia
# centreline pressure profile along x
probeline(sim.flow.p, :x, NY÷2, NZ÷2)

# vertical α profile in a column
probeline(vof.α, :z, hull_xc_idx, hull_yc_idx)
```

# Attributes

- `color` — line colour, default `:steelblue`.
- `linewidth` — default 2.
"""
Makie.@recipe(ProbeLine, field) do scene
    Makie.Attributes(
        along           = :x,
        fixed_indices  = nothing,
        color          = :steelblue,
        linewidth      = 2.0,
    )
end

function Makie.plot!(p::ProbeLine)
    f_obs = p[1]
    line_obs = Makie.lift(f_obs, p.along, p.fixed_indices) do f, ax, fi
        ax_i = _axis_idx(ax)
        D = ndims(f)
        if D == 3
            fi === nothing && throw(ArgumentError("3D field requires fixed_indices = (i, j)"))
            i1, i2 = fi
            return ax_i == 1 ? collect(@view f[:, i1, i2]) :
                   ax_i == 2 ? collect(@view f[i1, :, i2]) :
                               collect(@view f[i1, i2, :])
        elseif D == 2
            fi === nothing && throw(ArgumentError("2D field requires fixed_indices = (i,)"))
            i1 = fi isa Tuple ? fi[1] : fi
            return ax_i == 1 ? collect(@view f[:, i1]) :
                               collect(@view f[i1, :])
        elseif D == 1
            return collect(f)
        else
            throw(ArgumentError("probeline supports 1D/2D/3D fields, got D=$D"))
        end
    end
    xs_obs = Makie.lift(f_obs, p.along) do f, ax
        cell_centres(size(f, _axis_idx(ax)))
    end
    Makie.lines!(p, xs_obs, line_obs;
        color     = p.color,
        linewidth = p.linewidth,
    )
    return p
end
