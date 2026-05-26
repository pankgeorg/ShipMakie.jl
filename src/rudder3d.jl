"""
    rudder3d(chord, span; center = (0, 0, 0), δ = 0,
                       rudder_axis = (0, 0, 1))

Plot a 3D wireframe of a flat-plate `Rudder`. The plate is `chord` ×
`span`, centred at `center`, with its span aligned with `rudder_axis`
(default +z, i.e. vertical). The plate is rotated by `δ` radians
about `rudder_axis` to show the deflection angle.

# Attributes

- `δ` — rudder angle in radians, default 0.
- `center` — world-frame centre point.
- `rudder_axis` — span direction (unit vector, default +z).
- `flow_axis` — chord direction at δ=0 (default +x).
- `color` — outline colour, default `:lime`.
- `linewidth` — default 2.

# Example

```julia
rudder3d!(ax, rudder.chord, rudder.span;
    center = (rud_xc, rud_yc, rud_zc), δ = deg2rad(10))
```
"""
Makie.@recipe(Rudder3d, chord, span) do scene
    Makie.Attributes(
        δ           = 0.0,
        center      = (0.0, 0.0, 0.0),
        rudder_axis = (0.0, 0.0, 1.0),
        flow_axis   = (1.0, 0.0, 0.0),
        color       = :lime,
        linewidth   = 2.5,
    )
end

function Makie.plot!(p::Rudder3d)
    chord_obs = p[1]
    span_obs  = p[2]
    quad_obs = Makie.lift(chord_obs, span_obs, p.center, p.rudder_axis,
                          p.flow_axis, p.δ) do chord, span, c, sax, fax, δ
        sax = Float64.(sax); fax = Float64.(fax); c = Float64.(c)
        # Build orthonormal frame: span = sax, "chord at δ=0" = fax
        # Rotation about sax by δ rotates fax in the plane perpendicular
        # to sax.
        # n_sax × fax gives the perpendicular component.
        cδ, sδ = cos(δ), sin(δ)
        chord_dir = (
            cδ * fax[1] + sδ * (sax[2] * fax[3] - sax[3] * fax[2]),
            cδ * fax[2] + sδ * (sax[3] * fax[1] - sax[1] * fax[3]),
            cδ * fax[3] + sδ * (sax[1] * fax[2] - sax[2] * fax[1]),
        )
        half_c = chord / 2
        half_s = span / 2
        # Four corners of the plate
        corners = [
            (c[1] - half_c * chord_dir[1] - half_s * sax[1],
             c[2] - half_c * chord_dir[2] - half_s * sax[2],
             c[3] - half_c * chord_dir[3] - half_s * sax[3]),
            (c[1] + half_c * chord_dir[1] - half_s * sax[1],
             c[2] + half_c * chord_dir[2] - half_s * sax[2],
             c[3] + half_c * chord_dir[3] - half_s * sax[3]),
            (c[1] + half_c * chord_dir[1] + half_s * sax[1],
             c[2] + half_c * chord_dir[2] + half_s * sax[2],
             c[3] + half_c * chord_dir[3] + half_s * sax[3]),
            (c[1] - half_c * chord_dir[1] + half_s * sax[1],
             c[2] - half_c * chord_dir[2] + half_s * sax[2],
             c[3] - half_c * chord_dir[3] + half_s * sax[3]),
        ]
        # Close the loop
        push!(corners, corners[1])
        return corners
    end
    Makie.lines!(p, quad_obs;
        color = p.color, linewidth = p.linewidth)
    return p
end
