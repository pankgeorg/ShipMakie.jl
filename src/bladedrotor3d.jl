"""
    bladedrotor3d(N_blades, R, R_hub; center = (0,0,0), axis = (1,0,0))

Plot a 3D wireframe of a `BladedRotor`'s panel structure: `N_blades`
radial spokes from `R_hub` to `R`, each rotated 2π/N_blades around
`axis`, anchored at `center`.

The recipe takes loose inputs (not a `BladedRotor` struct directly)
so the package doesn't need a LiftingSurfaces dependency. Pass the
fields you have:

```julia
bladedrotor3d(rotor.N_blades, rotor.R, rotor.R_hub;
              center = (prop_xc, prop_yc, prop_zc),
              axis = (1, 0, 0))
```

# Attributes

- `center` — world-frame rotor centre, default origin.
- `axis` — rotor axis unit vector, default `(1, 0, 0)`.
- `color` — wire colour, default `:black`.
- `linewidth` — default 2.
"""
Makie.@recipe(BladedRotor3d, N_blades, R, R_hub) do scene
    Makie.Attributes(
        center     = (0.0, 0.0, 0.0),
        rotor_axis = (1.0, 0.0, 0.0),
        color      = :black,
        linewidth  = 2.0,
    )
end

# Build an orthonormal basis (e1, e2) ⊥ axis.
@inline function _ortho_basis(a::NTuple{3, T}) where T
    n = sqrt(a[1]^2 + a[2]^2 + a[3]^2)
    ax = (a[1] / n, a[2] / n, a[3] / n)
    # Pick e1 from {(1,0,0), (0,1,0)} that's least aligned with ax
    if abs(ax[1]) < 0.9
        e1 = (one(T) - ax[1] * ax[1], -ax[1] * ax[2], -ax[1] * ax[3])
    else
        e1 = (-ax[2] * ax[1], one(T) - ax[2] * ax[2], -ax[2] * ax[3])
    end
    n1 = sqrt(e1[1]^2 + e1[2]^2 + e1[3]^2)
    e1 = (e1[1] / n1, e1[2] / n1, e1[3] / n1)
    # e2 = ax × e1
    e2 = (
        ax[2] * e1[3] - ax[3] * e1[2],
        ax[3] * e1[1] - ax[1] * e1[3],
        ax[1] * e1[2] - ax[2] * e1[1],
    )
    return ax, e1, e2
end

function Makie.plot!(p::BladedRotor3d)
    Nb_obs = p[1]
    R_obs  = p[2]
    Rh_obs = p[3]
    # Build line segments: for each blade, from hub to tip.
    segs_obs = Makie.lift(Nb_obs, R_obs, Rh_obs, p.center, p.rotor_axis) do Nb, R, Rh, c, a
        ax, e1, e2 = _ortho_basis(Float64.(a))
        cc = Float64.(c)
        pts = Tuple{Float64, Float64, Float64}[]
        for b in 0:Nb-1
            θ = 2π * b / Nb
            cθ, sθ = cos(θ), sin(θ)
            r̂ = (cθ * e1[1] + sθ * e2[1],
                  cθ * e1[2] + sθ * e2[2],
                  cθ * e1[3] + sθ * e2[3])
            hub_pt = (cc[1] + Rh * r̂[1],
                      cc[2] + Rh * r̂[2],
                      cc[3] + Rh * r̂[3])
            tip_pt = (cc[1] + R * r̂[1],
                      cc[2] + R * r̂[2],
                      cc[3] + R * r̂[3])
            push!(pts, hub_pt, tip_pt, (NaN, NaN, NaN))
        end
        return pts
    end
    Makie.lines!(p, segs_obs;
        color     = p.color,
        linewidth = p.linewidth,
    )
    return p
end
