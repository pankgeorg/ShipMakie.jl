"""
    streamlines3d(u, seeds; nsteps = 100, dt = 0.5, color_by = :speed)

Plot 3D streamlines integrated forward from a list of `seeds` (each a
`Point3f` or `(x, y, z)` tuple) through the WaterLily face-staggered
velocity array `u::Array{T, 4}`. Uses RK4 integration with trilinear
interpolation onto the cell-centre velocity.

The streamline is truncated when it leaves the domain interior.

# Attributes

- `seeds` — `Vector{Point3f}` of starting positions in cell-centre
  coords. Use `random_seeds(grid_size, n)` (in `utils.jl`) for a
  scatter, or hand-craft for specific features.
- `nsteps` — max number of integration steps per seed (default 100).
- `dt` — integration step size in cell-time units (default 0.5).
- `color_by` — `:speed` (default), `:index` (per-line constant), or a
  fixed `Colorant`.
- `colormap` — default `:viridis`.
- `linewidth` — default 1.5.

# Example

```julia
seeds = [Point3f(prop_xc + 0.1, prop_yc + r*cos(θ), prop_zc + r*sin(θ))
         for r in 0.3:0.3:1.5, θ in 0:π/4:2π][:]
streamlines3d!(ax, sim.flow.u, seeds; nsteps = 200, dt = 0.4)
```
"""
Makie.@recipe(Streamlines3d, u, seeds) do scene
    Makie.Attributes(
        nsteps    = 100,
        dt        = 0.5,
        color_by  = :speed,
        colormap  = :viridis,
        linewidth = 1.5,
    )
end

# Trilinear interpolation of a face-staggered velocity component at a
# cell-centre world coordinate `p = (x, y, z)`. Returns Float32.
@inline function _trilinear_u(u::AbstractArray{T, 4}, p::NTuple{3, Float64},
                              comp::Int) where T
    nx, ny, nz, _ = size(u)
    gx = clamp(Float32(p[1] + 1.5), 1f0, Float32(nx - 1))
    gy = clamp(Float32(p[2] + 1.5), 1f0, Float32(ny - 1))
    gz = clamp(Float32(p[3] + 1.5), 1f0, Float32(nz - 1))
    i = clamp(floor(Int, gx), 1, nx - 1)
    j = clamp(floor(Int, gy), 1, ny - 1)
    k = clamp(floor(Int, gz), 1, nz - 1)
    fx, fy, fz = gx - i, gy - j, gz - k
    @inbounds c00 = u[i,   j,   k,   comp] * (1 - fx) + u[i+1, j,   k,   comp] * fx
    @inbounds c01 = u[i,   j,   k+1, comp] * (1 - fx) + u[i+1, j,   k+1, comp] * fx
    @inbounds c10 = u[i,   j+1, k,   comp] * (1 - fx) + u[i+1, j+1, k,   comp] * fx
    @inbounds c11 = u[i,   j+1, k+1, comp] * (1 - fx) + u[i+1, j+1, k+1, comp] * fx
    c0 = c00 * (1 - fy) + c10 * fy
    c1 = c01 * (1 - fy) + c11 * fy
    return Float32(c0 * (1 - fz) + c1 * fz)
end

@inline function _velocity_at(u::AbstractArray, p::NTuple{3, Float64})
    return (_trilinear_u(u, p, 1),
            _trilinear_u(u, p, 2),
            _trilinear_u(u, p, 3))
end

function _integrate_streamline(u::AbstractArray, seed::NTuple{3, Float64},
                               nsteps::Int, dt::Float64)
    nx, ny, nz, _ = size(u)
    pts = Tuple{Float64, Float64, Float64}[seed]
    speeds = Float64[]
    p = seed
    for _ in 1:nsteps
        # RK4 step
        v1 = _velocity_at(u, p)
        if !all(isfinite, v1); break; end
        p2 = p .+ (0.5 * dt) .* v1
        v2 = _velocity_at(u, p2)
        p3 = p .+ (0.5 * dt) .* v2
        v3 = _velocity_at(u, p3)
        p4 = p .+ dt .* v3
        v4 = _velocity_at(u, p4)
        dp = (
            dt / 6 * (v1[1] + 2 * v2[1] + 2 * v3[1] + v4[1]),
            dt / 6 * (v1[2] + 2 * v2[2] + 2 * v3[2] + v4[2]),
            dt / 6 * (v1[3] + 2 * v2[3] + 2 * v3[3] + v4[3]),
        )
        p = (p[1] + dp[1], p[2] + dp[2], p[3] + dp[3])
        # Bail if out of bounds
        if p[1] < 0 || p[2] < 0 || p[3] < 0 ||
           p[1] > nx - 2 || p[2] > ny - 2 || p[3] > nz - 2
            break
        end
        push!(pts, p)
        push!(speeds, sqrt(v1[1]^2 + v1[2]^2 + v1[3]^2))
    end
    return pts, speeds
end

function Makie.plot!(p::Streamlines3d)
    u_obs = p[1]
    seeds_obs = p[2]
    lines_obs = Makie.lift(u_obs, seeds_obs, p.nsteps, p.dt, p.color_by) do u, seeds, ns, δt, cb
        all_pts = Tuple{Float64, Float64, Float64}[]
        all_cols = Float64[]
        for (idx, s) in enumerate(seeds)
            seed = s isa NTuple{3} ? (Float64(s[1]), Float64(s[2]), Float64(s[3])) :
                                      (Float64(s[1]), Float64(s[2]), Float64(s[3]))
            pts, speeds = _integrate_streamline(u, seed, ns, δt)
            for k in 1:length(pts)
                push!(all_pts, pts[k])
                if cb === :speed
                    push!(all_cols, k ≤ length(speeds) ? Float64(speeds[k]) :
                          (length(speeds) > 0 ? Float64(speeds[end]) : 0.0))
                elseif cb === :index
                    push!(all_cols, Float64(idx))
                else
                    push!(all_cols, 0.0)
                end
            end
            # NaN-break between streamlines
            push!(all_pts, (NaN, NaN, NaN))
            push!(all_cols, NaN)
        end
        return (all_pts, all_cols)
    end
    Makie.lines!(p,
        Makie.lift(lc -> lc[1], lines_obs);
        color     = Makie.lift(lc -> lc[2], lines_obs),
        colormap  = p.colormap,
        linewidth = p.linewidth,
    )
    return p
end

"""
    random_seeds(grid_size::NTuple{3, Int}, n::Int; rng = nothing)

Convenience: `n` random seeds uniformly distributed in the interior
of a domain of size `grid_size`. Returns `Vector{NTuple{3, Float64}}`.
"""
function random_seeds(grid_size::NTuple{3, Int}, n::Int; rng = nothing)
    Nx, Ny, Nz = grid_size
    seeds = Vector{NTuple{3, Float64}}(undef, n)
    for i in 1:n
        seeds[i] = (
            (Nx - 4) * rand() + 1,
            (Ny - 4) * rand() + 1,
            (Nz - 4) * rand() + 1,
        )
    end
    return seeds
end
