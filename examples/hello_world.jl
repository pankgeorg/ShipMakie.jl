#!/usr/bin/env julia
#
# hello_world.jl — the smallest possible ShipMakie example.
# Uses synthetic data (no WaterLily) so it runs anywhere.
#
# Builds a single Figure with one panel per recipe, on a 16×16×16
# grid with hand-crafted fields. Useful as a recipe-by-recipe sanity
# check and as a quickstart for new users.

import Pkg
Pkg.activate(; temp=true)
Pkg.develop(Pkg.PackageSpec(path = joinpath(@__DIR__, "..")); io = devnull)
Pkg.add("CairoMakie"; io = devnull)

using ShipMakie
using CairoMakie

# ---------------------------------------------------------------------------
# Synthetic data
# ---------------------------------------------------------------------------
const NX, NY, NZ = 16, 16, 12

# α: half-water-half-air with a sinusoidal surface
α = zeros(Float32, NX, NY, NZ)
for k in 1:NZ, j in 1:NY, i in 1:NX
    z = k - 1.5f0
    η_target = NZ/2 + 0.5f0 * sinpi(2f0 * (i + j) / NX)
    α[i, j, k] = z ≤ η_target ? 1f0 : 0f0
end

# u: uniform inflow in +x plus a swirl about the y-axis
u = zeros(Float32, NX, NY, NZ, 3)
for k in 1:NZ, j in 1:NY, i in 1:NX
    x_c = i - NX/2; z_c = k - NZ/2
    r = sqrt(x_c^2 + z_c^2)
    θ = atan(z_c, x_c)
    v = 0.4f0 * exp(-r/5)
    u[i, j, k, 1] = 1f0 - v * sin(θ)
    u[i, j, k, 3] = +v * cos(θ)
end

# p: pressure field with high-pressure stagnation at the bow
p = zeros(Float32, NX, NY, NZ)
for k in 1:NZ, j in 1:NY, i in 1:NX
    x_c = i - 4f0; y_c = j - NY/2; z_c = k - NZ/2
    r = sqrt(x_c^2 + y_c^2 + z_c^2)
    p[i, j, k] = 1.5f0 * exp(-r/3)
end

# Sphere SDF for the hull
sphere_sdf(p) = sqrt((p[1] - 8)^2 + (p[2] - NY/2)^2 + (p[3] - NZ/2)^2) - 3.0

# ---------------------------------------------------------------------------
# Composite figure
# ---------------------------------------------------------------------------
fig = Figure(size = (1400, 900), backgroundcolor = :white)

# Row 1: 2D recipes
let ax = Axis(fig[1, 1]; title = "etaheatmap", aspect = DataAspect())
    etaheatmap!(ax, α; waterline_z = NZ/2)
end
let ax = Axis(fig[1, 2]; title = "velocityslice (u_x, y-mid)", aspect = DataAspect())
    velocityslice!(ax, u; slice_axis = :y, component = 1)
end
let ax = Axis(fig[1, 3]; title = "streamslice (y-mid)", aspect = DataAspect())
    streamslice!(ax, u; slice_axis = :y, density = 0.6)
end
let ax = Axis(fig[1, 4]; title = "vectorslice (y-mid)", aspect = DataAspect())
    vectorslice!(ax, u; slice_axis = :y, stride = 2)
end

# Row 2: 1D + 3D recipes
let ax = Axis(fig[2, 1]; title = "probeline (α along z at i=8, j=8)")
    probeline!(ax, α; along = :z, fixed_indices = (8, 8))
end
let ax = Axis(fig[2, 2]; title = "hullsilhouette (y-mid)", aspect = DataAspect())
    hullsilhouette!(ax, sphere_sdf; grid_size = (NX, NY, NZ), slice_axis = :y)
end
let ax = Axis3(fig[2, 3]; title = "vorticityvolume (Q, MIP)")
    vorticityvolume!(ax, u; field = :q_criterion, algorithm = :mip)
end
let ax = Axis3(fig[2, 4]; title = "hullmesh3d (sphere)")
    hullmesh3d!(ax, sphere_sdf; grid_size = (NX, NY, NZ))
end

Label(fig[0, 1:4]; text = "ShipMakie hello-world  •  synthetic data, every recipe",
    fontsize = 18, halign = :center)

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "hello_world"))
mkpath(OUTDIR)
out = joinpath(OUTDIR, "hello.png")
save(out, fig)
println("Saved $out")
