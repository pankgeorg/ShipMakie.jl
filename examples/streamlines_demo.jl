#!/usr/bin/env julia
#
# streamlines_demo: 3D streamlines around a Wigley hull. Particles are
# seeded in a transverse plane upstream of the bow and integrated
# forward through the flow field. Coloured by speed.

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path = joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path = joinpath(@__DIR__, "..", "..", "VoF.jl")),
    Pkg.PackageSpec(path = joinpath(@__DIR__, "..", "..", "ShipShapes.jl")),
    Pkg.PackageSpec(path = joinpath(@__DIR__, "..", "..", "LiftingSurfaces.jl")),
    Pkg.PackageSpec(path = joinpath(@__DIR__, "..")),
]; io = devnull)
Pkg.add(["VortexLattice", "GLMakie"]; io = devnull)

using WaterLily, VoF, ShipShapes, LiftingSurfaces, ShipMakie
using ShipShapes: StaticArrays
const SVector = StaticArrays.SVector
using Printf
using GLMakie
GLMakie.activate!()

const NX, NY, NZ = 96, 48, 32
const NSTEPS = 30
const L_c = 28f0; const B_c = 6f0; const T_c = 4f0
const ρ_w = 10f0; const ρ_a = 1f0
const U∞ = 1f0
const Fr = 0.30f0
const G_c = U∞^2 / (Fr^2 * L_c)
const Re = 5000f0
const ν_w_c = U∞ * L_c / Re
const μ_w_c = ρ_w * ν_w_c
const μ_a_c = ρ_a * 18 * ν_w_c
const H_w_c = Float32(NZ / 2)
const hull_xc = Float32(NX / 5); const hull_yc = Float32(NY / 2); const hull_zc = H_w_c
const R_prop  = 2.5f0
const prop_xc = Float32(hull_xc + L_c / 2 + T_c)
const prop_yc = hull_yc
const prop_zc = Float32(H_w_c - T_c / 2)
const J_op = 0.32
const Ω    = Float64(π) * U∞ / (J_op * R_prop)

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1f0 : 0f0
rotor = BladedRotor(; N_blades = 3, R = R_prop, R_hub = 0.5,
    chord = (0.8, 0.4), twist = (deg2rad(35.0), deg2rad(15.0)),
    ns = 12, nc = 4)
r_rot = rotor_forces(rotor, U∞, Ω)
Sref_r = π * R_prop^2
thrust = abs(r_rot.CT * 0.5 * U∞^2 * Sref_r)
torque = r_rot.CQ * 0.5 * U∞^2 * Sref_r * R_prop

vof = VoFFlow((NX, NY, NZ);
    α₀ = α₀, ρ_w = ρ_w, ρ_a = ρ_a, μ_w = μ_w_c, μ_a = μ_a_c, T = Float32)
function vof_pois_ctor(flow)
    L = similar(flow.μ₀)
    @inbounds for I in CartesianIndices(L)
        L[I] = flow.μ₀[I] * vof.L[I]
    end
    WaterLily.MultiLevelPoisson(flow.p, L, flow.σ; perdir = (2,))
end
sim = WaterLily.Simulation((NX, NY, NZ),
    (U∞, 0f0, 0f0), L_c;
    T = Float32, ν = vof.ν,
    g = (i, x, t) -> i == 3 ? -G_c : 0f0,
    Δt = 0.25f0, body = hull, ϵ = 1, perdir = (2,), exitBC = true,
    pois_ctor = vof_pois_ctor, U = U∞,
)
function vlm_udf(flow, t; kwargs...)
    smear_blades!(flow.f, thrust, torque,
        SVector(prop_xc, prop_yc, prop_zc),
        SVector(1.0, 0.0, 0.0),
        Float64(R_prop), 0.2 * Float64(R_prop);
        N_blades = 3, N_sections = 4, ε = 1.5)
    return nothing
end

@info "Running $NSTEPS steps to develop the flow field…"
for _ in 1:NSTEPS
    WaterLily.mom_step!(sim.flow, sim.pois; udf = vlm_udf,
        pois_tol = 1f-6, pois_itmx = 50)
    step_vof_mules!(vof, sim; dt = sim.flow.Δt[end-1], perdir = (2,))
end

# ---------------------------------------------------------------------------
# Seed pattern: grid of particles in a y-z plane upstream of the bow
# ---------------------------------------------------------------------------
x_seed = hull_xc - L_c / 2 - 3.0    # 3 cells upstream of the bow
ny_seeds, nz_seeds = 6, 5
seeds = NTuple{3, Float64}[]
for i in 1:ny_seeds, j in 1:nz_seeds
    y = hull_yc - B_c * 1.3 + 2 * (B_c * 1.3) * (i - 0.5) / ny_seeds
    z = hull_zc - T_c * 1.0 + 2 * (T_c * 1.0) * (j - 0.5) / nz_seeds
    push!(seeds, (Float64(x_seed), Float64(y), Float64(z)))
end
@info "Seeded $(length(seeds)) streamlines upstream of the bow"

hull_sdf_world(p) = ShipShapes.wigley_sdf(
    SVector(p[1] - hull_xc, p[2] - hull_yc, p[3] - hull_zc),
    Float64(L_c), Float64(B_c), Float64(T_c))

fig = Figure(size = (1400, 900), backgroundcolor = :gray12)
ax = Axis3(fig[1, 1];
    title = "Streamlines around Wigley + rotor",
    titlecolor = :white,
    aspect  = :data,
    azimuth = 0.35π, elevation = 0.20π,
    xlabelcolor = :white, ylabelcolor = :white, zlabelcolor = :white,
    xticklabelcolor = :white, yticklabelcolor = :white, zticklabelcolor = :white,
    backgroundcolor = :gray12,
)
# Hull (opaque)
hullmesh3d!(ax, hull_sdf_world;
    grid_size = (NX, NY, NZ), color = :navy)
# Rotor wireframe
bladedrotor3d!(ax, 3, R_prop, 0.5 * R_prop;
    center = (prop_xc, prop_yc, prop_zc),
    rotor_axis = (1.0, 0.0, 0.0),
    color = :orange, linewidth = 2.5)
# Streamlines
streamlines3d!(ax, sim.flow.u, seeds;
    nsteps = 250, dt = 0.4,
    colormap = :plasma, linewidth = 2.0)

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "streamlines_demo"))
mkpath(OUTDIR)
out = joinpath(OUTDIR, "streamlines.png")
save(out, fig)
@printf "Saved %s\n" out
