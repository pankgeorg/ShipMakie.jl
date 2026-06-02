#!/usr/bin/env julia
#
# Particle tracer animation: ~200 particles seeded around the bow
# are advected by the live flow field. Each frame: advance the sim,
# advance the particles, draw a snapshot. Stunning visualization of
# how the flow wraps around the hull and feeds into the rotor race.

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
const NFRAMES = parse(Int, get(ENV, "NFRAMES", "100"))
const BURNIN  = parse(Int, get(ENV, "BURNIN", "15"))
const NPART   = parse(Int, get(ENV, "NPART", "200"))
const TRAIL   = parse(Int, get(ENV, "TRAIL", "12"))  # number of past positions per particle

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
    T = Float32, ν = VoF.viscosity(vof),
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

@info "Burn-in $BURNIN steps…"
for _ in 1:BURNIN
    WaterLily.mom_step!(sim.flow, sim.pois; udf = vlm_udf,
        pois_tol = 1f-6, pois_itmx = 50)
    step_vof_mules!(vof, sim; dt = sim.flow.Δt[end-1], perdir = (2,))
end

# ---------------------------------------------------------------------------
# Particle initialisation: seeded in a transverse plane upstream of the bow
# ---------------------------------------------------------------------------
x_seed = hull_xc - L_c / 2 - 3.0
particles = NTuple{3, Float32}[]
for _ in 1:NPART
    y = hull_yc - B_c * 1.4 + 2 * (B_c * 1.4) * rand(Float32)
    z = hull_zc - T_c * 1.1 + 2 * (T_c * 1.1) * rand(Float32)
    push!(particles, (Float32(x_seed), y, z))
end

# Trails: a (TRAIL × NPART) ring buffer
trail_buffer = fill((NaN32, NaN32, NaN32), TRAIL, NPART)
trail_head = 1

# --- Helper: advance a particle one step ---------------------------
function advance_particle!(particles::Vector{NTuple{3, Float32}},
                           idx::Int, u, dt::Float32)
    p = particles[idx]
    p_dbl = (Float64(p[1]), Float64(p[2]), Float64(p[3]))
    # RK4 step using the streamlines3d helper
    v1 = ShipMakie._velocity_at(u, p_dbl)
    p2 = p_dbl .+ (0.5 * Float64(dt)) .* v1
    v2 = ShipMakie._velocity_at(u, p2)
    p3 = p_dbl .+ (0.5 * Float64(dt)) .* v2
    v3 = ShipMakie._velocity_at(u, p3)
    p4 = p_dbl .+ Float64(dt) .* v3
    v4 = ShipMakie._velocity_at(u, p4)
    dp = (
        Float64(dt) / 6 * (v1[1] + 2 * v2[1] + 2 * v3[1] + v4[1]),
        Float64(dt) / 6 * (v1[2] + 2 * v2[2] + 2 * v3[2] + v4[2]),
        Float64(dt) / 6 * (v1[3] + 2 * v2[3] + 2 * v3[3] + v4[3]),
    )
    p_new = (Float32(p_dbl[1] + dp[1]),
             Float32(p_dbl[2] + dp[2]),
             Float32(p_dbl[3] + dp[3]))
    speed = Float32(sqrt(v1[1]^2 + v1[2]^2 + v1[3]^2))
    # Wrap or recycle particles that leave the domain
    if p_new[1] > Float32(NX - 3) || p_new[2] < 1 || p_new[2] > Float32(NY - 2) ||
       p_new[3] < 1 || p_new[3] > Float32(NZ - 2)
        # Recycle: re-seed upstream of the bow
        y = hull_yc - B_c * 1.4f0 + 2f0 * (B_c * 1.4f0) * rand(Float32)
        z = hull_zc - T_c * 1.1f0 + 2f0 * (T_c * 1.1f0) * rand(Float32)
        particles[idx] = (Float32(x_seed), y, z)
        return 0f0
    end
    particles[idx] = p_new
    return speed
end

# Hull SDF for the static hull mesh
hull_sdf_world(p) = ShipShapes.wigley_sdf(
    SVector(p[1] - hull_xc, p[2] - hull_yc, p[3] - hull_zc),
    Float64(L_c), Float64(B_c), Float64(T_c))

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "particle_tracer", "frames"))
mkpath(OUTDIR)
@info "Recording $NFRAMES frames into $OUTDIR …"
t0 = time()

# Speed colour per particle (updated each frame)
speeds = fill(0f0, NPART)

for frame in 1:NFRAMES
    # 1. advance the simulation one step
    WaterLily.mom_step!(sim.flow, sim.pois; udf = vlm_udf,
        pois_tol = 1f-6, pois_itmx = 50)
    step_vof_mules!(vof, sim; dt = sim.flow.Δt[end-1], perdir = (2,))

    # 2. advance the particles
    for idx in 1:NPART
        speeds[idx] = advance_particle!(particles, idx, sim.flow.u, 0.5f0)
    end

    # 3. push current positions to the trail buffer
    for idx in 1:NPART
        trail_buffer[trail_head, idx] = particles[idx]
    end
    global trail_head = mod1(trail_head + 1, TRAIL)

    # 4. render
    fig = Figure(size = (1400, 900), backgroundcolor = :gray12)
    ax = Axis3(fig[1, 1];
        title       = @sprintf("particle tracer — frame %d/%d", frame, NFRAMES),
        titlecolor  = :white,
        aspect      = :data,
        azimuth     = 0.35π + 0.3π * sin(2π * frame / NFRAMES),
        elevation   = 0.20π,
        xlabelcolor = :white, ylabelcolor = :white, zlabelcolor = :white,
        xticklabelcolor = :white, yticklabelcolor = :white, zticklabelcolor = :white,
        backgroundcolor = :gray12,
    )
    # Hull (opaque, dark) — only render once per frame for now
    hullmesh3d!(ax, hull_sdf_world;
        grid_size = (NX, NY, NZ), color = :navy)
    # Rotor wireframe
    bladedrotor3d!(ax, 3, R_prop, 0.5 * R_prop;
        center = (prop_xc, prop_yc, prop_zc),
        rotor_axis = (1.0, 0.0, 0.0),
        color = :orange, linewidth = 3)
    # Particles
    pts = [Point3f(p[1], p[2], p[3]) for p in particles]
    Makie.scatter!(ax, pts;
        color = speeds, colormap = :plasma, colorrange = (0f0, 3f0),
        markersize = 7, transparency = false)
    # Trails: draw faint lines from each particle's trail buffer
    trail_pts = NTuple{3, Float32}[]
    for idx in 1:NPART
        for off in 0:TRAIL-1
            k = mod1(trail_head - 1 - off, TRAIL)
            push!(trail_pts, trail_buffer[k, idx])
        end
        push!(trail_pts, (NaN32, NaN32, NaN32))  # break between particles
    end
    Makie.lines!(ax, trail_pts;
        color = (:white, 0.30), linewidth = 0.7)

    save(joinpath(OUTDIR, @sprintf("frame_%05d.png", frame)), fig)
    if frame % 10 == 0
        @printf "frame %3d / %d  elapsed=%.1fs\n" frame NFRAMES (time() - t0)
        flush(stdout)
    end
end

@printf "Done: %d frames in %s\n" NFRAMES OUTDIR

gif_out = abspath(joinpath(@__DIR__, "..", "runs", "particle_tracer", "tracer.gif"))
try
    cmd = `convert -delay 8 -loop 0 -resize 50% $(joinpath(OUTDIR, "frame_*.png")) $gif_out`
    run(pipeline(cmd, devnull))
    @printf "GIF: %s (%.2f MB)\n" gif_out stat(gif_out).size / 1e6
catch err
    @info "Skipped GIF" err
end
