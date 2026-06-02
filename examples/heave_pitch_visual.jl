#!/usr/bin/env julia
#
# Heave + pitch 2-DOF visualization: actually show the moving ship.
# The body rotates and translates each frame; hullmesh3d picks up
# the new SDF, freesurface3d picks up the new α field. Companion to
# RESULTS-heave-pitch-containership.md (U1).

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path = joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path = joinpath(@__DIR__, "..", "..", "VoF.jl")),
    Pkg.PackageSpec(path = joinpath(@__DIR__, "..", "..", "ShipShapes.jl")),
    Pkg.PackageSpec(path = joinpath(@__DIR__, "..")),
]; io = devnull)
Pkg.add("GLMakie"; io = devnull)

using WaterLily, VoF, ShipShapes, ShipMakie
using ShipShapes: StaticArrays
const SVector = StaticArrays.SVector
using Printf
using GLMakie
GLMakie.activate!()

# Containership params (per U1 — converges nicely)
const NX, NY, NZ = 128, 64, 32
const NFRAMES = parse(Int, get(ENV, "NFRAMES", "60"))
const L_c = 36.0; const B_c = 8.0; const T_c = 5.0
const ρ_w = 10.0; const ρ_a = 1.0
const U∞  = 1.0
const Fr  = 0.25
const G_c = U∞^2 / (Fr^2 * L_c)
const Re  = 5000.0
const ν_w_c = U∞ * L_c / Re
const μ_w_c = ρ_w * ν_w_c
const μ_a_c = ρ_a * 18 * ν_w_c
const H_w_c = NZ/2
const hull_xc = NX/5; const hull_yc = NY/2; const hull_zc0 = H_w_c

const V0      = 0.75 * L_c * B_c * T_c
const M_ship  = ρ_w * V0
const I_pitch = M_ship * L_c^2 / 12.0

z_h    = Ref(0.0); zdot_h = Ref(0.0)
θ      = Ref(0.0); θdot   = Ref(0.0)

hull_map = (x, t) -> begin
    dx = x[1] - hull_xc
    dz = x[3] - (hull_zc0 + z_h[])
    cθ, sθ = cos(θ[]), sin(θ[])
    SVector(cθ * dx + sθ * dz, x[2] - hull_yc, -sθ * dx + cθ * dz)
end
hull = ShipShapes.Containership(; L = L_c, B = B_c, T = T_c,
    par_frac = 0.5, deck_h = T_c / 2, map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1.0 : 0.0

vof = VoFFlow((NX, NY, NZ);
    α₀ = α₀, ρ_w = ρ_w, ρ_a = ρ_a, μ_w = μ_w_c, μ_a = μ_a_c, T = Float64)
function vof_pois_ctor(flow)
    L = similar(flow.μ₀)
    @inbounds for I in CartesianIndices(L)
        L[I] = flow.μ₀[I] * vof.L[I]
    end
    WaterLily.MultiLevelPoisson(flow.p, L, flow.σ; perdir = (2,))
end
sim = WaterLily.Simulation((NX, NY, NZ), (U∞, 0.0, 0.0), L_c;
    T = Float64, ν = VoF.viscosity(vof),
    g = (i, x, t) -> i == 3 ? -G_c : 0.0,
    Δt = 0.25, body = hull, ϵ = 1, perdir = (2,), exitBC = true,
    pois_ctor = vof_pois_ctor, U = U∞,
)

# Containership SDF in world frame, time-varying through z_h, θ refs
function hull_sdf_world_dynamic(p)
    dx = p[1] - hull_xc
    dz = p[3] - (hull_zc0 + z_h[])
    cθ, sθ = cos(θ[]), sin(θ[])
    body_p = SVector(cθ * dx + sθ * dz, p[2] - hull_yc, -sθ * dx + cθ * dz)
    ShipShapes.containership_sdf(body_p, L_c, B_c, T_c, 0.5, T_c / 2)
end
mask_inside(i, j) =
    abs(Float64(i - 1.5) - hull_xc) ≤ L_c / 2 &&
    abs(Float64(j - 1.5) - hull_yc) ≤ B_c / 2

const GRAVITY_RAMP_STEPS = NFRAMES ÷ 4
x_CG = [hull_xc, hull_yc, hull_zc0]
OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "heave_pitch_visual"))
mkpath(joinpath(OUTDIR, "frames"))
@info "Rendering $NFRAMES frames of moving Containership…"
t0 = time()

for frame in 1:NFRAMES
    WaterLily.measure!(sim, sum(sim.flow.Δt))
    WaterLily.mom_step!(sim.flow, sim.pois; pois_tol = 1e-6, pois_itmx = 50)
    step_vof_mules!(vof, sim; dt = sim.flow.Δt[end-1], perdir = (2,))

    g_scale = min(1.0, frame / GRAVITY_RAMP_STEPS)
    g_now   = G_c * g_scale
    Fp = WaterLily.pressure_force(sim); Fv = WaterLily.viscous_force(sim)
    Mp = WaterLily.pressure_moment(x_CG, sim)
    F_buoy = -(Fp[3] + Fv[3])
    F_grav = -M_ship * g_now
    M_y_bdim = -Mp[2]
    K_pitch  = ρ_w * g_now * B_c * L_c^3 / 12   # Containership wp
    M_y      = M_y_bdim - K_pitch * θ[]
    β_h = 0.05 / sim.flow.Δt[end-1]
    β_p = 0.2  / sim.flow.Δt[end-1]
    z_ddot = (F_buoy + F_grav) / M_ship - β_h * zdot_h[]
    θ_ddot = M_y / I_pitch - β_p * θdot[]
    zdot_h[] += z_ddot * sim.flow.Δt[end-1]
    z_h[]    += zdot_h[] * sim.flow.Δt[end-1]
    θdot[]   += θ_ddot * sim.flow.Δt[end-1]
    θ[]      += θdot[] * sim.flow.Δt[end-1]
    z_h[] = clamp(z_h[], -T_c, T_c)
    θ[]   = clamp(θ[], -deg2rad(15), deg2rad(15))

    # Render
    fig = Figure(size = (1300, 850), backgroundcolor = :gray12)
    ax = Axis3(fig[1, 1];
        title = @sprintf("Heave + pitch 2-DOF  •  z_h=%+.2f cells  θ=%+.2f°  (frame %d/%d)",
                        z_h[], θ[]*180/π, frame, NFRAMES),
        titlecolor = :white,
        aspect = :data, azimuth = 0.40π, elevation = 0.20π,
        xlabelcolor = :white, ylabelcolor = :white, zlabelcolor = :white,
        xticklabelcolor = :white, yticklabelcolor = :white, zticklabelcolor = :white,
        backgroundcolor = :gray12,
    )
    hullmesh3d!(ax, hull_sdf_world_dynamic;
        grid_size = (NX, NY, NZ), color = :navy)
    freesurface3d!(ax, Float32.(vof.α);
        waterline_z = H_w_c, mask = mask_inside,
        colorrange = (-0.5, 0.5), alpha = 0.55)
    vorticityvolume!(ax, Float32.(sim.flow.u);
        field = :omega_mag, algorithm = :mip,
        colormap = :algae, absorption = 4.0)

    save(joinpath(OUTDIR, "frames", @sprintf("frame_%05d.png", frame)), fig)
    if frame % 5 == 0
        @printf "frame %3d / %d  z_h=%+.2f  θ=%+.2f°  elapsed=%.1fs\n" frame NFRAMES z_h[] θ[]*180/π (time()-t0)
        flush(stdout)
    end
end

gif_out = joinpath(OUTDIR, "heave_pitch.gif")
try
    cmd = `convert -delay 8 -loop 0 -resize 55% $(joinpath(OUTDIR, "frames", "frame_*.png")) $gif_out`
    run(pipeline(cmd, devnull))
    @printf "GIF: %s (%.2f MB)\n" gif_out stat(gif_out).size / 1e6
catch err
    @info "Skipped GIF" err
end
