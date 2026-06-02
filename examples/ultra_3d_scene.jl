#!/usr/bin/env julia
#
# "Ultra" 3D scene: ship hull mesh + free surface + vortex volume +
# rotor wireframe, all in one Axis3, rotating camera.
#
# Single-axis showcase for the full ShipMakie recipe set in 3D.

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
const NFRAMES = parse(Int, get(ENV, "NFRAMES", "60"))
const BURNIN  = parse(Int, get(ENV, "BURNIN", "12"))
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
const hull_xc = Float32(NX / 5)
const hull_yc = Float32(NY / 2)
const hull_zc = H_w_c

const R_prop  = 2.5f0
const prop_xc = Float32(hull_xc + L_c / 2 + T_c)
const prop_yc = hull_yc
const prop_zc = Float32(H_w_c - T_c / 2)
const J_op = 0.32
const Ω    = Float64(π) * U∞ / (J_op * R_prop)

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1f0 : 0f0

rotor  = BladedRotor(; N_blades = 3, R = R_prop, R_hub = 0.5,
    chord = (0.8, 0.4), twist = (deg2rad(35.0), deg2rad(15.0)),
    ns = 12, nc = 4)
r_rot  = rotor_forces(rotor, U∞, Ω)
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

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "ultra_3d_scene", "frames"))
mkpath(OUTDIR)

mask_inside(i, j) =
    abs(Float32(i - 1.5) - hull_xc) ≤ L_c / 2 &&
    abs(Float32(j - 1.5) - hull_yc) ≤ B_c / 2

hull_sdf_world(p) = ShipShapes.wigley_sdf(
    SVector(p[1] - hull_xc, p[2] - hull_yc, p[3] - hull_zc),
    Float64(L_c), Float64(B_c), Float64(T_c))

@info "Rendering $NFRAMES frames into $OUTDIR …"
t0 = time()
for frame in 1:NFRAMES
    WaterLily.mom_step!(sim.flow, sim.pois; udf = vlm_udf,
        pois_tol = 1f-6, pois_itmx = 50)
    step_vof_mules!(vof, sim; dt = sim.flow.Δt[end-1], perdir = (2,))

    # Rotating camera: az sweeps over the run
    az = 0.30π + 0.40π * sin(2π * frame / NFRAMES)
    el = 0.15π + 0.10π * sin(2π * frame / NFRAMES * 2)

    fig = Figure(size = (1400, 1000), backgroundcolor = :gray12)
    ax = Axis3(fig[1, 1];
        title       = (@sprintf "ShipMakie • frame %d/%d (Wigley + rotor + vortices)" frame NFRAMES),
        titlecolor  = :white,
        aspect      = :data,
        azimuth     = az, elevation = el,
        xlabelcolor = :white, ylabelcolor = :white, zlabelcolor = :white,
        xticklabelcolor = :white, yticklabelcolor = :white,
        zticklabelcolor = :white,
        xtickcolor  = :white, ytickcolor = :white, ztickcolor = :white,
        backgroundcolor = :gray12,
    )
    # 1. hull mesh (deep blue, opaque)
    hullmesh3d!(ax, hull_sdf_world;
        grid_size = (NX, NY, NZ), level = 0.0,
        color = :navy)
    # 2. free surface (translucent blue tint)
    freesurface3d!(ax, vof.α;
        waterline_z = H_w_c, mask = mask_inside,
        colorrange = (-0.4, 0.4), alpha = 0.55)
    # 3. ω_mag volume (semi-transparent)
    vorticityvolume!(ax, sim.flow.u;
        field = :omega_mag, algorithm = :mip,
        colormap = :algae, absorption = 5.0)
    # 4. rotor wireframe (red so it pops against navy/algae)
    bladedrotor3d!(ax, 3, R_prop, 0.5 * R_prop;
        center = (prop_xc, prop_yc, prop_zc),
        rotor_axis = (1.0, 0.0, 0.0),
        color = :red, linewidth = 3)

    save(joinpath(OUTDIR, @sprintf("frame_%05d.png", frame)), fig)
    if frame % 5 == 0
        @printf "frame %3d / %d  elapsed=%.1fs\n" frame NFRAMES (time() - t0)
        flush(stdout)
    end
end
@printf "Done: %d frames\n" NFRAMES

gif_out = abspath(joinpath(@__DIR__, "..", "runs", "ultra_3d_scene", "ultra.gif"))
try
    cmd = `convert -delay 8 -loop 0 -resize 55% $(joinpath(OUTDIR, "frame_*.png")) $gif_out`
    run(pipeline(cmd, devnull))
    @printf "GIF: %s (%.2f MB)\n" gif_out stat(gif_out).size / 1e6
catch err
    @info "Skipped GIF" err
end
