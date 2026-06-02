#!/usr/bin/env julia
#
# Demo: moving ship with vortices, using only ShipMakie recipes.
#
# Runs a small integrated WaterLily + VoF + Wigley + BladedRotor +
# Rudder simulation and renders each frame with:
#   - Top-left:  η free-surface heatmap (etaheatmap + hullsilhouette)
#   - Top-right: centreline u_x slice (velocityslice + hullsilhouette)
#   - Bottom:    3D vorticity volume + bladedrotor3d wireframe
#
# Output: per-frame PNGs in runs/ship_vortices_demo/frames/.

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path = joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path = joinpath(@__DIR__, "..", "..", "VoF.jl")),
    Pkg.PackageSpec(path = joinpath(@__DIR__, "..", "..", "ShipShapes.jl")),
    Pkg.PackageSpec(path = joinpath(@__DIR__, "..", "..", "LiftingSurfaces.jl")),
    Pkg.PackageSpec(path = joinpath(@__DIR__, "..")),    # ShipMakie itself
]; io = devnull)
const USE_GL = get(ENV, "USE_GL", "1") != "0"
Pkg.add(["VortexLattice", USE_GL ? "GLMakie" : "CairoMakie"]; io = devnull)

using WaterLily, VoF, ShipShapes, LiftingSurfaces, ShipMakie
using ShipShapes: StaticArrays
const SVector = StaticArrays.SVector
using Printf
if USE_GL
    using GLMakie
    GLMakie.activate!()
else
    using CairoMakie
    CairoMakie.activate!()
end

# ---------------------------------------------------------------------------
# Simulation setup (kept deliberately small for fast demo)
# ---------------------------------------------------------------------------
const NX, NY, NZ = 96, 48, 32
const NFRAMES = parse(Int, get(ENV, "NFRAMES", "40"))
const BURNIN  = parse(Int, get(ENV, "BURNIN", "10"))
const L_c = 28f0; const B_c = 6f0; const T_c = 4f0
const ρ_w = 10f0; const ρ_a = 1f0
const U∞  = 1f0
const Fr  = 0.30f0
const G_c = U∞^2 / (Fr^2 * L_c)
const Re  = 5000f0
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
const rud_xc = Float32(prop_xc + R_prop * 0.5)
const rud_yc = hull_yc
const rud_zc = Float32(H_w_c - T_c / 2)

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1f0 : 0f0

rotor  = BladedRotor(; N_blades = 3, R = R_prop, R_hub = 0.5,
    chord = (0.8, 0.4), twist = (deg2rad(35.0), deg2rad(15.0)),
    ns = 12, nc = 4)
rudder = Rudder(; chord = 3.0, span = 4.0, ns = 12, nc = 6)
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
    r_rud = rudder_forces(rudder, deg2rad(8.0), U∞;
        inflow = trilinear_inflow(flow.u))
    side  = r_rud.CL * 0.5 * U∞^2 * (rudder.chord * rudder.span)
    drag  = r_rud.CD * 0.5 * U∞^2 * (rudder.chord * rudder.span)
    smear_force!(flow.f, SVector(-Float32(drag), Float32(side), 0f0),
        SVector(rud_xc, rud_yc, rud_zc); ε = 2.0f0)
    return nothing
end

# ---------------------------------------------------------------------------
# Burn-in
# ---------------------------------------------------------------------------
@info "Burn-in $BURNIN steps…"
for _ in 1:BURNIN
    WaterLily.mom_step!(sim.flow, sim.pois; udf = vlm_udf,
        pois_tol = 1f-6, pois_itmx = 50)
    step_vof_mules!(vof, sim; dt = sim.flow.Δt[end-1], perdir = (2,))
end

# ---------------------------------------------------------------------------
# Frame rendering
# ---------------------------------------------------------------------------
OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "ship_vortices_demo", "frames"))
mkpath(OUTDIR)
@info "Rendering $NFRAMES frames into $OUTDIR …"

# Plan-view hull box (for overlays)
hull_box = (hull_xc, hull_yc, L_c, B_c)
# Mask function: hide the hull's plan-view footprint on the η heatmap
mask_inside(i, j) =
    abs(Float32(i - 1.5) - hull_xc) ≤ L_c/2 &&
    abs(Float32(j - 1.5) - hull_yc) ≤ B_c/2

# Closure for the SDF in world frame (for hullsilhouette)
hull_sdf_world(p) = ShipShapes.wigley_sdf(SVector(
        p[1] - hull_xc, p[2] - hull_yc, p[3] - hull_zc),
    Float64(L_c), Float64(B_c), Float64(T_c))

t0 = time()
for frame in 1:NFRAMES
    WaterLily.mom_step!(sim.flow, sim.pois; udf = vlm_udf,
        pois_tol = 1f-6, pois_itmx = 50)
    step_vof_mules!(vof, sim; dt = sim.flow.Δt[end-1], perdir = (2,))

    fig = Figure(size = (1300, 850), backgroundcolor = :white)

    # ---- Panel 1: free-surface η heatmap, plan view --------------------
    ax_eta = Axis(fig[1, 1]; aspect = DataAspect(),
        title = @sprintf("frame %d: free-surface η  (top view)", frame),
        xlabel = "x (cells)", ylabel = "y (cells)")
    etaheatmap!(ax_eta, vof.α;
        waterline_z = H_w_c,
        mask        = mask_inside,
        hull_box    = hull_box,
        colorrange  = (-0.4, 0.4),
    )
    # Rotor (orange) and rudder (lime) markers
    Makie.scatter!(ax_eta, [Point2f(prop_xc, prop_yc)];
        color = :orange, marker = :cross, markersize = 16, strokewidth = 2)
    Makie.scatter!(ax_eta, [Point2f(rud_xc, rud_yc)];
        color = :lime, marker = :dtriangle, markersize = 14, strokewidth = 2)

    # ---- Panel 2: centreline u_x slice, side view ----------------------
    ax_side = Axis(fig[1, 2]; aspect = DataAspect(),
        title = "centreline u_x  (side view)",
        xlabel = "x (cells)", ylabel = "z (cells)")
    velocityslice!(ax_side, sim.flow.u;
        slice_axis = :y, index = NY ÷ 2, component = 1,
        colorrange = (-1.5, 3.0))
    # Overlay hull silhouette on the side slice
    hullsilhouette!(ax_side, hull_sdf_world;
        grid_size = (NX, NY, NZ), slice_axis = :y, slice_index = NY ÷ 2)
    Makie.scatter!(ax_side, [Point2f(prop_xc, prop_zc)];
        color = :orange, marker = :cross, markersize = 16, strokewidth = 2)
    Makie.scatter!(ax_side, [Point2f(rud_xc, rud_zc)];
        color = :lime, marker = :dtriangle, markersize = 14, strokewidth = 2)

    # ---- Panel 3: 3D vorticity volume + rotor wireframe ----------------
    ax_3d = Axis3(fig[2, 1:2];
        title = "vorticity volume (ω_mag, MIP) + rotor",
        xlabel = "x", ylabel = "y", zlabel = "z",
        aspect = :data)
    vorticityvolume!(ax_3d, sim.flow.u;
        field = :omega_mag, algorithm = :mip,
        colormap = :algae)
    bladedrotor3d!(ax_3d, 3, R_prop, 0.5 * R_prop;
        center = (prop_xc, prop_yc, prop_zc),
        rotor_axis = (1.0, 0.0, 0.0),
        color = :black, linewidth = 2.5)

    save(joinpath(OUTDIR, @sprintf("frame_%05d.png", frame)), fig)
    if frame % 5 == 0
        @printf "frame %3d / %d  |u|=%.3f  elapsed=%.1fs\n" frame NFRAMES maximum(abs, sim.flow.u) (time() - t0)
        flush(stdout)
    end
end

@printf "Done: %d frames written to %s\n" NFRAMES OUTDIR

# Assemble GIF via ImageMagick (best-effort; skip if convert missing)
gif_out = abspath(joinpath(@__DIR__, "..", "runs", "ship_vortices_demo", "ship_vortices.gif"))
try
    cmd = `convert -delay 10 -loop 0 -resize 60% $(joinpath(OUTDIR, "frame_*.png")) $gif_out`
    run(pipeline(cmd, devnull))
    @printf "Assembled GIF: %s (%.2f MB)\n" gif_out stat(gif_out).size / 1e6
catch err
    @info "ImageMagick `convert` not available — skipping GIF; frames are in $OUTDIR" err
end
