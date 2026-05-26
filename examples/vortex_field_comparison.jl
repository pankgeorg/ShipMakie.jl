#!/usr/bin/env julia
#
# Vortex-field comparison: ‖ω‖ vs λ₂ vs Q-criterion on the same flow.
# Three Axis3 side-by-side, each rendered with the corresponding
# vorticityvolume `field` option. Different methods highlight
# different vortical structures; this demo lets you see at a glance
# which one suits your data.

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

@info "Running $NSTEPS steps to develop wake…"
for _ in 1:NSTEPS
    WaterLily.mom_step!(sim.flow, sim.pois; udf = vlm_udf,
        pois_tol = 1f-6, pois_itmx = 50)
    step_vof_mules!(vof, sim; dt = sim.flow.Δt[end-1], perdir = (2,))
end

hull_sdf_world(p) = ShipShapes.wigley_sdf(
    SVector(p[1] - hull_xc, p[2] - hull_yc, p[3] - hull_zc),
    Float64(L_c), Float64(B_c), Float64(T_c))

# Three-panel figure: each panel uses a different vortex-detection field.
fig = Figure(size = (1800, 700), backgroundcolor = :gray10)

function make_axis(parent_pos, title)
    Axis3(parent_pos;
        title = title, titlecolor = :white, titlesize = 16,
        aspect = :data, azimuth = 0.40π, elevation = 0.18π,
        xlabelcolor = :white, ylabelcolor = :white, zlabelcolor = :white,
        xticklabelcolor = :white, yticklabelcolor = :white, zticklabelcolor = :white,
        backgroundcolor = :gray10,
    )
end

# Panel 1: ‖ω‖ MIP
ax1 = make_axis(fig[1, 1], "‖∇ × u‖  (MIP)")
hullmesh3d!(ax1, hull_sdf_world; grid_size = (NX, NY, NZ), color = :navy)
vorticityvolume!(ax1, sim.flow.u;
    field = :omega_mag, algorithm = :mip,
    colormap = :algae, absorption = 4.0)
bladedrotor3d!(ax1, 3, R_prop, 0.5 * R_prop;
    center = (prop_xc, prop_yc, prop_zc), color = :orange, linewidth = 2)

# Panel 2: λ₂ iso
ax2 = make_axis(fig[1, 2], "λ₂  (iso, Jeong & Hussain 1995)")
hullmesh3d!(ax2, hull_sdf_world; grid_size = (NX, NY, NZ), color = :navy)
vorticityvolume!(ax2, sim.flow.u;
    field = :lambda2, algorithm = :iso,
    colormap = :magma, absorption = 3.0)
bladedrotor3d!(ax2, 3, R_prop, 0.5 * R_prop;
    center = (prop_xc, prop_yc, prop_zc), color = :orange, linewidth = 2)

# Panel 3: Q-criterion iso
ax3 = make_axis(fig[1, 3], "Q  (iso, Hunt-Wray-Moin 1988)")
hullmesh3d!(ax3, hull_sdf_world; grid_size = (NX, NY, NZ), color = :navy)
vorticityvolume!(ax3, sim.flow.u;
    field = :q_criterion, algorithm = :iso,
    colormap = :plasma, absorption = 3.0)
bladedrotor3d!(ax3, 3, R_prop, 0.5 * R_prop;
    center = (prop_xc, prop_yc, prop_zc), color = :orange, linewidth = 2)

Label(fig[0, 1:3];
    text = "Vortex-field comparison: same flow, three detection methods",
    color = :white, fontsize = 22)

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "vortex_field_comparison"))
mkpath(OUTDIR)
out = joinpath(OUTDIR, "comparison.png")
save(out, fig)
@printf "Saved %s\n" out
