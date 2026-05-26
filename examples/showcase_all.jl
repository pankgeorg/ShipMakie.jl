#!/usr/bin/env julia
#
# Showcase: every ShipMakie recipe in a single multi-panel figure.
# This is the "all recipes in one picture" demo, intended for the
# README and for visual regression spotting.

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
const rud_xc = Float32(prop_xc + R_prop * 0.5)
const rud_yc = hull_yc
const rud_zc = Float32(H_w_c - T_c / 2)

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1f0 : 0f0

rotor  = BladedRotor(; N_blades = 3, R = R_prop, R_hub = 0.5,
    chord = (0.8, 0.4), twist = (deg2rad(35.0), deg2rad(15.0)),
    ns = 12, nc = 4)
rudder_chord, rudder_span = 3.0, 4.0

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

@info "Running $NSTEPS steps…"
for _ in 1:NSTEPS
    WaterLily.mom_step!(sim.flow, sim.pois; udf = vlm_udf,
        pois_tol = 1f-6, pois_itmx = 50)
    step_vof_mules!(vof, sim; dt = sim.flow.Δt[end-1], perdir = (2,))
end

hull_sdf_world(p) = ShipShapes.wigley_sdf(
    SVector(p[1] - hull_xc, p[2] - hull_yc, p[3] - hull_zc),
    Float64(L_c), Float64(B_c), Float64(T_c))
mask_inside(i, j) =
    abs(Float32(i - 1.5) - hull_xc) ≤ L_c / 2 &&
    abs(Float32(j - 1.5) - hull_yc) ≤ B_c / 2
hull_box = (hull_xc, hull_yc, L_c, B_c)

# ---------------------------------------------------------------------------
# Showcase figure: 6 panels covering every recipe
# ---------------------------------------------------------------------------
fig = Figure(size = (1800, 1200), backgroundcolor = :white)

# 1. etaheatmap (top-left)
ax1 = Axis(fig[1, 1]; aspect = DataAspect(),
    title = "etaheatmap (plan view)", xlabel = "x", ylabel = "y")
etaheatmap!(ax1, vof.α;
    waterline_z = H_w_c, mask = mask_inside,
    hull_box = hull_box, colorrange = (-0.4, 0.4))

# 2. velocityslice (top-mid)
ax2 = Axis(fig[1, 2]; aspect = DataAspect(),
    title = "velocityslice (u_x, side)", xlabel = "x", ylabel = "z")
velocityslice!(ax2, sim.flow.u;
    slice_axis = :y, index = NY ÷ 2, component = 1,
    colorrange = (-1.5, 3.0))
hullsilhouette!(ax2, hull_sdf_world;
    grid_size = (NX, NY, NZ), slice_axis = :y, slice_index = NY ÷ 2,
    color = :black, linewidth = 1.5)

# 3. streamslice (top-right)
ax3 = Axis(fig[1, 3]; aspect = DataAspect(),
    title = "streamslice (centreline)", xlabel = "x", ylabel = "z")
streamslice!(ax3, sim.flow.u;
    slice_axis = :y, index = NY ÷ 2, density = 0.7)
hullsilhouette!(ax3, hull_sdf_world;
    grid_size = (NX, NY, NZ), slice_axis = :y, slice_index = NY ÷ 2)

# 4. probeline (mid-left)
ax4 = Axis(fig[2, 1];
    title = "probeline: centreline α(z) at midship",
    xlabel = "z", ylabel = "α")
probeline!(ax4, vof.α;
    along = :z, fixed_indices = (round(Int, hull_xc), round(Int, hull_yc)),
    color = :steelblue)

# 5. 3D ultra-scene (mid-mid + mid-right)
ax5 = Axis3(fig[2, 2:3];
    title = "3D: hullmesh3d + freesurface3d + vorticityvolume + bladedrotor3d + rudder3d",
    aspect = :data, azimuth = 0.40π, elevation = 0.20π)
hullmesh3d!(ax5, hull_sdf_world;
    grid_size = (NX, NY, NZ), color = :navy)
freesurface3d!(ax5, vof.α;
    waterline_z = H_w_c, mask = mask_inside,
    colorrange = (-0.4, 0.4), alpha = 0.55)
vorticityvolume!(ax5, sim.flow.u;
    field = :omega_mag, algorithm = :mip,
    colormap = :algae, absorption = 5.0)
bladedrotor3d!(ax5, 3, R_prop, 0.5 * R_prop;
    center = (prop_xc, prop_yc, prop_zc),
    rotor_axis = (1.0, 0.0, 0.0),
    color = :orange, linewidth = 3)
rudder3d!(ax5, rudder_chord, rudder_span;
    center = (rud_xc, rud_yc, rud_zc),
    rudder_axis = (0.0, 0.0, 1.0),
    δ = deg2rad(8.0),
    color = :lime, linewidth = 3)

# 6. pressureisosurface (bot-left)
ax6 = Axis3(fig[3, 1];
    title = "pressureisosurface",
    aspect = :data, azimuth = 0.40π, elevation = 0.20π)
pressureisosurface!(ax6, Array(sim.flow.p))

# 7. vorticityvolume λ₂ (bot-mid)
ax7 = Axis3(fig[3, 2];
    title = "vorticityvolume (λ₂, iso)",
    aspect = :data, azimuth = 0.40π, elevation = 0.20π)
vorticityvolume!(ax7, sim.flow.u;
    field = :lambda2, algorithm = :iso,
    colormap = :magma, absorption = 3.0)

# 8. bladedrotor3d + rudder3d isolated (bot-right)
ax8 = Axis3(fig[3, 3];
    title = "bladedrotor3d + rudder3d",
    aspect = :data, azimuth = 0.40π, elevation = 0.20π)
bladedrotor3d!(ax8, 3, R_prop, 0.5 * R_prop;
    center = (prop_xc, prop_yc, prop_zc),
    rotor_axis = (1.0, 0.0, 0.0),
    color = :darkorange, linewidth = 3.5)
rudder3d!(ax8, rudder_chord, rudder_span;
    center = (rud_xc, rud_yc, rud_zc),
    rudder_axis = (0.0, 0.0, 1.0),
    δ = deg2rad(20.0),
    color = :forestgreen, linewidth = 3.5)

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "showcase_all"))
mkpath(OUTDIR)
out = joinpath(OUTDIR, "showcase.png")
save(out, fig)
@printf "Showcase saved to %s\n" out
