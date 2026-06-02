#!/usr/bin/env julia
#
# Kelvin wake demo: pure Wigley + free surface, no rotor. Renders
# the classic V-shaped wave pattern (half-angle ≈ 19.47°). Long
# domain in x so the wake has room to develop.

import Pkg
Pkg.activate(; temp=true)
Pkg.develop([
    Pkg.PackageSpec(path = joinpath(@__DIR__, "..", "..", "WaterLily")),
    Pkg.PackageSpec(path = joinpath(@__DIR__, "..", "..", "VoF.jl")),
    Pkg.PackageSpec(path = joinpath(@__DIR__, "..", "..", "ShipShapes.jl")),
    Pkg.PackageSpec(path = joinpath(@__DIR__, "..")),
]; io = devnull)
Pkg.add("CairoMakie"; io = devnull)

using WaterLily, VoF, ShipShapes, ShipMakie
using ShipShapes: StaticArrays
const SVector = StaticArrays.SVector
using Printf
using CairoMakie
CairoMakie.activate!()

# Long-in-x domain to show wake
const NX, NY, NZ = 240, 96, 32
const NSTEPS = parse(Int, get(ENV, "NSTEPS", "180"))
const L_c = 30f0; const B_c = 6f0; const T_c = 4f0
const ρ_w = 10f0; const ρ_a = 1f0
const U∞ = 1f0
const Fr = parse(Float32, get(ENV, "WL_FR", "0.30"))
const G_c = U∞^2 / (Fr^2 * L_c)
const Re = 5000f0
const ν_w_c = U∞ * L_c / Re
const μ_w_c = ρ_w * ν_w_c
const μ_a_c = ρ_a * 18 * ν_w_c
const H_w_c = Float32(NZ / 2)
const hull_xc = Float32(NX / 4); const hull_yc = Float32(NY / 2); const hull_zc = H_w_c

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1f0 : 0f0

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

@info "Running $NSTEPS steps at $NX×$NY×$NZ, Fr=$Fr…"
t0 = time()
for s in 1:NSTEPS
    WaterLily.mom_step!(sim.flow, sim.pois; pois_tol = 1f-6, pois_itmx = 50)
    step_vof_mules!(vof, sim; dt = sim.flow.Δt[end-1], perdir = (2,))
    if s % 30 == 0
        @printf "step %3d/%d  |u|=%.3f  elapsed=%.1fs\n" s NSTEPS maximum(abs, sim.flow.u) (time() - t0)
    end
end

mask_inside(i, j) =
    abs(Float32(i - 1.5) - hull_xc) ≤ L_c / 2 &&
    abs(Float32(j - 1.5) - hull_yc) ≤ B_c / 2

# Plan-view figure showing the wake
fig = Figure(size = (1800, 700), backgroundcolor = :white)
ax = Axis(fig[1, 1]; aspect = DataAspect(),
    title = @sprintf("Kelvin wake at Fr = %.2f  (classical half-angle = 19.47°)", Fr),
    xlabel = "x (cells)", ylabel = "y (cells)")
etaheatmap!(ax, vof.α;
    waterline_z = H_w_c, mask = mask_inside,
    hull_box = (hull_xc, hull_yc, L_c, B_c),
    colorrange = (-0.3, 0.3))

# Overlay the analytical Kelvin half-angle lines emanating from the bow
θ_kelvin = deg2rad(19.47)
x_bow = hull_xc + L_c / 2
x_tail = NX - 2.0
dx = x_tail - x_bow
y_off = dx * tan(θ_kelvin)
Makie.lines!(ax,
    [x_bow, x_tail], [hull_yc, hull_yc + y_off];
    color = :red, linestyle = :dash, linewidth = 1.5)
Makie.lines!(ax,
    [x_bow, x_tail], [hull_yc, hull_yc - y_off];
    color = :red, linestyle = :dash, linewidth = 1.5)
Makie.text!(ax, x_bow + 0.4 * dx, hull_yc + 0.5 * y_off;
    text = "19.47°", color = :red, fontsize = 14)

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "kelvin_wake"))
mkpath(OUTDIR)
out = joinpath(OUTDIR, @sprintf("kelvin_Fr%.2f.png", Fr))
save(out, fig)
@printf "Saved %s\n" out
