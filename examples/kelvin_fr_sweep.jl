#!/usr/bin/env julia
#
# Kelvin Fr sweep: render the η heatmap at 4 different Fr values to
# show the Fr-invariance of the Kelvin half-angle (19.47°). Composite
# the 4 panels into a single figure for direct comparison.

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

const NX, NY, NZ = 200, 96, 32
const NSTEPS = parse(Int, get(ENV, "NSTEPS", "150"))
const FR_LIST = parse.(Float32, split(get(ENV, "WL_FR_LIST", "0.25,0.30,0.35,0.40"), ","))
const L_c = 30f0; const B_c = 6f0; const T_c = 4f0
const ρ_w = 10f0; const ρ_a = 1f0
const U∞ = 1f0
const Re = 5000f0
const ν_w_c = U∞ * L_c / Re
const μ_w_c = ρ_w * ν_w_c
const μ_a_c = ρ_a * 18 * ν_w_c
const H_w_c = Float32(NZ / 2)
const hull_xc = Float32(NX / 4); const hull_yc = Float32(NY / 2); const hull_zc = H_w_c

hull_map = (x, t) -> SVector(x[1] - hull_xc, x[2] - hull_yc, x[3] - hull_zc)
hull = ShipShapes.Wigley(; L = L_c, B = B_c, T = T_c, map = hull_map)
α₀(_i, x_cell) = (x_cell[3] ≤ H_w_c) ? 1f0 : 0f0
mask_inside(i, j) =
    abs(Float32(i - 1.5) - hull_xc) ≤ L_c / 2 &&
    abs(Float32(j - 1.5) - hull_yc) ≤ B_c / 2

function run_fr(Fr)
    G_c = U∞^2 / (Fr^2 * L_c)
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
    t0 = time()
    for _ in 1:NSTEPS
        WaterLily.mom_step!(sim.flow, sim.pois; pois_tol = 1f-6, pois_itmx = 50)
        step_vof_mules!(vof, sim; dt = sim.flow.Δt[end-1], perdir = (2,))
    end
    @printf "Fr=%.2f done in %.1fs\n" Fr (time() - t0)
    return vof.α
end

# Composite figure with 4 panels (one per Fr)
fig = Figure(size = (1500, 1200), backgroundcolor = :white)

for (idx, Fr) in enumerate(FR_LIST)
    @info "Running Fr=$Fr"
    α = run_fr(Fr)
    row = (idx - 1) ÷ 2 + 1
    col = (idx - 1) % 2 + 1
    ax = Axis(fig[row, col]; aspect = DataAspect(),
        title = @sprintf("Fr = %.2f", Fr),
        xlabel = "x (cells)", ylabel = "y (cells)")
    etaheatmap!(ax, α;
        waterline_z = H_w_c, mask = mask_inside,
        hull_box = (hull_xc, hull_yc, L_c, B_c),
        colorrange = (-0.4, 0.4))
    # Overlay Kelvin wedge for reference
    θ_k = deg2rad(19.47)
    x_bow = hull_xc + L_c / 2
    x_tail = NX - 2.0
    dx = x_tail - x_bow
    y_off = dx * tan(θ_k)
    Makie.lines!(ax, [x_bow, x_tail], [hull_yc, hull_yc + y_off];
        color = :red, linestyle = :dash, linewidth = 1.5)
    Makie.lines!(ax, [x_bow, x_tail], [hull_yc, hull_yc - y_off];
        color = :red, linestyle = :dash, linewidth = 1.5)
end

Label(fig[0, 1:2];
    text = "Kelvin wake — Fr sweep (Wigley, no rotor)\n" *
           "Red dashed lines: classical 19.47° half-angle (Fr-invariant)",
    fontsize = 16, halign = :center)

OUTDIR = abspath(joinpath(@__DIR__, "..", "runs", "kelvin_fr_sweep"))
mkpath(OUTDIR)
out = joinpath(OUTDIR, "fr_sweep.png")
save(out, fig)
@printf "Saved %s\n" out
