using Test
using ShipMakie
using CairoMakie
using StaticArrays

@testset "ShipMakie" begin

    @testset "etaheatmap on a synthetic α field" begin
        # Build a fake VoF α field with a sinusoidal free surface.
        NX, NY, NZ = 32, 16, 16
        α = zeros(Float32, NX, NY, NZ)
        for k in 1:NZ, j in 1:NY, i in 1:NX
            z = k - 1.5f0
            η_target = NZ/2 + 0.5f0 * sinpi(2f0 * i / NX)
            α[i, j, k] = z ≤ η_target ? 1f0 : 0f0
        end
        fig = Figure(size=(400, 200))
        ax  = Axis(fig[1, 1])
        etaheatmap!(ax, α; waterline_z = NZ/2)
        @test fig isa Figure
    end

    @testset "etaheatmap mask + hull_box" begin
        NX, NY, NZ = 32, 16, 16
        α = ones(Float32, NX, NY, NZ)
        α[:, :, NZ÷2+1:end] .= 0f0
        fig = Figure(size=(400, 200))
        ax  = Axis(fig[1, 1])
        mask = (i, j) -> 10 ≤ i ≤ 22 && 5 ≤ j ≤ 12
        etaheatmap!(ax, α; mask = mask, hull_box = (16.0, 8.0, 12.0, 8.0))
        @test fig isa Figure
    end

    @testset "velocityslice 3D u field" begin
        NX, NY, NZ = 16, 8, 8
        u = zeros(Float32, NX, NY, NZ, 3)
        u[:, :, :, 1] .= 1f0
        u[:, :, :, 3] .= 0.1f0
        for ax in (:x, :y, :z), comp in 1:3
            fig = Figure(size=(400, 200))
            ax_obj = Axis(fig[1, 1])
            velocityslice!(ax_obj, u; slice_axis = ax, component = comp)
            @test fig isa Figure
        end
    end

    @testset "hullsilhouette on a sphere SDF" begin
        sphere_sdf(p) = sqrt(p[1]^2 + p[2]^2 + p[3]^2) - 4.0
        fig = Figure(size=(400, 200))
        ax  = Axis(fig[1, 1])
        hullsilhouette!(ax, sphere_sdf;
            grid_size = (16, 16, 16), slice_axis = :z)
        @test fig isa Figure
    end

    @testset "probeline 1D / 2D / 3D" begin
        # 3D
        f3 = rand(Float32, 16, 8, 8)
        fig3 = Figure()
        ax3  = Axis(fig3[1, 1])
        probeline!(ax3, f3; along = :x, fixed_indices = (4, 4))
        @test fig3 isa Figure
        # 2D
        f2 = rand(Float32, 16, 8)
        fig2 = Figure()
        ax2  = Axis(fig2[1, 1])
        probeline!(ax2, f2; along = :y, fixed_indices = 8)
        @test fig2 isa Figure
        # 1D
        f1 = rand(Float32, 16)
        fig1 = Figure()
        ax1  = Axis(fig1[1, 1])
        probeline!(ax1, f1; along = :x)
        @test fig1 isa Figure
    end

    @testset "vorticityvolume omega_mag + lambda2" begin
        NX, NY, NZ = 16, 16, 16
        u = zeros(Float32, NX, NY, NZ, 3)
        # Add a rotational core around z-axis to produce vorticity
        for k in 1:NZ, j in 1:NY, i in 1:NX
            x = i - NX/2; y = j - NY/2
            r = sqrt(x*x + y*y)
            θ = atan(y, x)
            v = exp(-r/3)
            u[i, j, k, 1] = -v * sin(θ)
            u[i, j, k, 2] = +v * cos(θ)
        end
        # omega_mag
        fig = Figure(size=(400, 400))
        ax = Axis3(fig[1, 1])
        vorticityvolume!(ax, u; field = :omega_mag, algorithm = :mip)
        @test fig isa Figure
        # lambda2
        fig2 = Figure(size=(400, 400))
        ax2 = Axis3(fig2[1, 1])
        vorticityvolume!(ax2, u; field = :lambda2, algorithm = :mip)
        @test fig2 isa Figure
    end

    @testset "bladedrotor3d wireframe" begin
        fig = Figure(size=(400, 400))
        ax = Axis3(fig[1, 1])
        bladedrotor3d!(ax, 3, 2.0, 0.4;
            center = (10.0, 5.0, 5.0), rotor_axis = (1.0, 0.0, 0.0))
        @test fig isa Figure
    end

    @testset "freesurface3d builds a Figure" begin
        NX, NY, NZ = 16, 16, 16
        α = zeros(Float32, NX, NY, NZ)
        for k in 1:NZ, j in 1:NY, i in 1:NX
            z = k - 1.5f0
            η_target = NZ/2 + 0.3f0 * sinpi(2f0 * i / NX)
            α[i, j, k] = z ≤ η_target ? 1f0 : 0f0
        end
        fig = Figure(size=(400, 400))
        ax = Axis3(fig[1, 1])
        freesurface3d!(ax, α; waterline_z = NZ/2)
        @test fig isa Figure
    end

    @testset "streamslice builds a Figure" begin
        NX, NY, NZ = 16, 8, 8
        u = zeros(Float32, NX, NY, NZ, 3)
        u[:, :, :, 1] .= 1f0  # constant flow in +x
        fig = Figure(size=(400, 200))
        ax = Axis(fig[1, 1])
        streamslice!(ax, u; slice_axis = :y, density = 0.6)
        @test fig isa Figure
    end

    @testset "ShipMakie.eta_from_alpha lifts surface correctly" begin
        # A flat surface at z = 5
        NX, NY, NZ = 8, 8, 16
        α = zeros(Float32, NX, NY, NZ)
        for k in 1:NZ, j in 1:NY, i in 1:NX
            z = k - 1.5f0
            α[i, j, k] = z ≤ 5f0 ? 1f0 : 0f0
        end
        η = ShipMakie.eta_from_alpha(α, NZ/2)
        # All columns should report the same surface, ≈ -3 cells from
        # the midline at NZ/2 = 8 (5 - 8 = -3)
        ηfin = filter(isfinite, vec(η))
        @test all(abs.(ηfin .- (5f0 - NZ/2)) .< 0.55f0)
    end

end # @testset
