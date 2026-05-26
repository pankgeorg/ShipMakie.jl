"""
    vorticityvolume(u::AbstractArray{T, 4}; field = :omega_mag, …)

Plot a 3D volume rendering of a vorticity-related field derived from
a WaterLily face-staggered velocity array `u::Array{T, 4}` (size
`(NX, NY, NZ, 3)`).

# Attributes

- `field` — which scalar to render:
    - `:omega_mag` (default) — `‖∇ × u‖`, classic vortex visual
    - `:lambda2`             — Jeong & Hussain 1995 vortex-core
                               criterion (the second eigenvalue of
                               `S² + Ω²`, negated; positive lobes
                               are vortex cores)
    - `:q_criterion`         — Hunt, Wray & Moin 1988
                               `Q = ½(‖Ω‖² − ‖S‖²)`; positive
                               indicates rotation-dominated regions.
- `algorithm` — `:iso` (default) or `:mip` (maximum intensity).
- `isovalue` — for `:iso`, the iso-level. Default 0.3 of the field's
                95th percentile.
- `isorange` — `:iso` thickness. Default 0.05·isovalue.
- `colormap` — default `:algae` (matches WaterLily-Examples).

# Notes

The vorticity is computed at cell centres via central differences.
For the staggered velocity, `∂u_i/∂x_j` at cell centre uses 4-point
averages over the two relevant face-staggered cell pairs. We compute
this in-place into a workspace array of shape `(NX, NY, NZ)`.

# Example

```julia
vorticityvolume(sim.flow.u; field = :omega_mag, algorithm = :iso)
```
"""
Makie.@recipe(VorticityVolume, u) do scene
    Makie.Attributes(
        field      = :omega_mag,
        algorithm  = :iso,
        isovalue   = Makie.automatic,
        isorange   = Makie.automatic,
        colormap   = :algae,
        colorrange = Makie.automatic,
        absorption = 4.0,
    )
end

# Cell-centred ∂u_i/∂x_j on a face-staggered grid (WaterLily convention).
# Stencil: face value at (I, i) lives at world position I[i] - 1.
# For ∂(i, j, I, u) at cell I, use central difference along axis j of
# the i-component, averaged over the two faces of cell I.
@inline function _dui_dxj(u, I, i, j)
    # Forward and back neighbours along j
    Iplus  = Base.setindex(I, I[j] + 1, j)
    Iminus = Base.setindex(I, I[j] - 1, j)
    return 0.5 * (u[Iplus..., i] - u[Iminus..., i])
end

# Compute ω = ∇ × u and its magnitude at each interior cell.
function _omega_mag!(out::Array{T, 3}, u::Array{T, 4}) where T
    nx, ny, nz, _ = size(u)
    @inbounds for k in 2:nz-1, j in 2:ny-1, i in 2:nx-1
        I = (i, j, k)
        # ω_x = ∂u_z/∂y - ∂u_y/∂z
        ωx = _dui_dxj(u, I, 3, 2) - _dui_dxj(u, I, 2, 3)
        ωy = _dui_dxj(u, I, 1, 3) - _dui_dxj(u, I, 3, 1)
        ωz = _dui_dxj(u, I, 2, 1) - _dui_dxj(u, I, 1, 2)
        out[i, j, k] = sqrt(ωx * ωx + ωy * ωy + ωz * ωz)
    end
    # Pad boundaries with zero (we don't compute there)
    out[1, :, :]  .= zero(T); out[end, :, :] .= zero(T)
    out[:, 1, :]  .= zero(T); out[:, end, :] .= zero(T)
    out[:, :, 1]  .= zero(T); out[:, :, end] .= zero(T)
    return out
end

# Q-criterion (Hunt, Wray & Moin 1988): Q = (Ω_ij·Ω_ij - S_ij·S_ij) / 2.
# Positive Q indicates regions where rotation dominates strain (vortex cores).
function _q_criterion!(out::Array{T, 3}, u::Array{T, 4}) where T
    nx, ny, nz, _ = size(u)
    @inbounds for k in 2:nz-1, j in 2:ny-1, i in 2:nx-1
        I = (i, j, k)
        g11 = _dui_dxj(u, I, 1, 1); g12 = _dui_dxj(u, I, 1, 2); g13 = _dui_dxj(u, I, 1, 3)
        g21 = _dui_dxj(u, I, 2, 1); g22 = _dui_dxj(u, I, 2, 2); g23 = _dui_dxj(u, I, 2, 3)
        g31 = _dui_dxj(u, I, 3, 1); g32 = _dui_dxj(u, I, 3, 2); g33 = _dui_dxj(u, I, 3, 3)
        # S, Ω components
        S11 = g11; S22 = g22; S33 = g33
        S12 = 0.5 * (g12 + g21); S13 = 0.5 * (g13 + g31); S23 = 0.5 * (g23 + g32)
        Ω12 = 0.5 * (g12 - g21); Ω13 = 0.5 * (g13 - g31); Ω23 = 0.5 * (g23 - g32)
        ΩΩ = Ω12*Ω12 + Ω13*Ω13 + Ω23*Ω23     # antisymmetric: only off-diag
        SS = S11*S11 + S22*S22 + S33*S33 +
             2 * (S12*S12 + S13*S13 + S23*S23)
        out[i, j, k] = ΩΩ - 0.5 * SS         # Q = (||Ω||² - ||S||²)/2 = ΩΩ - SS/2
    end
    out[1, :, :]  .= zero(T); out[end, :, :] .= zero(T)
    out[:, 1, :]  .= zero(T); out[:, end, :] .= zero(T)
    out[:, :, 1]  .= zero(T); out[:, :, end] .= zero(T)
    return out
end

# λ₂ vortex-core criterion (Jeong & Hussain 1995). Returns -λ₂ so that
# positive values indicate vortex cores (intuitive for iso rendering).
function _lambda2!(out::Array{T, 3}, u::Array{T, 4}) where T
    nx, ny, nz, _ = size(u)
    @inbounds for k in 2:nz-1, j in 2:ny-1, i in 2:nx-1
        I = (i, j, k)
        # Build velocity-gradient tensor
        g11 = _dui_dxj(u, I, 1, 1); g12 = _dui_dxj(u, I, 1, 2); g13 = _dui_dxj(u, I, 1, 3)
        g21 = _dui_dxj(u, I, 2, 1); g22 = _dui_dxj(u, I, 2, 2); g23 = _dui_dxj(u, I, 2, 3)
        g31 = _dui_dxj(u, I, 3, 1); g32 = _dui_dxj(u, I, 3, 2); g33 = _dui_dxj(u, I, 3, 3)
        # Symmetric S and antisymmetric Ω
        S11 = g11; S22 = g22; S33 = g33
        S12 = 0.5 * (g12 + g21); S13 = 0.5 * (g13 + g31); S23 = 0.5 * (g23 + g32)
        Ω12 = 0.5 * (g12 - g21); Ω13 = 0.5 * (g13 - g31); Ω23 = 0.5 * (g23 - g32)
        # M = S² + Ω² (symmetric)
        M11 = S11*S11 + S12*S12 + S13*S13 - Ω12*Ω12 - Ω13*Ω13
        M22 = S12*S12 + S22*S22 + S23*S23 - Ω12*Ω12 - Ω23*Ω23
        M33 = S13*S13 + S23*S23 + S33*S33 - Ω13*Ω13 - Ω23*Ω23
        M12 = S11*S12 + S12*S22 + S13*S23 + Ω13*Ω23
        M13 = S11*S13 + S12*S23 + S13*S33 - Ω12*Ω23
        M23 = S12*S13 + S22*S23 + S23*S33 + Ω12*Ω13
        # Eigenvalues of a 3×3 symmetric matrix — middle eigenvalue.
        # Closed form via the characteristic polynomial.
        tr   = M11 + M22 + M33
        I2   = M11*M22 + M22*M33 + M33*M11 - M12*M12 - M13*M13 - M23*M23
        I3   = M11*(M22*M33 - M23*M23) - M12*(M12*M33 - M23*M13) +
               M13*(M12*M23 - M22*M13)
        # Cardano (Wikipedia: Eigenvalues of symmetric 3×3 matrix)
        q = tr / 3
        p2 = (M11 - q)^2 + (M22 - q)^2 + (M33 - q)^2 +
             2 * (M12^2 + M13^2 + M23^2)
        p_ = sqrt(p2 / 6)
        if p_ < eps(T)
            λ2 = q
        else
            B11 = (M11 - q) / p_; B22 = (M22 - q) / p_; B33 = (M33 - q) / p_
            B12 = M12 / p_; B13 = M13 / p_; B23 = M23 / p_
            detB = B11*(B22*B33 - B23*B23) - B12*(B12*B33 - B23*B13) +
                   B13*(B12*B23 - B22*B13)
            r = clamp(detB / 2, -one(T), one(T))
            φ = acos(r) / 3
            λ2 = q + 2 * p_ * cos(φ + 2π / 3)  # middle eigenvalue
        end
        out[i, j, k] = -λ2  # positive in vortex cores
    end
    out[1, :, :]  .= zero(T); out[end, :, :] .= zero(T)
    out[:, 1, :]  .= zero(T); out[:, end, :] .= zero(T)
    out[:, :, 1]  .= zero(T); out[:, :, end] .= zero(T)
    return out
end

function Makie.plot!(p::VorticityVolume)
    u_obs = p[1]
    field_obs = Makie.lift(u_obs, p.field) do u, fld
        out = zeros(eltype(u), size(u, 1), size(u, 2), size(u, 3))
        if fld === :omega_mag
            _omega_mag!(out, u)
        elseif fld === :lambda2
            _lambda2!(out, u)
        elseif fld === :q_criterion
            _q_criterion!(out, u)
        else
            throw(ArgumentError(
                "field must be :omega_mag, :lambda2, or :q_criterion; got $(repr(fld))"))
        end
        return out
    end
    # Default isovalue: 30% of the 95th percentile of the field.
    iso_obs = Makie.lift(field_obs, p.isovalue) do f, iv
        iv === Makie.automatic ? sort_quantile(f, 0.95) * 0.3 : iv
    end
    range_obs = Makie.lift(iso_obs, p.isorange) do iv, ir
        ir === Makie.automatic ? iv * 0.05 : ir
    end
    Makie.volume!(p, field_obs;
        algorithm  = p.algorithm,
        isovalue   = iso_obs,
        isorange   = range_obs,
        colormap   = p.colormap,
        colorrange = p.colorrange,
        absorption = p.absorption,
    )
    return p
end

# Cheap quantile (avoids sorting the full array). Linear scan for the
# 95th percentile threshold using a partial sort approach.
@inline function sort_quantile(a::AbstractArray, q::Real)
    v = vec(a)
    n = length(v)
    k = max(1, Int(round(n * q)))
    return partialsort(copy(v), k)
end
