# --- Internal utilities ---------------------------------------------------

"""
    cell_centres(N::Int)

Return the cell-centre world coordinates for a `N`-cell axis.
WaterLily convention: cell `I` has centre at `I − 1.5` (one ghost cell
at each end, so interior `2:N−1` maps to world coords `0.5 : N − 2.5`).
"""
@inline cell_centres(N::Int) = (1:N) .- 1.5

"""
    interior_view(arr::AbstractArray{T, D})

Return a view of `arr` without the one-cell ghost layer on every axis.
"""
@inline function interior_view(arr::AbstractArray{T, D}) where {T, D}
    sz = size(arr)
    @inbounds view(arr, ntuple(d -> 2:sz[d]-1, D)...)
end

"""
    eta_from_alpha(α::AbstractArray{T,3}, waterline_z::Real;
                   mask_inside::Union{Nothing, Function}=nothing) -> Matrix

Extract the free-surface elevation `η(x, y)` from a VoF colour function
`α`. For each `(i, j)` column, finds the highest `k` where
`α[i, j, k] ≥ 0.5` and the next cell up has `α < 0.5`; linearly
interpolates the crossing.

Returns a 2D `Matrix{T}` with `NaN` where the column has no surface
(fully dry or fully wet).

`mask_inside(i, j)` (optional, both args in cell index space) can be
passed to NaN-out columns inside the hull silhouette for cleaner plots.
"""
function eta_from_alpha(α::AbstractArray{T, 3}, waterline_z::Real;
                        mask_inside::Union{Nothing, Function} = nothing) where T
    nx, ny, nz = size(α)
    η = fill(T(NaN), nx, ny)
    @inbounds for j in 1:ny, i in 1:nx
        mask_inside !== nothing && mask_inside(i, j) && continue
        for k in (nz - 1):-1:2
            if α[i, j, k] ≥ T(0.5) && α[i, j, k + 1] < T(0.5)
                Δα = α[i, j, k] - α[i, j, k + 1]
                t  = abs(Δα) > T(1e-9) ? (T(0.5) - α[i, j, k + 1]) / Δα : T(0.5)
                z  = T(k + 1 - 1.5) - t
                η[i, j] = z - T(waterline_z)
                break
            end
        end
    end
    return η
end
