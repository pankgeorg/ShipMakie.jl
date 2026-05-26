# --- Curated colormaps for ship-CFD --------------------------------------
# These are *named pointers* to existing colormaps; we don't ship new
# colour data. They encode the per-quantity convention used across the
# stack's RESULTS plots so a downstream caller doesn't have to remember
# which colormap to pair with which field.

"""
    ship_colormaps[]

Const-binding of per-quantity colormap names. Use like:

```julia
heatmap(η; colormap = ShipMakie.ship_colormaps[].eta)
```

| Key       | Quantity                | Colormap  | Notes |
|-----------|-------------------------|-----------|-------|
| `eta`     | Free-surface elevation  | `:RdBu`   | symmetric around 0 |
| `velocity`| Velocity component       | `:viridis`| 0…U∞ range  |
| `pressure`| Pressure deviation       | `:RdBu`   | symmetric around 0 |
| `vorticity`| Vorticity magnitude     | `:algae`  | matches WaterLily-Examples |
| `vof_alpha`| α colour function      | `:Blues`  | 0=air, 1=water |
"""
const ship_colormaps = Ref((
    eta       = :RdBu,
    velocity  = :viridis,
    pressure  = :RdBu,
    vorticity = :algae,
    vof_alpha = :Blues,
))
