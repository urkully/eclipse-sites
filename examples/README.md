# Examples

Output from real runs, kept so you can see what the tool produces without
downloading a DEM first. The DEM these came from is not in the repository. It
is the PNOA 5 m LiDAR DTM for Asturias, several hundred megabytes, and
`terrain/manifest.json` records how it was merged.

## `candamo-pool.png`, `candamo-pool.json`

The pool at San Roman de Candamo, 43.4524508, -6.0705374. This is the
regression anchor the project is checked against, so it is the one output here
that has a known right answer: **+8.6 degrees of clearance** at totality, with
the horizon at +1.8 degrees formed by a 753 m massif 18.7 km away.

```bash
python -m eclipse_sites.cli site \
    --lat 43.4524508 --lon -6.0705374 --elev 96.7 \
    --dem terrain/dtm05.tif --range 40000 \
    --view examples/candamo-pool.png --json examples/candamo-pool.json
```

The elevation is bare terrain from the DTM, not a surface model and not phone
GPS. Both of those would be wrong here for reasons the main README explains.

Note `"ray_truncated": true` in the JSON. The ray ran out of raster at 37.8 km
of the 40 km asked for, so the last 2 km in that direction is unknown rather
than clear. The tool says so instead of quietly rounding it to a clear view.

## `cruz-de-candamo.png`

A second point about 350 m away, 43.45539, -6.07192, which comes out at +9.1
degrees. Useful mainly as a contrast: half a degree of clearance separates two
spots you could walk between in five minutes, which is the whole argument for
checking a coordinate rather than a village.

The labels collide slightly in this render where the totality marker meets the
20:15 tick. It is a cosmetic flaw in the renderer, not a bad result.

## `candamo-canopy.json`

A canopy profile measured by fitting a photograph to the terrain, which is what
`calibrate --out-canopy` writes. The source photograph is not in the
repository. It covers azimuths 235 to 323 degrees and kept 29 bins, dropping 20
others that could not be honestly converted to a height.

Read `meta` before trusting the numbers. It records the fitted camera at 265.22
degrees bearing and 14.65 degrees pitch, a residual of 0.447 degrees rms over
158 anchor columns, and exactly why each dropped bin was dropped. A profile
with a poor residual is not worth using, because canopy inherits every error
the camera pose has.

Feed it back in with:

```bash
python -m eclipse_sites.cli site --lat .. --lon .. --dem d.tif \
    --canopy-profile examples/candamo-canopy.json
```
