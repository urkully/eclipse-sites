# eclipse-sites

Find places with a clear view of a solar eclipse, using real terrain rather than
guesswork. Built for Asturias and the 12 August 2026 totality, but nothing in it
is specific to either.

Two questions, both answered from a coordinate:

- **Is it viable?** Contact times, Sun position, horizon height in the Sun's
  direction, and how much clear sky sits under it.
- **What will it look like?** A rendered view of the eclipse over the skyline
  as seen from that exact point. No photograph needed: the skyline comes from
  the DEM, so it works for places nobody has ever stood with a camera.

Plus **search**, which sweeps a whole region and ranks the best spots, and
**calibrate**, which fits a photograph to the terrain and measures whatever the
survey cannot see.

## Install

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

No ephemeris download. `pyephem` carries its own theory, so `site` works
offline immediately.

## Quick check, no DEM needed

```bash
python -m eclipse_sites.cli site --lat 43.4524508 --lon -6.0705374 --elev 117.6
```

```
  eclipse type      total
  partial begins    19:31:02
  TOTALITY begins   20:26:50
  maximum           20:27:48
  TOTALITY ends     20:28:45
  duration          115.4 s
  partial ends      21:20:59
  Sun at maximum    altitude 10.46deg  azimuth 280.60deg
```

## Getting terrain

You need a DEM in a **projected CRS in metres**, with **square pixels**.

| Source | Resolution | Notes |
|---|---|---|
| IGN Spain MDT05 | 5 m | best for confirming a shortlist |
| IGN Spain MDT25 | 25 m | right for the regional sweep |
| PNOA LiDAR | ~1 m | surface model, so it includes trees |
| Copernicus EU-DEM | 25 m | whole of Europe, easiest to start with |

IGN products come as EPSG:25829/25830 tiles. Merge and reproject:

```bash
gdalbuildvrt asturias.vrt *.asc
gdalwarp -t_srs EPSG:25830 -tr 25 25 -r bilinear asturias.vrt coarse.tif
```

### Let the tool do the merging

Point `prepare` at whatever came out of the download. Zips, subfolders, mixed
products and mixed UTM zones are all fine.

```bash
python -m eclipse_sites.cli prepare --in ~/Downloads/cnig --out ./terrain
```

It sorts by product, works out the CRS, merges each group and writes one clean
GeoTIFF per product into the directory you name, plus a `manifest.json`. Output
is `terrain/dtm05.tif`, `terrain/dtm25.tif`, `terrain/ndsm_vegetation.tif` and
so on, named after the product rather than the sheet numbers.

IGN ASCII grids often carry no projection. The UTM zone is read from the file
name (HU29, HU30) where possible; if that fails, pass `--src-crs EPSG:25829`.
Candamo is in zone 29, most of the rest of Asturias is zone 30.

Two things it handles that hand-rolled GDAL commands usually get wrong.
Normalised products such as `ndsm_vegetation` are resampled with MAX rather than
bilinear, because averaging a treeline loses it and a horizon needs the tallest
thing in the cell. And gaps between tiles are written as declared nodata, never
as zero: a gap read as 0 m ground looks like a cliff and fakes a clear view.

If you would rather drive GDAL yourself, `gdalinfo`, `gdalwarp` and
`gdalbuildvrt` come from the GDAL command line tools, which are a separate
install (`brew install gdal`, `apt install gdal-bin`). Rasterio bundles the
libraries but not the binaries. Its own `rio info`, `rio warp` and `rio merge`
are already on your path.

### The buffer rule, which bites people

A cell's horizon is only knowable if its ray stays inside the raster. With
`--range 40000`, every cell within 40 km of the edge has an **unknown** horizon
and is returned as nodata, not as a wide-open view. So your DEM must extend
`--range` beyond the area you care about. Either download a generous margin or
lower `--range`.

Do not reason that a short range is enough because nothing nearby is tall
enough to hide the Sun. That is true and it is the wrong test. Fully blocking a
10.5° Sun at 15 km does take a ~2900 m mountain, but *clearance* is the number
this tool ranks sites on, and clearance only needs the horizon to move. A 753 m
hill at 18.7 km stands at +1.8°, hides nothing, and still costs 1.2° of
clearance. That is precisely how the Candamo anchor came to be 1.2° too
optimistic for months. Set `--range` from the clearance precision you want, not
from what could eclipse the Sun outright.

To find out what radius a site can honestly be judged on before you trust any
horizon from it:

```bash
python -m eclipse_sites.cli coverage \
    --lat 43.4524508 --lon -6.0705374 --dem terrain/dtm05.tif --range 40000
```

It reports the free bounds test first, then marches every azimuth for the real
answer, which also catches nodata inside the bounds. `contiguous_m` stops at the
first gap rather than the last data, because ground seen beyond a hole does not
vouch for the hole.

## What it will look like

```bash
python -m eclipse_sites.cli site \
    --lat 43.4524508 --lon -6.0705374 --elev 117.6 \
    --dem mdt05.tif --range 15000 \
    --view candamo.png
```

Renders the Sun's whole track from first to fourth contact over the terrain
skyline, with quarter-hour marks, contact markers, a corona ring at totality,
and a verdict banner. Ridges are shaded by the distance at which they form the
horizon, near dark and close, far pale and hazy, so it reads as a landscape
rather than a chart.

`--fov` sets the field of view if you want it tighter or wider than the
default, which frames the whole track plus a margin.

## Measuring obstructions from a photo

DEMs are bare earth and LiDAR flights are years old. A photo knows about the
trees, the new barn and the hedge that got away. `calibrate` extracts them.

```bash
python -m eclipse_sites.cli calibrate \
    --image pool.jpg --lat 43.4524508 --lon -6.0705374 \
    --dem mdt05.tif --range 15000 \
    --mask 0,0,560,1848 \
    --out-csv canopy.csv --out-overlay check.jpg
```

The order matters and is the whole point. Terrain fixes the camera geometry
first; only then is the leftover attributed to obstructions. Ask a photo for
the camera angles and the obstruction height simultaneously and you have one
equation with two unknowns, so bearing error silently absorbs the obstruction
and you get a confident wrong answer.

Anchoring is what keeps them separate. By default the fit uses only azimuths
where the horizon forms beyond `--anchor-min-range` (2 km), because 20 m of
canopy is worth 0.33 deg at 3.5 km and 14 deg at 80 m. Fit over near wooded
ground instead and the canopy is absorbed into the camera angles, destroying the
measurement. There is a test for exactly this.

`--mask x0,y0,x1,y1` removes foreground clutter and is repeatable. Always check
`--out-overlay`: yellow is the traced skyline, cyan is what the DEM predicts
given the fit, and the green bar marks the anchor sectors.

### Read the warnings

Two are worth stopping for. A **rival solution** means the skyline is too
repetitive to pin the bearing; supply `--bearing-prior` with a tighter
`--bearing-window`, or use a photo containing more distinctive distant relief.
A **high residual** means the trace is probably wrong, so look at the overlay
before believing any canopy number.

Sanity-check the output yourself. A 3 deg residual at 400 m is 21 m of canopy,
which is a tree. The same residual attributed to a 3 km horizon is 157 m, which
is not a tree; it means the fit is wrong.

There is also a **trace check**. A photograph cannot see ground below the
ground, so a traced skyline sitting under bare earth is a defect, and the usual
cause is pale distant terrain being read as sky. The trace then drops to nearer
treetops and every angle in that sector is wrong. `extract_skyline` also ends
the sky wherever the colour steps by more than `--sky-tol` from the sky just
above, which on the overcast Candamo frame catches the hazy ridges the other
two thresholds let through. If the check still fires, raise `--sky-lum`, lower
`--sky-sat` or lower `--sky-tol` until it clears. Before `--sky-tol` existed,
tuning the first two on that frame moved the fit from 1.18 to 0.45 deg rms.

### Feeding it back

```bash
python -m eclipse_sites.cli calibrate ... --out-canopy canopy.json
python -m eclipse_sites.cli site --lat .. --lon .. --dem d.tif \
    --canopy-profile canopy.json
```

`--out-canopy` writes the measured cover per azimuth, and `site` applies it
instead of the flat `--canopy` scalar, falling back to that scalar outside the
arc the photograph covered. The clearance line says when a profile contributed.

Only bins that can honestly be converted are written: the horizon must form
within 5 km, since angle-to-metres scales with range and at 18.7 km a fifth of
a degree is 65 m of imaginary tree. Bins failing that, or implying a negative or
absurd height, are dropped, and a dropped bin is **not** rebuilt by
interpolating its neighbours. A poor fit refuses to write a profile at all,
because canopy inherits every error the camera pose has.

This is still a crude model: one height per azimuth, applied along the whole
ray. The real fix is the `z_surface` second raster, which is designed and not
built.

### The version that needs none of this

Stand on the spot at the right minute on a clear evening near the date and
photograph the Sun itself against the skyline. No bearing, no pitch, no fit, no
residual. For one site that beats everything above. Calibration earns its keep
when you want to check twenty candidate sites without visiting each one.

## Regional search

```bash
python -m eclipse_sites.cli search \
    --dem coarse.tif --fine mdt05.tif \
    --range 15000 --top 40 \
    --out-raster clearance.tif --out-csv sites.csv
```

Staged on purpose. The coarse pass sweeps everything, then only the shortlist is
re-checked on the fine DEM with the exact ray march. You never process a 5 m
raster over a whole region.

Output: `clearance.tif` is degrees of sky between the Sun and the skyline at
maximum eclipse. Positive means visible. Open it in QGIS. `sites.csv` is the
ranked shortlist with coordinates.

Eclipse circumstances are computed on a lattice across the raster and
interpolated per cell, so the Sun's azimuth drifts correctly across a large
region instead of being pinned to one town's value.

## Options worth knowing

| Flag | Default | Why |
|---|---|---|
| `--eye` | 1.6 | standing eye height |
| `--range` | 40000 | horizon search distance, see the buffer rule |
| `--canopy` | 0 | metres of tree cover added to ground heights, on both passes |
| `--az-step` | 2.5 | azimuth resolution of the sweep |
| `--min-sep` | 20 | keeps results from clustering in one field |

## Limitations, in the order likely to bite you

**Tree canopy.** DEMs are bare earth. At 80 m a 20 m tree is worth 14°, which is
enough to flip a verdict on its own. `--canopy` is a blunt global offset. If
close obstructions matter, use a LiDAR surface model instead.

Note what "global" means. The offset raises every cell, including the ones
immediately around you, so on flat ground `--canopy 20` puts 20 m of trees in
the next cell and the horizon comes out near vertical. That is the honest
answer for someone standing inside a forest, and it is why the flag is a
screening tool rather than a measurement. To find out what is actually there,
use `calibrate` against a photograph.

**Near-field quantisation in the grid engine.** The sweep snaps rays to whole
pixels, so azimuths within a few cells are approximate. This is why `--fine`
exists: `horizon_point` marches the true ray with bilinear sampling. Trust the
sweep for shortlisting, never for a final answer.

**Vertical datum.** DEMs are geoid-referenced, phone GPS is ellipsoidal, and the
separation in northwest Spain is roughly 50 m. Never mix the two. Take your
observer elevation from the same DEM you are querying.

**Sea horizons.** Ray marching over water returns nothing and the cell is masked.
A coastal site with a genuine sea horizon needs the dip formula,
`0.0293 * sqrt(height_m)` degrees below level. Not implemented.

**Contact times.** Computed from pyephem's truncated theory. C1 and C4 land
within about a minute of published values. Totality *duration* is the weak
number: it depends on the small difference between two apparent radii, and this
code gives 115 s where local sources say 109 s. Use it for planning, not for
printing on a poster.

**Refraction.** A fixed coefficient, k = 0.13. Fine at 10° altitude. Near the
horizon the real value is variable and worth about half a degree.

## Tests

```bash
python -m pytest tests/ -v
```

63 tests, no network, no DEM. Terrain cases are synthetic with analytic answers:
a wall at a known distance, azimuth convention, grid against point engine,
horizon falling as you retreat from a cone. The ephemeris case is the real 2026
eclipse checked against published local timings. There is also a brute-force
check of the array shift over 255 offsets, which caught a real wrapping bug where
oversized shifts silently returned wrong data.

Calibration is proved by round trip: invent a camera, synthesise the skyline it
would see, and check the fit recovers it. Clean data comes back at 283.00 deg
bearing and 14.68 deg pitch against a true 283.0 and 14.7. Injecting 8 deg of
canopy on near ground and anchoring correctly still recovers the camera;
anchoring on the trees instead throws the bearing 13 deg out, which is the
documented failure mode and is asserted as such. Canopy height comes back at
11.2 m where 11.2 m was injected.

The renderer is checked against a stored horizon fan for the Candamo pool, which
pins the render maths but is **not** independent evidence about that site. It
was once described here as corroborated by a photogrammetric route agreeing to
0.08 degrees. That was wrong twice over: `fit_camera` solves a camera against
whichever DEM fan you hand it, so it cannot disagree with its own input, and the
fan itself stopped looking at 14 km and missed the 753 m massif at 18.7 km that
actually forms the horizon. The live answer for that site is +8.6 degrees, not
+9.8. See `FINDINGS-anchor.md`.

## Layout

```
eclipse_sites/
  ephemeris.py   Sun and Moon, contacts, obscuration
  horizon.py     horizon_point (exact), horizon_grid (fast),
                 coverage_point and edge_reach (what radius is knowable)
  render.py      synthetic view of the eclipse over the skyline
  calibrate.py   fit a photo to terrain, measure canopy from the residual
  prepare.py     sort, merge and reproject a folder of downloaded tiles
  search.py      DEM IO, regional sweep, ranking, confirmation
  cli.py         command line
tests/
```

`horizon.py` and `ephemeris.py` have no geospatial dependencies and take plain
numpy arrays, so they are easy to reuse or test in isolation.
