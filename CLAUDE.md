# eclipse-sites

# Claude instructions

Behavioral guidelines to reduce common LLM coding mistakes. Merge with project-specific instructions as needed.

**Tradeoff:** These guidelines bias toward caution over speed. For trivial tasks, use judgment.

## 1. Think Before Coding

**Don't assume. Don't hide confusion. Surface tradeoffs.**

Before implementing:
- State your assumptions explicitly. If uncertain, ask.
- If multiple interpretations exist, present them - don't pick silently.
- If a simpler approach exists, say so. Push back when warranted.
- If something is unclear, stop. Name what's confusing. Ask.

## 2. Simplicity First

**Minimum code that solves the problem. Nothing speculative.**

- No features beyond what was asked.
- No abstractions for single-use code.
- No "flexibility" or "configurability" that wasn't requested.
- No error handling for impossible scenarios.
- If you write 200 lines and it could be 50, rewrite it.

Ask yourself: "Would a senior engineer say this is overcomplicated?" If yes, simplify.

## 3. Surgical Changes

**Touch only what you must. Clean up only your own mess.**

When editing existing code:
- Don't "improve" adjacent code, comments, or formatting.
- Don't refactor things that aren't broken.
- Match existing style, even if you'd do it differently.
- If you notice unrelated dead code, mention it - don't delete it.

When your changes create orphans:
- Remove imports/variables/functions that YOUR changes made unused.
- Don't remove pre-existing dead code unless asked.

The test: Every changed line should trace directly to the user's request.

## 4. Goal-Driven Execution

**Define success criteria. Loop until verified.**

Transform tasks into verifiable goals:
- "Add validation" → "Write tests for invalid inputs, then make them pass"
- "Fix the bug" → "Write a test that reproduces it, then make it pass"
- "Refactor X" → "Ensure tests pass before and after"

For multi-step tasks, state a brief plan:
```
1. [Step] → verify: [check]
2. [Step] → verify: [check]
3. [Step] → verify: [check]
```

Strong success criteria let you loop independently. Weak criteria ("make it work") require constant clarification.

---

**These guidelines are working if:** fewer unnecessary changes in diffs, fewer rewrites due to overcomplication, and clarifying questions come before implementation rather than after mistakes.

## Mission

Given a latitude and longitude, answer two questions:

1. **Is the eclipse viable from here?** Contact times, Sun position, terrain
   horizon in the Sun's direction, clearance.
2. **What will it look like?** A rendered view from that exact point.

Plus regional search for better sites, photo calibration to measure what the
terrain survey cannot see, and data preparation.

## Commands

```bash
python -m pytest tests/ -q                      # 63 tests, no network, no DEM
python -m eclipse_sites.cli site --lat .. --lon ..          # works with no DEM
python -m eclipse_sites.cli coverage --lat .. --lon .. --dem d.tif --range ..
python -m eclipse_sites.cli prepare --in <dump> --out <dir>
python -m eclipse_sites.cli site --lat .. --lon .. --dem d.tif --view v.png
python -m eclipse_sites.cli search --dem coarse.tif --fine fine.tif
python -m eclipse_sites.cli calibrate --image p.jpg --lat .. --lon .. --dem d.tif \
    --out-canopy canopy.json   # then: site --canopy-profile canopy.json
```

## Architecture

| Module | Role | Dependencies |
|---|---|---|
| `ephemeris.py` | Sun and Moon, contacts, obscuration | pyephem only |
| `horizon.py` | `horizon_point` (exact ray march), `horizon_grid` (fast), `coverage_point` and `edge_reach` (what radius is even knowable) | numpy only |
| `render.py` | synthetic view over the skyline | + Pillow |
| `calibrate.py` | fit a photo to terrain, canopy from residual | + Pillow |
| `search.py` | DEM IO, regional sweep, ranking | + rasterio, pyproj |
| `prepare.py` | sort, merge, reproject downloaded tiles | + rasterio |

`ephemeris.py` and `horizon.py` deliberately have no geospatial dependencies and
take plain numpy arrays. Keep it that way; it is what makes them testable.

## Invariants. Every one of these is a bug that was found and fixed

Do not regress these. Each has a test.

1. **Unknown horizon must never read as open.** Cells whose ray leaves the raster
   return NaN, not the `-90` sentinel. The first regional sweep ranked its best
   sites at +96 degrees clearance because edge cells kept the sentinel.
1b. **A cell with no elevation of its own is also unknown.** Where the observer
   sits on nodata, `eye` is NaN, every angle along the ray is NaN, and
   `np.fmax` ignores NaN, so the sentinel survives. Checking only the far end
   of the ray does not catch it. This ranked a nodata strip as the eight best
   sites in Asturias at +100 degrees.
1c. **`horizon_point` reports `reach_m`.** The grid engine masks what it could
   not see; the point engine is what produces every final answer, so it says
   how far it actually looked. A ray that ran out of DEM at 3 km has not shown
   that the next 37 km are clear, it has only failed to look.
2. **Tile gaps must be declared nodata.** A gap filled with 0 reads as ground at
   sea level, looks like a cliff, and fakes a clear view. Optimistic, plausible,
   silent.
3. **The camera fit must be penalised for low coverage.** Without it the solver
   slides the frame until only a sliver overlaps the anchor sectors and fits
   that perfectly. A real run returned 0.05 degrees rms from 24 columns and a
   bearing 56 degrees wrong.
4. **Calibration anchors only on distant terrain.** 20 m of canopy is worth 0.33
   degrees at 3.5 km and 14 degrees at 80 m. Fit over near wooded ground and the
   canopy is absorbed into the camera angles, destroying the measurement. Costs
   13 degrees of bearing when violated.
5. **Observer elevation comes from bare terrain, never a surface model.** A DSM
   will place the eye on top of a tree.
6. **Normalised surface products resample with MAX, not bilinear.** Averaging a
   treeline loses it. Note MAX is warp-only and `merge()` rejects it, so warp
   each tile to the target grid first, then merge with nearest.
7. **Array shifts must clamp bounds.** Python treats a negative slice stop as an
   index from the end, so oversized shifts silently returned wrong terrain.
8. **Bisection brackets must be ordered.** A reversed bracket made C3 and C4
   silently return the end of the search window instead of a contact.
9. **The bearing search sweeps densely before refining.** Coarse-to-fine alone
   latches onto false minima on a repetitive ridgeline.
10. **`canopy_m` is metres, everywhere.** It goes to the horizon engine to be
   added to sampled ground, never onto an angle. `sweep` used to add it to the
   horizon in degrees, so `--canopy 20` meant 20 m on the fine pass and a flat
   20 degrees on the coarse one, and the two clearance columns in the same
   output table disagreed by tens of degrees.
11. **The reported residual is a true rms in degrees.** The coverage penalty
   that stops the sliver overfit is a search cost, not an angle. Reporting the
   penalised number as degrees inflates it exactly when coverage is poor,
   which is when the operator most needs it to mean what it says.
12. **Azimuth is a direction, not a number.** Binning, the rival-bearing test,
   and the render frame all compare on the circle. A full-circle horizon fan
   is a ring: 0 to 359.5 has a hole in it as far as `np.interp` is concerned.
13. **A horizon is only as good as the radius that produced it.** Coverage
   stops at the first gap, not the last data: ground seen beyond a hole does
   not vouch for the hole. `coverage_point` reports `contiguous_m` for this,
   separately from the looser `reach_m` both engines share. `edge_reach` is
   the free bounds test and is necessary, never sufficient, because it cannot
   see nodata inside the bounds. Cross-checking a horizon against a second
   route that consumes the same fan proves only that the arithmetic agrees.

## Testing conventions

Synthetic terrain with analytic answers: a wall at known distance, azimuth
convention, grid engine against point engine, horizon falling as you retreat
from a cone. Ephemeris is checked against the real 2026-08-12 eclipse with
published local timings. Calibration is proved by round trip: invent a camera,
synthesise the skyline it would see, recover it.

**Test fixtures in this project have been wrong twice, and both times the code
was right.** A 500 m raster asked to march 2000 m, and a nodata test using
perfectly adjacent tiles with no gaps. When a test fails, check the fixture
before changing the code.

## Regression anchor

The pool at San Roman de Candamo, 43.4524508, -6.0705374, eye 1.65 m above DEM
ground. Totality clearance is **+8.6 degrees** on the PNOA 5 m LiDAR DTM. The
horizon at the Sun's bearing is +1.8 degrees and is formed by a 753 m massif at
43.48893, -6.29577, 18.7 km out. If a change moves this, the change is wrong.

The number is insensitive to eye height, which is why the convention above is
safe: the blocker is 18.7 km away, so 96 m to 140 m of eye elevation moves the
answer by 0.1 degrees. You would have to stand at 500 m to see +9.8.

An earlier anchor of +9.8 was wrong and is kept here as a warning. It came from
an EU-DEM fan whose furthest entry was 14 km, so it never saw the massif, and
which also read 0.8 degrees low in the near field where 25 m posts smooth a
crest at 3.3 km. It was believed because a photogrammetric route "agreed" to
0.08 degrees. It could not have disagreed: `fit_camera` takes the DEM fan as
its reference and solves the camera against it, so it was handed the same wrong
fan and gave it back. Two routes sharing an input are one route. See
`FINDINGS-anchor.md`, and invariant 13.

Contacts there: C1 19:31, C2 20:26:50, C4 21:20:59 CEST.

## Known limitations, in order of how likely they are to matter

- **Canopy.** DEMs are bare earth. `--canopy` is a crude global offset. Second
  coverage LiDAR vegetation data is 2015-2021, so treat it as a floor.
  A second-raster path (`z_surface` separate from `z_ground`) is designed but
  not built.
- **Near-field quantisation** in `horizon_grid`; rays snap to whole pixels. Use
  `horizon_point` for final answers, which is why `--fine` exists.
- **Sea horizons** are unhandled. Ray marching over water returns nothing and
  the cell is masked. Needs the dip formula, `0.0293 * sqrt(height_m)` degrees.
- **Totality duration** comes out ~6% long (115 s where local sources say 109 s)
  because it depends on a small difference between two apparent radii. Fine for
  planning, not for publishing.
- **Vertical datum.** DEMs are geoid-referenced, GPS is ellipsoidal, roughly
  50 m apart in northwest Spain. Never mix sources.

## Working style

Plain, direct sentences of varied length. No em dashes. Flag uncertainty
explicitly rather than sounding confident. Ask rather than guess. Prefer showing
a draft to react to over open questions.
