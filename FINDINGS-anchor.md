# Anchor investigation, 2026-08-11

Question: the regression anchor says totality clearance at the San Roman de
Candamo pool is +9.8 degrees. A live run against the PNOA 5 m DTM gives +8.6.
Which is right?

Answer: **+8.6**. The +9.8 is wrong. Below is the evidence and the two defects
that produced it.

## 1. What blocks the sun

A real massif, crest 753 m, at **43.48893, -6.29577**, 18.7 km out on the sun
bearing of 280.6 deg. It stands at +1.84 deg from the pool. The profile rises
smoothly from 286 m at 17.2 km to 719 m at 18.7 km and falls away after, with
no nodata anywhere in a 400 m box around the crest. This is not an artifact.

## 2. It is not an eye-elevation artifact

The anchor coordinate lands on a slope where the DEM spans 74 to 124 m within
50 m, so eye elevation was the obvious suspect. It is not the cause.

| eye elev | horizon @280.6 | clearance |
|---|---|---|
| 96.6 m | +1.84 | +8.62 |
| 117.6 m | +1.78 | +8.68 |
| 140 m | +1.71 | +8.75 |
| 500 m | +0.60 | +9.86 |

You would have to stand at 500 m to recover +9.8. Across every plausible eye
elevation the answer moves 11 cm of angle, because the blocker is 18.7 km away.

Fitting eye elevation to the old EU-DEM fan on distant sectors only also fails
to reconcile: best rms 0.80 deg at an implausible 140 m, and the residual stays
structured (+1.41 at az 260, -0.98 at 240) instead of flattening. Eye height
moves all distant sectors coherently. This does not.

## 3. The photograph confirms the massif

`20260730_201008.jpg`, Galaxy S23 ultrawide, 13 mm equivalent, 108 deg of
horizontal field. Fitted with this repo's own `calibrate` path at eye 1.65 m:
bearing 263.78, pitch +15.99, 1.18 deg rms.

| source | horizon at az 280.6 |
|---|---|
| traced skyline, see caveat | +2.34 |
| PNOA 5 m DTM | **+1.84** |
| old EU-DEM fan | +0.64 |

**Caveat on the traced number.** +2.34 comes from `extract_skyline`, which is
defective in this exact sector, see section 6. It is the top of the near
treeline, not the terrain skyline. Treat the photograph as qualitative
confirmation good to about a degree, not as a precision route. Do not quote a
clearance from it. The terrain answer of +8.62 does not depend on any of this.

What the photograph does establish: crop the overlay to az 269.5-291.8 and
there is layered, pale, blue-grey distant terrain behind dark sharp near trees,
with the massif crossing the near ridge around az 279-282. That is atmospheric
perspective. At 3.5 km on an overcast evening you do not get it. At 18.7 km you
do. An open horizon at +0.64 deg with the next ground at 3.5 km cannot look
like this.

Note the pale layer at az 273-278 is *not* the massif. PNOA's horizon there
forms at ~3.2 km. The massif only tops the near ridge from about az 279.

## 4. Defect one: the fan never looked far enough

The stored fan in `tests/test_eclipse_sites.py:227` has 14 km as its furthest
entry. The massif is at 18.7 km. It could not have been seen.

Note the `rng` column of that fan is the range at which the horizon *formed*,
not a search cap. Capping the search at those values does not reproduce the fan
(az 320 capped at 80 m gives 8.02 where the fan says 12.86), so the fan did
search past them. It just did not search past ~14 km.

There is a second, smaller error in the near field. Within 3.5 km at az 280 the
fan says 0.64 and PNOA says 1.45. That 0.8 deg is the direction you expect from
25 m posts smoothing a ridge crest at 3.3 km, and the 5 m LiDAR is better data.
So the gap splits roughly 0.8 deg of near-field smoothing plus 0.4 deg of a
summit nobody looked at.

## 5. Defect two: the two routes were never independent

`fit_camera(skyline_y, intr, dem_az, dem_alt, dem_range, ...)` at
`calibrate.py:219` takes the DEM horizon fan as its reference and solves camera
bearing and pitch against it. The photogrammetric route does not measure
clearance independently. It aligns a photograph to whatever fan you hand it.

Both routes were handed the same EU-DEM fan. So +9.82 and +9.74 agreeing to
0.08 deg is shared-input agreement, not corroboration. That is why the error
survived, and why CLAUDE.md's claim of two independent routes does not hold.

## 6. NEW BUG: extract_skyline loses hazy distant terrain

Found while checking the photograph. `extract_skyline` classifies pixels as sky
by luminance and saturation (`sky_lum=115.0, sky_sat=0.32`). A distant ridge
under atmospheric haze against an overcast sky is pale and desaturated, so it
reads as sky. The trace then follows the dark near treeline *below* the ridge.

Measured on the crest of the pale massif at five points:

| az | pale crest | PNOA bare earth | diff |
|---|---|---|---|
| 273.07 | 4.42 | 2.93 | +1.49 |
| 274.88 | 4.40 | 2.47 | +1.93 |
| 276.66 | 4.19 | 2.41 | +1.78 |
| 278.42 | 3.98 | 2.13 | +1.85 |
| 280.14 | 3.71 | 1.90 | +1.81 |

A consistent ~1.8 deg offset across the sectors the fit is supposed to anchor
on.

**The 1.8 deg is unresolved, and is not a measured bias figure.** It cannot be
pitch error of that size: the median (trace - bare earth) over 1262 distant
columns is +0.26 deg, and a pitch 1.8 deg high would put the median trace about
1.5 deg *below* bare earth, which the non-negative-canopy argument forbids. That
caps any pitch bias near 0.3 to 0.5 deg. The rest is some mix of crest points
eyeballed off a 3x upscaled hazy JPEG (worth +/- 0.4 to 0.5 deg on its own) and
the bearing instability of 13.7 deg seen across the eye scan. A cloud bank
sitting on the ridgeline is also not fully excluded from a single overcast
frame.

What survives on the visual evidence alone is the defect itself: the trace
demonstrably runs below a feature that reads as terrain. Its angular cost is
not yet measured.

This is the same failure class as the rest of the invariant list: it makes the
horizon look *lower* than it is, so it is optimistic, plausible and silent.

It also supplies a mechanism for how the photogrammetric route reached +9.82
even though the massif is plainly visible in the frame: the extractor discarded
it, and the too-low trace matched the too-low EU-DEM fan. Inference, not proof,
since the original run cannot be reproduced.

## Status

Applied since this was written, as a separate piece of work on coverage:
`coverage_point`, `Coverage`, `edge_reach` and `edge_reach_point` in
`horizon.py`, a `coverage` CLI command, 7 tests, invariant 13, and the CLAUDE.md
command list. None of that changes any horizon number.

Still **not** applied, awaiting a decision:

- the CLAUDE.md regression anchor, which still says +9.8 and still claims two
  independent routes. The repo currently contradicts itself: that section says
  "if a change moves it, the change is wrong" while this file says +9.8 is dead.
- the docstring of `test_render_view_candamo_matches_independent_analysis`,
  which repeats the independence claim.
- the README's Tests section, which repeats it again.
