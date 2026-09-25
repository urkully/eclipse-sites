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

**Superseded by the refit in section 6.** The camera above and the +2.34 below
were fitted against the defective trace. The current numbers are bearing
265.12, pitch +13.83, roll -1.63, 0.23 deg rms, and +1.70 at az 280.6.

**Caveat on the traced number.** +2.34 comes from `extract_skyline`, which was
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

## 6. FIXED: extract_skyline lost hazy distant terrain

Fixed in PR #1, merged as 3fd6fca. The diagnosis below is kept as
written; the fix and the refit it forced are recorded after it.

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

### The fix

A pixel that passes the absolute luminance and saturation test must now also
stay within `--sky-tol` of the running mean of the sky in the `run` pixels
above it. The reference follows the sky down the column, so a gradient passes
and a ridge does not. Measured on this frame:

| quantity | levels |
|---|---|
| hazy ridge step down from the sky | ~90 |
| cloud texture, 99.999th percentile | ~27 |
| cloud texture, worst | 50 |
| **default `--sky-tol`** | **60** |

That is only about 10 levels of margin over the worst cloud step, on one
overcast frame. Clear skies, sunsets and strong gradients are untested. On
this frame the new trace sits 1.2 to 2.8 deg above the old one depending on
sector.

### The refit, 2026-09-25

The camera in section 3 was fitted against the defective trace, so it had to be
redone before anything from the photograph could be quoted. Same frame, same
eye 1.65 m, `--hfov 108`, `--mask 0,1100,830,1850`. The mask was checked on the
overlay: the trace is suppressed across both parasols and resumes at x=830, and
it stays above the flag tops.

| free parameters | bearing | pitch | roll | focal | rms | below bare earth | az 280.6 |
|---|---|---|---|---|---|---|---|
| old fit, defective trace | 263.78 | +15.99 | - | - | 1.18 | - | +2.34 |
| bearing, pitch | 266.49 | +14.93 | - | - | 0.49 | 36% | +2.70 |
| + roll | 266.79 | +14.80 | -1.51 | - | 0.34 | 20% | +2.22 |
| + focal | 265.30 | +14.13 | - | x1.078 | 0.42 | 30% | +2.60 |
| **+ roll + focal** | **265.12** | **+13.83** | **-1.63** | **x1.098** | **0.23** | **5%** | **+1.70** |

PNOA bare earth at az 280.6 is +1.83, so the best fit lands 0.14 deg *below*
it, against +0.87 with bearing and pitch alone. The worry going in was that the
corrected trace sat 1.3 to 1.6 deg above PNOA at the five crest points above
while still wearing the old camera. Roll and focal together close that.

Both extra parameters are physically plausible rather than slack being
absorbed. Roll -1.63 deg is a handheld phone. Focal x1.098 turns the assumed
108 deg field into 102.8 deg, so the `--hfov 108` guess for a Galaxy S23
ultrawide with no EXIF focal length was about 5 deg too wide. Fitting either
one alone leaves the trace-check warning standing; fitting both clears it, and
the share of distant columns below bare earth falls from 36% to 5%.

The recommended command for this frame is therefore:

```bash
python -m eclipse_sites.cli calibrate \
    --image 20260730_201008.jpg \
    --lat 43.4524508 --lon -6.0705374 --dem <pnoa_dtm.tif> \
    --eye 1.65 --hfov 108 --mask 0,1100,830,1850 \
    --fit-roll --fit-focal
```

**This is still not an independent route**, and the caveat in section 3 stands
for the same reason as ever: `fit_camera` solves the camera against the PNOA
fan, so agreement at az 280.6 is not confirmation. What the refit buys is that
the photograph no longer *contradicts* PNOA, and that the residual is small
enough for the canopy measurement to mean something. The terrain answer of
+8.62 does not depend on any of it.

Median gap is -0.05 deg, still marginally negative where canopy should make it
positive. Small enough to be trace noise on a hazy frame, not chased further.

### How much margin `--sky-tol` actually has, on this frame

Swept with the best camera free (roll and focal), everything else fixed:

| `--sky-tol` | bearing | pitch | roll | rms | below bare earth | az 280.6 |
|---|---|---|---|---|---|---|
| 20 | 289.75 | +11.61 | +5.65 | 1.72 | 63% | +0.73 |
| 30 | 265.51 | +13.97 | -1.60 | 0.50 | 11% | +2.12 |
| 40 | 265.23 | +13.95 | -1.69 | 0.24 | 7% | +1.75 |
| 50 | 265.13 | +13.83 | -1.67 | 0.23 | 4% | +1.75 |
| **60** | **265.12** | **+13.83** | **-1.63** | **0.23** | **5%** | **+1.70** |
| 70 | 264.97 | +13.79 | -1.68 | 0.28 | 5% | +1.58 |
| 85 | 264.59 | +13.75 | -1.82 | 0.67 | 11% | +1.69 |
| 100 | 263.91 | +14.29 | -1.91 | 1.22 | 34% | +2.86 |
| 130 | 264.60 | +15.01 | -0.73 | 1.18 | 34% | +2.59 |
| 200 | 264.60 | +15.01 | -0.73 | 1.18 | 34% | +2.59 |

There is a stable plateau from about 40 to 70 where the fit barely moves. The
default of 60 sits near its top, so on this frame the margin is asymmetric:
roughly 20 levels of room on the tight side and 10 on the loose side. **50
would be a slightly better centred default**, but on one frame that is not
enough to justify moving it.

Both failure modes are visible and they fail differently. Too loose degrades
gracefully into the original bug: 130 and 200 give identical results, meaning
the test has saturated and stopped doing anything, and the fit lands at bearing
264.60 pitch +15.01, close to the old defective 263.78 / +15.99. Too tight
fails hard and early: at 20 the bearing is 289.75, nearly 25 deg wrong, with
roll swinging to +5.65 and 63% of distant columns below bare earth. A too-tight
tolerance is therefore the more dangerous setting, because it does not creep,
it jumps.

### Bug found during the refit

`--fit-focal` without `--fit-roll` crashed with an `IndexError`. The solver
packs free parameters in order, so the focal scale is at index 2 when roll is
fixed and index 3 when it is not, and `fit_camera` read a fixed `p[3]` in two
places. Had the lengths ever lined up it would have silently read the roll as a
focal scale instead of crashing. Fixed, with a round-trip test that synthesises
through a longer lens than the fit is told about and recovers the scale.

## Status

Applied since this was written, as a separate piece of work on coverage:
`coverage_point`, `Coverage`, `edge_reach` and `edge_reach_point` in
`horizon.py`, a `coverage` CLI command, 7 tests, invariant 13, and the CLAUDE.md
command list. None of that changes any horizon number.

All three items that were once outstanding are now applied. The CLAUDE.md
regression anchor reads +8.6 and carries +9.8 as a warning; the docstring of
`test_render_view_candamo_matches_independent_analysis` and the README's Tests
section both state that the route was never independent.

Section 6 is fixed and the camera refitted. Nothing here is awaiting a
decision.

The one open item is `--sky-tol` on frames other than this one. The sweep above
maps the safe band on the 2026-07-30 frame, but it is still one overcast
evening. Clear skies, sunsets and strong gradients are untried. Testing them
needs a second frame with known coordinates inside a DEM, which is not
something more analysis of this frame can supply.
