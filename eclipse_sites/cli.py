"""Command line interface.

  eclipse-sites site   --lat .. --lon .. [--dem tif]   one location
  eclipse-sites search --dem coarse.tif [--fine fine.tif]  whole region
"""
from __future__ import annotations

import argparse
import datetime as dt
import os
import json
import sys

DEFAULT_DATE = "2026-08-12"
DEFAULT_TZ = 2.0


def _bar(label):
    def f(i, n):
        if not n:
            return
        pct = int(100 * i / n)
        sys.stderr.write(f"\r{label}: {pct:3d}%  ({i}/{n})")
        sys.stderr.flush()
        if i == n:
            sys.stderr.write("\n")
    return f


def _local(t, tz):
    return None if t is None else (t + dt.timedelta(hours=tz)).strftime("%H:%M:%S")


def cmd_site(a):
    from .ephemeris import circumstances
    d = dt.date.fromisoformat(a.date)
    c = circumstances(a.lat, a.lon, a.elev, d)
    print(f"\n{a.lat:.6f}, {a.lon:.6f}   {a.date}   UTC{a.tz:+g}")
    print(f"  eclipse type      {c.kind}")
    if c.kind == "none":
        print("  no eclipse here on this date")
        return
    print(f"  partial begins    {_local(c.c1, a.tz)}")
    if c.c2:
        print(f"  TOTALITY begins   {_local(c.c2, a.tz)}")
        print(f"  maximum           {_local(c.cmax, a.tz)}")
        print(f"  TOTALITY ends     {_local(c.c3, a.tz)}")
        print(f"  duration          {c.duration_s:.1f} s")
    else:
        print(f"  maximum           {_local(c.cmax, a.tz)}"
              f"   obscuration {c.max_obscuration*100:.1f}%")
    print(f"  partial ends      {_local(c.c4, a.tz)}")
    print(f"  Sun at maximum    altitude {c.sun_alt_at_max:.2f}deg  "
          f"azimuth {c.sun_az_at_max:.2f}deg")

    def write_json(extra=None):
        if not a.json:
            return
        out = c.as_dict()
        out.update(extra or {})
        with open(a.json, "w") as fh:
            json.dump(out, fh, indent=2)
        print(f"  wrote {a.json}")

    if not a.dem:
        print("\n  (no --dem given, so no horizon check)")
        write_json()
        return

    from pyproj import Transformer
    from .horizon import horizon_point
    from .search import read_dem
    z, tr, crs, cell = read_dem(a.dem)
    to_dem = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    x, y = to_dem.transform(a.lon, a.lat)
    col, row = (~tr) * (x, y)
    ny, nx = z.shape
    if not (0 <= row < ny and 0 <= col < nx):
        to_wgs = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
        c0 = to_wgs.transform(*(tr * (0, ny)))
        c1 = to_wgs.transform(*(tr * (nx, 0)))
        sys.exit(f"\n  that point is outside the DEM.\n"
                 f"  DEM covers lat {c0[1]:.4f} to {c1[1]:.4f}, "
                 f"lon {c0[0]:.4f} to {c1[0]:.4f}")
    floor_at = _skyline_floor(a)
    hz, hr, reach = horizon_point(
        z, cell, row - 0.5, col - 0.5, [c.sun_az_at_max],
        max_range_m=a.range, eye_h=a.eye, canopy_m=a.canopy)[c.sun_az_at_max]
    photo_hz, hz_bare = floor_at(c.sun_az_at_max), hz
    raised = photo_hz > hz
    hz = max(hz, photo_hz)
    clear = c.sun_alt_at_max - hz
    truncated = reach < a.range - cell
    print(f"\n  horizon at az {c.sun_az_at_max:.1f}deg   {hz:+.2f}deg "
          f"(formed at {hr:.0f} m)")
    print(f"  CLEARANCE         {clear:+.2f}deg")
    if getattr(a, "canopy_profile", None):
        # invariant 11: a number that quietly includes photo-measured cover is
        # not the same number, so it has to say so where it is printed
        said = (f"raised it to this from {hz_bare:+.2f}deg" if raised
                else "says nothing at this bearing")
        print(f"  (photograph {said}, from {a.canopy_profile})")
    print("  verdict           " + (
        "VISIBLE, comfortable" if clear > 5 else
        "VISIBLE, tight - confirm on the ground" if clear > 0 else
        "BLOCKED by terrain"))
    if truncated:
        print(f"\n  WARNING: the ray ran out of DEM at {reach:.0f} m of the "
              f"{a.range:.0f} m asked for.\n  Ground beyond that is unknown, "
              f"not clear. A distant ridge outside the\n  raster would not "
              f"appear above. Extend the DEM or lower --range to say\n  what "
              f"you actually checked.")

    if a.view:
        import numpy as np
        from .render import render_view
        azs = list(np.arange(c.sun_az_at_max - 45, c.sun_az_at_max + 45.1, 1.0))
        # horizon_point loops per azimuth anyway, so one scalar call each costs
        # the same as a vector API would and leaves the engine untouched
        h = horizon_point(z, cell, row - 0.5, col - 0.5, azs,
                          max_range_m=a.range, eye_h=a.eye, canopy_m=a.canopy)
        h = {q: v._replace(alt_deg=max(v.alt_deg, floor_at(q)))
             for q, v in h.items()}
        render_view(a.lat, a.lon, d, azs, [h[x].alt_deg for x in azs],
                    [h[x].range_m for x in azs], a.view, elev_m=a.elev,
                    tz=a.tz, fov=a.fov)
        print(f"  wrote {a.view}")

    write_json({"horizon_deg": hz, "clearance_deg": clear,
                "horizon_range_m": hr, "reach_m": reach,
                "ray_truncated": truncated})


def cmd_search(a):
    from .search import sweep, top_sites, confirm, write_geotiff, write_csv
    d = dt.date.fromisoformat(a.date)
    clearance, z, tr, crs, cell, circ = sweep(
        a.dem, d, eye_h=a.eye, max_range_m=a.range, az_step=a.az_step,
        canopy_m=a.canopy, progress=_bar("sweeping"))
    print(f"\ncoarse pass done on {z.shape[1]}x{z.shape[0]} cells at {cell:g} m")
    if circ and circ.c2:
        print(f"sample totality {circ.duration_s:.0f} s, "
              f"Sun alt {circ.sun_alt_at_max:.1f}deg az {circ.sun_az_at_max:.1f}deg")

    if a.out_raster:
        write_geotiff(a.out_raster, clearance, tr, crs)
        print(f"wrote {a.out_raster}")

    sites = top_sites(clearance, z, tr, crs, n=a.top, min_sep_cells=a.min_sep)
    if a.fine:
        print(f"confirming top {len(sites)} on the fine DEM ...")
        sites, dropped = confirm(a.fine, sites, d, eye_h=a.eye,
                                 canopy_m=a.canopy, max_range_m=a.range)
        if dropped:
            print(f"  {len(dropped)} of them are not covered by the fine DEM "
                  f"and could not be confirmed:")
            for s in dropped[:5]:
                print(f"    coarse rank {s.get('rank','?'):>3}  "
                      f"{s['lat']:.5f}, {s['lon']:.5f}")
            if len(dropped) > 5:
                print(f"    ... and {len(dropped) - 5} more")
        cut = [s for s in sites if s.get("ray_truncated")]
        if cut:
            print(f"  WARNING: {len(cut)} of the confirmed sites had the ray "
                  f"leave the fine DEM\n  before --range. Their horizons are "
                  f"only known as far as reach_m.")
    if a.out_csv:
        write_csv(a.out_csv, sites)
        print(f"wrote {a.out_csv}")

    key = "fine_clearance_deg" if a.fine else "clearance_deg"
    print(f"\n{'rank':>4} {'lat':>10} {'lon':>11} {'elev':>7} {'clear':>8}")
    for s in sites[:20]:
        print(f"{s['rank']:>4} {s['lat']:>10.5f} {s['lon']:>11.5f} "
              f"{s['elev_m']:>7.0f} {s.get(key, 0):>+8.2f}")


def cmd_prepare(a):
    from .prepare import prepare, scan
    found = scan(a.indir, a.src_crs)
    if not found:
        sys.exit(f"  no raster files under {a.indir}")
    print(f"  found {len(found)} files")
    seen = {}
    for t in found:
        seen[t.product] = seen.get(t.product, 0) + 1
    for k, v in sorted(seen.items()):
        print(f"    {k:20s} {v:4d} tiles")

    rep = prepare(a.indir, a.outdir, src_crs=a.src_crs,
                  target_crs=a.target_crs, only=a.only, res=a.res,
                  max_cells=a.max_cells, progress=_bar("preparing"))
    if not rep.written:
        print("\n  nothing written.")
    else:
        print(f"\n  target CRS {rep.target_crs}")
        for w in rep.written:
            print(f"    {w['path']}  {w['size'][0]}x{w['size'][1]} @ "
                  f"{w['res_m']:g} m  ({w['tiles']} tiles, {w['resampling']})")
    for s_ in rep.skipped[:8]:
        print(f"  skipped {os.path.basename(s_['path'])}: {s_['why']}")
    if len(rep.skipped) > 8:
        print(f"  ... and {len(rep.skipped) - 8} more skipped")
    for w in rep.warnings:
        print(f"  WARNING: {w}")
    dtm = [w for w in rep.written if w["product"].startswith("dtm")]
    if dtm:
        best = min(dtm, key=lambda w: w["res_m"])
        print(f"\n  next:  python -m eclipse_sites.cli site --lat .. --lon .. "
              f"--dem {best['path']} --range 15000 --view view.png")


def cmd_calibrate(a):
    import numpy as np
    from pyproj import Transformer
    from .calibrate import (intrinsics_from_exif, extract_skyline, fit_camera,
                            canopy_from_residual, diagnostic_overlay,
                            below_bare_earth, canopy_profile, pixel_to_altaz)
    from .ephemeris import circumstances
    from .horizon import horizon_point
    from .search import read_dem, write_csv

    c = circumstances(a.lat, a.lon, 0.0, dt.date.fromisoformat(a.date))
    intr = intrinsics_from_exif(a.image, f_px=a.f_px, hfov_deg=a.hfov)
    print(f"  intrinsics  {intr.width}x{intr.height}  f_px {intr.f_px:.1f}  "
          f"({intr.source})")

    masks = []
    for m in a.mask:
        masks.append(tuple(int(v) for v in m.split(",")))
    sky = extract_skyline(a.image, masks=masks, sky_lum=a.sky_lum,
                          sky_sat=a.sky_sat)
    print(f"  skyline     {int(np.isfinite(sky).sum())}/{intr.width} columns")

    z, tr, crs, cell = read_dem(a.dem)
    x, y = Transformer.from_crs("EPSG:4326", crs, always_xy=True).transform(
        a.lon, a.lat)
    col, row = (~tr) * (x, y)
    ny, nx = z.shape
    if not (0 <= row < ny and 0 <= col < nx):
        sys.exit(f"  the camera position {a.lat:.5f}, {a.lon:.5f} is outside "
                 f"the DEM ({nx}x{ny} cells)")
    azs = list(np.arange(0.0, 360.0, a.cal_az_step))
    h = horizon_point(z, cell, row - 0.5, col - 0.5, azs,
                      max_range_m=a.range, eye_h=a.eye)
    dem_alt = [h[q].alt_deg for q in azs]
    dem_rng = [h[q].range_m for q in azs]
    short = sum(1 for q in azs if h[q].reach_m < a.range - cell)
    print(f"  horizon     {len(azs)} azimuths at {a.cal_az_step}deg from the DEM")
    if short:
        print(f"              {short} of them ran out of DEM before "
              f"{a.range:.0f} m")

    fit = fit_camera(sky, intr, azs, dem_alt, dem_rng,
                     bearing_prior=a.bearing_prior,
                     bearing_window=a.bearing_window,
                     anchor_min_range_m=a.anchor_min_range,
                     fit_roll=a.fit_roll, fit_focal=a.fit_focal)
    print(f"\n  BEARING     {fit.bearing_deg:.2f} deg")
    print(f"  PITCH       {fit.pitch_deg:+.2f} deg")
    if a.fit_roll:
        print(f"  ROLL        {fit.roll_deg:+.2f} deg")
    if a.fit_focal:
        print(f"  FOCAL       x{fit.f_scale:.4f}")
    print(f"  residual    {fit.rms_deg:.2f} deg rms over {fit.n_anchor} anchor "
          f"columns  (search cost {fit.cost:.2f})")
    if np.isfinite(fit.ambiguity_deg):
        print(f"\n  WARNING: a rival solution {fit.ambiguity_deg:.0f} deg away scores "
              f"{fit.rival_cost:.2f}.\n  The skyline is too repetitive to pin the "
              f"bearing. Give --bearing-prior\n  and a tighter --bearing-window, or "
              f"use a photo with more distinct\n  distant relief.")
    if fit.rms_deg > 2.0:
        print("\n  WARNING: high residual. Check the skyline trace with "
              "--out-overlay before\n  trusting the canopy numbers.")

    # a photographed skyline cannot sit below bare earth; if it does, the sky
    # test has eaten pale distant terrain and the trace is riding nearer trees
    frac, med, ncol = below_bare_earth(sky, intr, fit, azs, dem_alt, dem_rng,
                                       anchor_min_range_m=a.anchor_min_range)
    if ncol:
        print(f"\n  trace check   {frac * 100:.0f}% of {ncol} distant columns "
              f"sit below bare earth\n                median gap "
              f"{med:+.2f} deg (canopy should make this positive)")
        if frac > 0.15:
            print(f"  WARNING: a photograph cannot see ground below the "
                  f"ground. Pale distant\n  terrain is being read as sky, so "
                  f"the trace has dropped to nearer\n  treetops. Raise "
                  f"--sky-lum (now {a.sky_lum:g}) or lower --sky-sat (now "
                  f"{a.sky_sat:g}),\n  then check --out-overlay again.")

    # the most direct feedback there is: what the photo puts where the Sun will
    sun_az = c.sun_az_at_max
    p_alt, p_az = pixel_to_altaz(
        np.arange(intr.width, dtype=float), sky, intr, fit.bearing_deg,
        fit.pitch_deg, fit.roll_deg, fit.f_scale)
    near = np.isfinite(sky) & (np.abs((p_az - sun_az + 180) % 360 - 180) < 0.5)
    print(f"\n  at the Sun's bearing {sun_az:.1f} deg:")
    if near.any():
        ph = float(np.median(p_alt[near]))
        dm = float(np.interp(sun_az, azs, dem_alt))
        rg = float(np.interp(sun_az, azs, dem_rng))
        print(f"    photographed skyline  {ph:+.2f} deg   (includes canopy)")
        print(f"    DEM bare earth        {dm:+.2f} deg   "
              f"(forms at {rg:.0f} m)")
        print(f"    difference            {ph - dm:+.2f} deg")
        if rg > 5000:
            print(f"    the horizon there forms at {rg / 1000:.1f} km, too far "
                  f"for the difference\n    to be canopy. It is fit error, "
                  f"haze or trace error, and is not\n    written to the "
                  f"profile.")
    else:
        print(f"    not in frame; the photograph does not cover that bearing")

    eye_elev = float(z[int(row), int(col)]) + a.eye
    rows_out = canopy_from_residual(sky, intr, fit, azs, dem_alt, dem_rng,
                                    eye_elev_m=eye_elev, az_bin=2.0)
    if a.out_csv:
        write_csv(a.out_csv, rows_out)
        print(f"\n  wrote {a.out_csv}")
    if a.out_canopy:
        try:
            prof = canopy_profile(rows_out, fit)
        except ValueError as e:
            sys.exit(f"\n  refusing to write a canopy profile: {e}")
        m = prof["meta"]
        with open(a.out_canopy, "w") as f:
            json.dump(prof, f, indent=2)
        print(f"\n  wrote {a.out_canopy}")
        print(f"    {m['bins_kept']} bins kept, dropped "
              + ", ".join(f"{v} {k.replace('_', ' ')}"
                          for k, v in m["bins_dropped"].items() if v))
        if m["bins_kept"]:
            print(f"    covers az {m['az_min']:.0f} to {m['az_max']:.0f}; "
                  f"apply with  site --canopy-profile {a.out_canopy}")
        else:
            print(f"    nothing survived the filters. Every horizon in frame "
                  f"forms too far\n    away for its residual to be canopy, so "
                  f"there is nothing to apply.")
    if a.out_overlay:
        diagnostic_overlay(a.image, sky, intr, fit, azs, dem_alt, dem_rng,
                           a.out_overlay, a.anchor_min_range)
        print(f"  wrote {a.out_overlay}")

    big = [r for r in rows_out if abs(r["extra_height_m"]) > 3][:12]
    if big:
        print(f"\n  {'az':>6} {'resid':>7} {'range':>8} {'extra':>8}")
        for r in big:
            print(f"  {r['az_deg']:>6.1f} {r['residual_deg']:>+7.2f} "
                  f"{r['horizon_range_m']:>7d}m {r['extra_height_m']:>+7.1f}m")


def _skyline_floor(a):
    """Measured skyline altitude per azimuth, as a floor on the DEM horizon.

    Returns a callable giving degrees, or -90 where the photograph says nothing.
    The photograph is applied as an angle, never as metres of canopy along the
    ray: canopy_m lands on the sample 2.5 m from the eye as well as the far
    ones, and 20 m of tree at 2.5 m is 83 degrees of imaginary wall.
    """
    path = getattr(a, "canopy_profile", None)
    if not path:
        return lambda az: -90.0
    import json
    from .calibrate import profile_to_azimuths
    with open(path) as f:
        prof = json.load(f)
    meta = prof.get("meta", {})
    n = meta.get("bins_kept", len(prof.get("canopy", [])))
    print(f"  canopy      {n} measured bins from {path}"
          + (f", az {meta['az_min']:.0f} to {meta['az_max']:.0f}"
             if meta.get("az_min") is not None else ""))
    if meta.get("rms_deg") is not None:
        print(f"              camera fitted to {meta['rms_deg']:.2f} deg rms; "
              f"canopy carries that error")
    print(f"              outside that arc, --canopy {a.canopy:g} m applies")
    return lambda az: float(profile_to_azimuths(
        [az], prof, fallback=-90.0, key="skyline_alt_deg")[0])


def cmd_coverage(a):
    """What radius can this site honestly be judged on?"""
    import numpy as np
    from pyproj import Transformer
    from .horizon import coverage_point, edge_reach_point
    from .search import read_dem

    z, tr, crs, cell = read_dem(a.dem)
    x, y = Transformer.from_crs("EPSG:4326", crs, always_xy=True).transform(
        a.lon, a.lat)
    col, row = (~tr) * (x, y)
    ny, nx = z.shape
    print(f"  raster      {nx}x{ny} cells at {cell:g} m "
          f"({nx * cell / 1000:.1f} x {ny * cell / 1000:.1f} km)")
    if not (0 <= row < ny and 0 <= col < nx):
        sys.exit(f"  {a.lat:.5f}, {a.lon:.5f} is outside the raster entirely")

    # free: four subtractions, no data read
    box = edge_reach_point(z.shape, cell, row - 0.5, col - 0.5)
    print(f"  radius asked  {a.range / 1000:.1f} km")
    print(f"  bounds allow  {box / 1000:.1f} km   (nearest raster edge)")
    if box < a.range:
        print(f"  the raster cannot answer a {a.range / 1000:.0f} km question "
              f"from here, whatever\n  the data inside it looks like.")

    azs = list(np.arange(0.0, 360.0, a.az_step))
    step = a.step or cell
    try:
        cov = coverage_point(z, cell, row - 0.5, col - 0.5, azs,
                             max_range_m=a.range, step_m=a.step)
    except ValueError:
        # invariant 1b: an observer on nodata has no elevation of its own, so
        # nothing can be computed from here at all. That is the answer, not a
        # crash, and it is exactly what this command exists to report.
        sys.exit(f"\n  the observer cell itself is nodata.\n"
                 f"  VERDICT: zero coverage. No horizon can be computed from "
                 f"this point.\n  Move to a cell with ground, or fill the gap "
                 f"in the DEM.")
    con = np.array([cov[q].contiguous_m for q in azs])
    rch = np.array([cov[q].reach_m for q in azs])
    # the last sample lands up to one step short of the range asked for, so a
    # fully covered ray need not reach it exactly
    full = int((con >= a.range - max(step, cell)).sum())
    print(f"\n  {len(azs)} azimuths marched at {step:g} m")
    print(f"  full coverage {full} of {len(azs)} azimuths")
    print(f"  contiguous    min {con.min() / 1000:.1f} km   "
          f"median {np.median(con) / 1000:.1f} km   "
          f"max {con.max() / 1000:.1f} km")
    holed = int((rch > con + step).sum())
    if holed:
        print(f"  {holed} azimuth{'s' if holed > 1 else ''} ha"
              f"{'ve' if holed > 1 else 's'} ground beyond a gap. That ground "
              f"was seen but the\n  gap before it was not, so the ray is only "
              f"good to its contiguous range.")

    worst = int(np.argmin(con))
    print(f"  worst bearing {azs[worst]:.1f} deg at "
          f"{con[worst] / 1000:.1f} km")
    if full == len(azs):
        print(f"\n  VERDICT: covered to {a.range / 1000:.0f} km in every "
              f"direction.")
    else:
        print(f"\n  VERDICT: NOT covered to {a.range / 1000:.0f} km. Any "
              f"horizon from here is\n  only known to "
              f"{con.min() / 1000:.1f} km in the worst direction. A low "
              f"horizon\n  is not evidence of a clear view beyond that.")

    if a.json:
        import json
        with open(a.json, "w") as f:
            json.dump({"lat": a.lat, "lon": a.lon,
                       "range_asked_m": a.range,
                       "edge_reach_m": box,
                       "contiguous_min_m": float(con.min()),
                       "contiguous_median_m": float(np.median(con)),
                       "azimuths_full": full, "azimuths": len(azs),
                       "complete": full == len(azs)}, f, indent=2)
        print(f"  wrote {a.json}")


def main(argv=None):
    p = argparse.ArgumentParser(prog="eclipse-sites",
                                description="Eclipse visibility from terrain")
    sub = p.add_subparsers(dest="cmd", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--date", default=DEFAULT_DATE)
    common.add_argument("--tz", type=float, default=DEFAULT_TZ,
                        help="hours east of UTC, for display only")
    common.add_argument("--eye", type=float, default=1.6, help="eye height, m")
    common.add_argument("--range", type=float, default=40000.0,
                        help="max horizon search range, m")
    common.add_argument("--canopy", type=float, default=0.0,
                        help="metres of tree cover to add to ground heights")

    s = sub.add_parser("site", parents=[common], help="one location")
    s.add_argument("--lat", type=float, required=True)
    s.add_argument("--lon", type=float, required=True)
    s.add_argument("--elev", type=float, default=0.0)
    s.add_argument("--dem")
    s.add_argument("--json")
    s.add_argument("--canopy-profile", dest="canopy_profile",
                   help="canopy JSON from calibrate --out-canopy; applies "
                        "measured cover per azimuth instead of a flat --canopy")
    s.add_argument("--view", help="render the predicted view to this PNG")
    s.add_argument("--fov", type=float, help="field of view in degrees")
    s.set_defaults(func=cmd_site)

    r = sub.add_parser("search", parents=[common], help="whole region")
    r.add_argument("--dem", required=True, help="coarse DEM, projected, metres")
    r.add_argument("--fine", help="fine DEM for confirming the shortlist")
    r.add_argument("--az-step", type=float, default=2.5)
    r.add_argument("--top", type=int, default=25)
    r.add_argument("--min-sep", type=int, default=20,
                   help="minimum separation between results, in cells")
    r.add_argument("--out-raster")
    r.add_argument("--out-csv")
    r.set_defaults(func=cmd_search)

    k = sub.add_parser("calibrate", parents=[common],
                       help="fit a photo to the terrain and measure canopy")
    k.add_argument("--image", required=True)
    k.add_argument("--lat", type=float, required=True)
    k.add_argument("--lon", type=float, required=True)
    k.add_argument("--dem", required=True)
    k.add_argument("--mask", action="append", default=[],
                   help="x0,y0,x1,y1 region to ignore; repeatable")
    k.add_argument("--f-px", type=float)
    k.add_argument("--hfov", type=float)
    k.add_argument("--bearing-prior", type=float)
    k.add_argument("--bearing-window", type=float, default=180.0)
    k.add_argument("--anchor-min-range", type=float, default=2000.0)
    k.add_argument("--az-step", dest="cal_az_step", type=float, default=0.5)
    k.add_argument("--fit-roll", action="store_true")
    k.add_argument("--fit-focal", action="store_true")
    k.add_argument("--sky-lum", dest="sky_lum", type=float, default=115.0,
                   help="luminance above which a pixel may be sky; raise it if "
                        "pale distant ridges are being eaten")
    k.add_argument("--sky-sat", dest="sky_sat", type=float, default=0.32,
                   help="saturation below which a pixel may be sky; lower it "
                        "for the same reason")
    k.add_argument("--out-csv")
    k.add_argument("--out-canopy",
                   help="write a canopy profile for site --canopy-profile")
    k.add_argument("--out-overlay")
    k.set_defaults(func=cmd_calibrate)

    v = sub.add_parser("coverage",
                       help="what radius can this site honestly be judged on")
    v.add_argument("--lat", type=float, required=True)
    v.add_argument("--lon", type=float, required=True)
    v.add_argument("--dem", required=True)
    v.add_argument("--range", type=float, default=40000.0,
                   help="radius you mean to ask for, m")
    v.add_argument("--az-step", dest="az_step", type=float, default=1.0)
    v.add_argument("--step", type=float,
                   help="march step, m; defaults to one cell")
    v.add_argument("--json")
    v.set_defaults(func=cmd_coverage)

    q = sub.add_parser("prepare",
                       help="merge and reproject a folder of downloaded tiles")
    q.add_argument("--in", dest="indir", required=True,
                   help="folder of downloaded files; zips and subfolders fine")
    q.add_argument("--out", dest="outdir", required=True,
                   help="directory to write the prepared rasters into")
    q.add_argument("--src-crs", help="CRS of inputs that carry none, "
                                     "e.g. EPSG:25829 for UTM zone 29")
    q.add_argument("--target-crs", help="defaults to the commonest input CRS")
    q.add_argument("--res", type=float, help="override output pixel size")
    q.add_argument("--only", action="append",
                   help="limit to a product, e.g. dtm05; repeatable")
    q.add_argument("--max-cells", type=float, default=4e8,
                   help="skip a product whose mosaic would exceed this many "
                        "cells; merge holds it all in memory")
    q.set_defaults(func=cmd_prepare)

    a = p.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    main()
