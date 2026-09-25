"""Fit a photograph to the terrain, then read the trees out of the residual.

The order matters. Earlier attempts asked a photo for the camera angles and the
obstruction height at the same time, which is one equation with two unknowns:
the bearing error silently absorbs the obstruction and you get a confident wrong
answer. Here the DEM fixes the geometry first.

  1. trace the skyline in the image
  2. fit camera bearing and pitch against the DEM, using ONLY sectors where the
     DEM is trustworthy, meaning the horizon forms far away where canopy is
     worth a fraction of a degree
  3. everything left over is what the survey does not know about: trees,
     buildings, anything newer than the LiDAR flight

Step 2's anchor choice is the whole ballgame. Fit over wooded near ground and
the canopy is absorbed into the camera angles, destroying the measurement.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

__all__ = ["Intrinsics", "CameraFit", "intrinsics_from_exif", "extract_skyline",
           "pixel_to_altaz", "fit_camera", "canopy_from_residual"]

R_EFF = 6371000.0 / 0.87


# ------------------------------------------------------------------ camera

@dataclass
class Intrinsics:
    width: int
    height: int
    f_px: float
    cx: float
    cy: float
    source: str = "exif"


def intrinsics_from_exif(path, native_aspect=(4, 3), f_px=None, hfov_deg=None):
    """Focal length in pixels, from EXIF unless overridden.

    A phone reports 35 mm equivalent focal length for its *native* frame. If the
    image is a crop of that frame the diagonal no longer matches, so the native
    aspect ratio is needed. When in doubt let the fit refine it: pass
    fit_focal=True to fit_camera.
    """
    from PIL import Image
    from PIL.ExifTags import TAGS

    im = Image.open(path)
    W, H = im.size
    if f_px is not None:
        return Intrinsics(W, H, float(f_px), W / 2, H / 2, "supplied f_px")
    if hfov_deg is not None:
        f = (W / 2) / math.tan(math.radians(hfov_deg) / 2)
        return Intrinsics(W, H, f, W / 2, H / 2, "supplied hfov")

    f35 = None
    ex = im.getexif() or {}
    ifds = [ex]
    try:
        ifds.append(ex.get_ifd(0x8769))     # phones put it in the Exif sub-IFD
    except Exception:
        pass
    for d in ifds:
        for k, v in (d or {}).items():
            if TAGS.get(k) == "FocalLengthIn35mmFilm" and v:
                f35 = float(v)
    if not f35:
        raise ValueError(
            "no FocalLengthIn35mmFilm in EXIF; pass --f-px or --hfov")

    # diagonal field of view of a 35 mm frame at this equivalent focal length
    diag_fov = 2 * math.atan(21.63 / f35)
    aw, ah = native_aspect
    native_h = W * ah / aw            # native frame height in this image's pixels
    native_diag = math.hypot(W, native_h)
    f = (native_diag / 2) / math.tan(diag_fov / 2)
    return Intrinsics(W, H, f, W / 2, H / 2, f"exif f35={f35:g}mm")


def pixel_to_altaz(x, y, intr, bearing_deg, pitch_deg, roll_deg=0.0,
                   f_scale=1.0):
    """Exact rectilinear mapping from pixel to (altitude, azimuth) in degrees."""
    f = intr.f_px * f_scale
    u = np.asarray(x, float) - intr.cx
    v = np.asarray(y, float) - intr.cy
    if roll_deg:
        r = math.radians(roll_deg)
        u, v = u * math.cos(r) + v * math.sin(r), -u * math.sin(r) + v * math.cos(r)

    b, p = math.radians(bearing_deg), math.radians(pitch_deg)
    sb, cb, sp, cp = math.sin(b), math.cos(b), math.sin(p), math.cos(p)
    fwd = np.array([sb * cp, cb * cp, sp])
    right = np.array([cb, -sb, 0.0])
    up = np.array([-sb * sp, -cb * sp, cp])

    d = (u[..., None] * right) + (-v[..., None] * up) + (f * fwd)
    d /= np.linalg.norm(d, axis=-1, keepdims=True)
    alt = np.degrees(np.arcsin(np.clip(d[..., 2], -1, 1)))
    az = np.degrees(np.arctan2(d[..., 0], d[..., 1])) % 360.0
    return alt, az


# ---------------------------------------------------------------- skyline

def extract_skyline(path, masks=(), sky_lum=115.0, sky_sat=0.32, run=20,
                    smooth=15):
    """Skyline row index per image column.

    masks: iterable of (x0, y0, x1, y1) boxes to ignore, for parasols, fences,
    poles and anything else in the foreground. Columns that are entirely masked
    or where no skyline is found come back as NaN.
    """
    from PIL import Image
    a = np.asarray(Image.open(path).convert("RGB"), dtype=np.float32)
    H, W, _ = a.shape
    lum = a.mean(axis=2)
    mx, mn = a.max(axis=2), a.min(axis=2)
    sat = (mx - mn) / np.maximum(mx, 1e-6)
    sky = (lum > sky_lum) & (sat < sky_sat)

    blocked = np.zeros((H, W), bool)
    for x0, y0, x1, y1 in masks:
        blocked[max(0, y0):min(H, y1), max(0, x0):min(W, x1)] = True

    prof = np.full(W, np.nan)
    for x in range(W):
        if blocked[:, x].all():
            continue
        col = sky[:, x] | blocked[:, x]     # masked pixels never end the sky run
        cnt = 0
        for y in range(H):
            if col[y]:
                cnt = 0
            else:
                cnt += 1
                if cnt >= run:
                    prof[x] = y - run + 1
                    break

    if smooth and smooth > 1:               # median filter, kills speckle
        out = prof.copy()
        h = smooth // 2
        for x in range(W):
            w = prof[max(0, x - h):x + h + 1]
            w = w[np.isfinite(w)]
            if w.size:
                out[x] = np.median(w)
        prof = out
    return prof


# -------------------------------------------------------------------- fit

def _ring(az, alt, rng):
    """Sort a horizon fan by azimuth, closing the seam if it is a full circle.

    A 0 to 359.5 fan has a 0.5 degree hole in it as far as np.interp is
    concerned, and everything in that hole clamps to an endpoint instead of
    wrapping.
    """
    az = np.asarray(az, float)
    o = np.argsort(az)
    az, alt, rng = az[o], np.asarray(alt, float)[o], np.asarray(rng, float)[o]
    if az[-1] - az[0] >= 350.0:
        az = np.concatenate([az[-1:] - 360.0, az, az[:1] + 360.0])
        alt = np.concatenate([alt[-1:], alt, alt[:1]])
        rng = np.concatenate([rng[-1:], rng, rng[:1]])
    return az, alt, rng


@dataclass
class CameraFit:
    bearing_deg: float
    pitch_deg: float
    roll_deg: float
    f_scale: float
    rms_deg: float          # true angular rms over the anchor columns
    n_anchor: int
    intr: Intrinsics
    ambiguity_deg: float = float("nan")   # bearing gap to the next rival basin
    rival_cost: float = float("nan")
    cost: float = float("nan")            # rms_deg penalised by coverage


def _sweep(cost2, b_lo, b_hi, p_lo, p_hi, b_step=1.0, p_step=1.0):
    """Dense 2D scan. Deliberately not clever: a coarse-to-fine search alone
    latches onto false minima because the cost surface over bearing is bumpy
    when a ridgeline has repeating structure."""
    bs = np.arange(b_lo, b_hi + b_step / 2, b_step)
    ps = np.arange(p_lo, p_hi + p_step / 2, p_step)
    surf = np.full((bs.size, ps.size), np.inf)
    for i, b in enumerate(bs):
        for j, p in enumerate(ps):
            surf[i, j] = cost2(b, p)
    return bs, ps, surf


def _refine(cost, centres, spans, rounds=5, steps=9):
    best = (cost(centres), list(centres))
    for _ in range(rounds):
        axes = [np.linspace(c - s, c + s, steps) if s > 0 else np.array([c])
                for c, s in zip(centres, spans)]
        for idx in np.ndindex(*[len(a) for a in axes]):
            p = [ax[i] for ax, i in zip(axes, idx)]
            c = cost(p)
            if c < best[0]:
                best = (c, list(p))
        centres = best[1]
        spans = [s / (steps // 2) for s in spans]
    return best[1], best[0]


def fit_camera(skyline_y, intr, dem_az, dem_alt, dem_range,
               bearing_prior=None, bearing_window=180.0,
               pitch_window=35.0, anchor_min_range_m=2000.0,
               anchor_sectors=None, fit_roll=False, fit_focal=False,
               col_step=8, min_anchor_frac=0.25):
    """Solve camera bearing and pitch against a DEM horizon profile.

    dem_az/dem_alt/dem_range describe the true horizon. Anchoring defaults to
    azimuths where the horizon forms beyond anchor_min_range_m, because 20 m of
    canopy is worth 0.33 deg at 3.5 km but 14 deg at 80 m. Override with
    anchor_sectors=[(az_lo, az_hi), ...] if you know the ground.
    """
    dem_az, dem_alt, dem_range = _ring(dem_az, dem_alt, dem_range)

    cols = np.arange(0, intr.width, col_step)
    ys = skyline_y[cols]
    ok = np.isfinite(ys)
    cols, ys = cols[ok], ys[ok]
    if cols.size < 20:
        raise ValueError("too few usable skyline columns to fit")

    def anchored(az):
        if anchor_sectors:
            m = np.zeros(az.shape, bool)
            for lo, hi in anchor_sectors:
                # hi < lo means the sector runs through north
                m |= ((az >= lo) | (az <= hi)) if hi < lo else \
                     ((az >= lo) & (az <= hi))
            return m
        return np.interp(az, dem_az, dem_range) >= anchor_min_range_m

    lo_az = dem_az[0]
    hi_az = dem_az[-1]
    n_cols = cols.size
    min_cols = max(20, int(min_anchor_frac * n_cols))

    def cost(p):
        b, pit = p[0], p[1]
        roll = p[2] if fit_roll else 0.0
        fs = p[3] if fit_focal else 1.0
        alt, az = pixel_to_altaz(cols, ys, intr, b, pit, roll, fs)
        m = (az >= lo_az) & (az <= hi_az) & anchored(az)
        n = int(m.sum())
        # A fit is not allowed to buy a low residual by explaining almost
        # nothing. Without this the solver slides the frame until only a
        # handful of columns overlap the anchor sectors and fits those
        # perfectly, which is how you get a 0.05 deg rms that means nothing.
        if n < min_cols:
            return 1e6
        r = np.interp(az[m], dem_az, dem_alt) - alt[m]
        return float(np.sqrt(np.mean(r * r))) * math.sqrt(n_cols / n)

    b0 = bearing_prior if bearing_prior is not None else float(np.mean(dem_az))

    def cost2(b, p):
        return cost([b, p] + ([0.0] if fit_roll else []) +
                    ([1.0] if fit_focal else []))

    bs, ps, surf = _sweep(cost2, b0 - bearing_window, b0 + bearing_window,
                          -pitch_window, pitch_window)
    flat = surf.min(axis=1)
    i_best = int(np.argmin(flat))
    b_best, p_best = bs[i_best], ps[int(np.argmin(surf[i_best]))]

    # Is there a genuinely different bearing that fits nearly as well? If so the
    # ridgeline is too repetitive to pin down and the user needs to know.
    # Compare on the circle: a rival 355 deg away is 5 deg away, not distant.
    amb, rival = float("nan"), float("nan")
    gap = np.abs((bs - b_best + 180.0) % 360.0 - 180.0)
    far = gap > 5.0
    if far.any() and np.isfinite(flat[far]).any():
        j = int(np.argmin(np.where(far, flat, np.inf)))
        if flat[j] < flat[i_best] * 2.0:
            amb, rival = float(gap[j]), float(flat[j])

    centres = [b_best, p_best]
    spans = [1.5, 1.5]
    if fit_roll:
        centres.append(0.0)
        spans.append(10.0)
    if fit_focal:
        centres.append(1.0)
        spans.append(0.2)
    p, best_cost = _refine(cost, centres, spans)

    b, pit = p[0] % 360.0, p[1]
    roll = p[2] if fit_roll else 0.0
    fs = p[3] if fit_focal else 1.0
    alt, az = pixel_to_altaz(cols, ys, intr, b, pit, roll, fs)
    m = (az >= lo_az) & (az <= hi_az) & anchored(az)
    n = int(m.sum())
    # Report the true angular rms, not the coverage-penalised search cost. The
    # penalty is what stops the solver buying a low residual with a sliver of
    # overlap, but it is not a number of degrees and must not be shown as one.
    if n:
        r = np.interp(az[m], dem_az, dem_alt) - alt[m]
        rms = float(np.sqrt(np.mean(r * r)))
    else:
        rms = float("nan")
    return CameraFit(b, pit, roll, fs, rms, n, intr, amb, rival, best_cost)


# ----------------------------------------------------------------- canopy

def canopy_from_residual(skyline_y, intr, fit, dem_az, dem_alt, dem_range,
                         eye_elev_m, col_step=8, az_bin=2.0):
    """Height of whatever the DEM is missing, in metres, binned by azimuth.

    Residual angle plus the range at which the horizon forms gives a height
    directly. Sanity check the output: a large residual attributed to a distant
    horizon yields an absurd height, which means the fit is wrong rather than
    the trees are tall.
    """
    dem_az, dem_alt, dem_range = _ring(dem_az, dem_alt, dem_range)

    cols = np.arange(0, intr.width, col_step)
    ys = skyline_y[cols]
    ok = np.isfinite(ys)
    cols, ys = cols[ok], ys[ok]
    alt, az = pixel_to_altaz(cols, ys, intr, fit.bearing_deg, fit.pitch_deg,
                             fit.roll_deg, fit.f_scale)

    inside = (az >= dem_az[0]) & (az <= dem_az[-1])
    alt, az = alt[inside], az[inside]
    d = np.interp(az, dem_az, dem_range)
    g = np.interp(az, dem_az, dem_alt)
    valid = d > 0

    # height of the photographed skyline above the eye, at that range
    h_photo = eye_elev_m + d * np.tan(np.radians(alt)) + d * d / (2 * R_EFF)
    h_dem = eye_elev_m + d * np.tan(np.radians(g)) + d * d / (2 * R_EFF)
    extra = h_photo - h_dem

    rows = []
    lo = math.floor(az.min() / az_bin) * az_bin
    while lo < az.max():
        m = valid & (az >= lo) & (az < lo + az_bin)
        if m.sum() >= 3:
            rows.append({
                "az_deg": round(lo + az_bin / 2, 2),
                "photo_alt_deg": round(float(np.median(alt[m])), 2),
                "dem_alt_deg": round(float(np.median(g[m])), 2),
                "residual_deg": round(float(np.median(alt[m] - g[m])), 2),
                "horizon_range_m": round(float(np.median(d[m]))),
                "extra_height_m": round(float(np.median(extra[m])), 1),
            })
        lo += az_bin
    return rows


def diagnostic_overlay(image_path, skyline_y, intr, fit, dem_az, dem_alt,
                       dem_range, out_path, anchor_min_range_m=2000.0):
    """Photo with the traced skyline, the DEM prediction, and anchor sectors."""
    from PIL import Image, ImageDraw
    im = Image.open(image_path).convert("RGB")
    d = ImageDraw.Draw(im)
    cols = np.arange(0, intr.width, 4)
    ys = skyline_y[cols]

    # traced skyline in yellow
    pts = [(int(x), int(y)) for x, y in zip(cols, ys) if np.isfinite(y)]
    if len(pts) > 1:
        d.line(pts, fill=(255, 220, 60), width=5)

    # where the DEM says the skyline should be, given this fit
    rows = np.arange(0, intr.height, 2.0)
    pred = []
    for x in cols[::2]:
        alt, az = pixel_to_altaz(np.full(rows.shape, float(x)), rows, intr,
                                 fit.bearing_deg, fit.pitch_deg, fit.roll_deg,
                                 fit.f_scale)
        h = np.interp(az, dem_az, dem_alt, left=np.nan, right=np.nan)
        s = np.where(np.diff(np.sign(alt - h)) != 0)[0]
        if s.size:
            pred.append((int(x), int(rows[s[0]])))
    for i in range(0, len(pred) - 1):
        d.line([pred[i], pred[i + 1]], fill=(120, 220, 255), width=4)

    # green band marks the sectors the fit was anchored on
    for x in cols:
        _, az = pixel_to_altaz(np.array([float(x)]), np.array([intr.cy]),
                               intr, fit.bearing_deg, fit.pitch_deg,
                               fit.roll_deg, fit.f_scale)
        r = np.interp(az[0], dem_az, dem_range, left=0, right=0)
        if r >= anchor_min_range_m:
            d.rectangle([x, intr.height - 40, x + 4, intr.height],
                        fill=(80, 230, 130))
    im.save(out_path, quality=90)
    return out_path


# ------------------------------------------------------------ trace check

def below_bare_earth(skyline_y, intr, fit, dem_az, dem_alt, dem_range,
                     anchor_min_range_m=2000.0, tol_deg=0.3, col_step=8):
    """How much of the traced skyline falls below bare earth, where it cannot.

    A photograph cannot show ground lower than the ground. Canopy, buildings
    and haze only ever add height, so a column sitting below the DEM by more
    than tol_deg is a defect rather than a measurement. The usual cause is the
    sky test eating pale distant terrain: a hazy ridge reads as sky, the trace
    drops to the nearer treetops behind it, and every angle in that sector is
    wrong. Raise --sky-lum or lower --sky-sat until this clears.

    Deliberately a check on the trace, run after the fit and never folded into
    it. Penalising negative residuals inside the cost would push pitch up until
    the broken columns looked physical, which inflates the apparent canopy on
    every sound column and corrupts the very numbers this module exists to
    produce. The constraint is only meaningful once the trace is believed.

    Returns (fraction_below, median_gap_deg, n_columns) over distant sectors.
    """
    dem_az, dem_alt, dem_range = _ring(dem_az, dem_alt, dem_range)
    cols = np.arange(0, intr.width, col_step)
    ys = skyline_y[cols]
    ok = np.isfinite(ys)
    cols, ys = cols[ok], ys[ok]
    if cols.size == 0:
        return float("nan"), float("nan"), 0
    alt, az = pixel_to_altaz(cols, ys, intr, fit.bearing_deg, fit.pitch_deg,
                             fit.roll_deg, fit.f_scale)
    inside = (az >= dem_az[0]) & (az <= dem_az[-1])
    alt, az = alt[inside], az[inside]
    if az.size == 0:
        return float("nan"), float("nan"), 0
    far = np.interp(az, dem_az, dem_range) >= anchor_min_range_m
    if not far.any():
        return float("nan"), float("nan"), 0
    gap = alt[far] - np.interp(az[far], dem_az, dem_alt)
    return (float((gap < -tol_deg).mean()), float(np.median(gap)),
            int(far.sum()))


# ---------------------------------------------------------- canopy profile

def canopy_profile(rows, fit, max_range_m=5000.0, max_height_m=60.0,
                   min_residual_deg=-0.3, max_rms_deg=2.0):
    """Reduce measured canopy rows to the ones that can honestly be applied.

    A residual is only canopy where the horizon forms near. The conversion from
    angle to metres scales with range, so at 18.7 km a fifth of a degree is 65 m
    of imaginary tree; bins beyond max_range_m are dropped rather than believed.
    Small negative extras are noise and clamp to zero. A large negative one
    means the trace or the fit is wrong, so the bin goes. Heights above
    max_height_m go too: an implausible height is a broken measurement, not a
    tall tree, and clamping it to the cap would invent cover nobody saw.

    Refuses outright on a poor fit. A canopy number carries the camera's errors
    with it, and there is no honest way to publish one from a pose you do not
    believe.
    """
    if not np.isfinite(fit.rms_deg) or fit.rms_deg > max_rms_deg:
        raise ValueError(
            f"fit residual is {fit.rms_deg:.2f} deg rms, above the "
            f"{max_rms_deg:.2f} deg limit for writing a canopy profile. "
            f"Fix the trace or the pose first; canopy inherits both.")

    kept, dropped = [], {"horizon_too_far": 0, "residual_negative": 0,
                         "height_implausible": 0}
    for r in rows:
        if r["horizon_range_m"] > max_range_m:
            dropped["horizon_too_far"] += 1
            continue
        if r["residual_deg"] < min_residual_deg:
            dropped["residual_negative"] += 1
            continue
        if r["extra_height_m"] > max_height_m:
            dropped["height_implausible"] += 1
            continue
        kept.append({"az_deg": r["az_deg"],
                     "canopy_m": round(max(float(r["extra_height_m"]), 0.0), 1),
                     # the angle is what gets applied; the metres are for reading
                     "skyline_alt_deg": max(float(r["photo_alt_deg"]),
                                            float(r["dem_alt_deg"])),
                     "dem_alt_deg": r["dem_alt_deg"],
                     "horizon_range_m": r["horizon_range_m"]})

    azs = [k["az_deg"] for k in kept]
    return {
        "meta": {
            "bearing_deg": round(float(fit.bearing_deg), 2),
            "pitch_deg": round(float(fit.pitch_deg), 2),
            "rms_deg": round(float(fit.rms_deg), 3),
            "n_anchor": int(fit.n_anchor),
            "az_min": min(azs) if azs else None,
            "az_max": max(azs) if azs else None,
            "bins_kept": len(kept),
            "bins_dropped": dropped,
            "filters": {"max_range_m": max_range_m,
                        "max_height_m": max_height_m,
                        "min_residual_deg": min_residual_deg},
        },
        "canopy": kept,
    }


def profile_to_azimuths(azimuths_deg, profile, fallback=0.0, max_gap_deg=4.0,
                        key="canopy_m"):
    """One measured field per azimuth, from a written profile.

    key selects the field. "skyline_alt_deg" is the one to apply to a horizon:
    it is the altitude the photograph actually saw, trees included. Do NOT feed
    "canopy_m" to horizon_point as its canopy_m argument. That adds the height
    to every ground sample along the ray including the one 2.5 m from the eye,
    where 20 m of tree subtends 83 degrees, and the horizon comes back as a wall.
    The photograph measured an angle; converting it to metres and back is what
    breaks it.

    Linear interpolation inside the arc the photograph covered, fallback
    outside it. A profile is an arc and not a ring: it describes only what was
    photographed, so unlike a full-circle horizon fan there is no seam to
    close. The arc is located by its largest gap, so one crossing north is
    handled without pretending the numbers are a sorted sequence.

    max_gap_deg is what stops a dropped bin from coming back. canopy_profile
    discards bins whose horizon forms too far away to convert to a height, and
    plain interpolation would rebuild them out of their neighbours, which is
    the discarded measurement returning under another name. An azimuth further
    than max_gap_deg from any surviving bin gets the fallback instead.
    """
    entries = profile.get("canopy", []) if isinstance(profile, dict) else profile
    q = np.asarray(azimuths_deg, dtype=float)
    out = np.full(q.shape, float(fallback))
    if not entries:
        return out

    az_p = np.array([e["az_deg"] for e in entries], dtype=float)
    can = np.array([e[key] for e in entries], dtype=float)
    if az_p.size == 1:
        return out

    s = np.sort(az_p)
    gaps = np.diff(np.concatenate([s, [s[0] + 360.0]]))
    start = s[(int(np.argmax(gaps)) + 1) % s.size]   # arc begins after the gap
    rel = (az_p - start) % 360.0
    order = np.argsort(rel)
    rel, can = rel[order], can[order]

    qr = (q - start) % 360.0
    # the test is the width of the hole the query falls in, not its distance to
    # the nearest bin: a query 3 deg from a surviving bin is still inside an
    # 8 deg hole, and interpolating it rebuilds the bins that were discarded
    j = np.clip(np.searchsorted(rel, qr), 1, rel.size - 1)
    spanned = (rel[j] - rel[j - 1]) <= max_gap_deg
    on_bin = np.min(np.abs(qr[:, None] - rel[None, :]), axis=1) < 1e-9
    inside = (qr <= rel[-1]) & (spanned | on_bin)
    out[inside] = np.interp(qr[inside], rel, can)
    return out
