"""Draw what the eclipse will look like from a given point.

No photograph needed. The skyline is built from the DEM, so this works for any
coordinate, including places nobody has ever stood with a camera.

Ridges are shaded by the distance at which they form the horizon, near dark and
close, far pale and hazy, which is what aerial perspective actually does and
makes the result readable as a landscape rather than a chart.
"""
from __future__ import annotations

import datetime as dt
import math
import os

import numpy as np
from PIL import Image, ImageDraw, ImageFont

from .ephemeris import circumstances, geometry

__all__ = ["horizon_from_dem", "render_view"]

SKY_HI = (86, 116, 168)
SKY_LO = (206, 214, 226)
NEAR = (28, 38, 32)
FAR = (128, 146, 168)
AMBER = (255, 180, 60)
PEARL = (242, 237, 228)
UMBRA = (24, 30, 52)
UNSURVEYED = (96, 92, 104)   # hatch for frame outside the horizon fan


# Scalable faces, in preference order, per platform. A hardcoded Linux path
# used to fall through to load_default(), which ignores the size argument, so
# every label on macOS rendered at 10 px instead of the 26 to 42 asked for.
_FONT_DIRS = (
    "/usr/share/fonts/truetype/dejavu",             # Debian, Ubuntu
    "/usr/share/fonts/dejavu",                      # Fedora, Arch
    "/System/Library/Fonts/Supplemental",           # macOS
    "/System/Library/Fonts",                        # macOS
    "C:/Windows/Fonts",                             # Windows
)
_FONT_FILES = {
    False: ("DejaVuSans.ttf", "Helvetica.ttc", "Arial.ttf", "SFNS.ttf",
            "arial.ttf"),
    True: ("DejaVuSans-Bold.ttf", "Arial Bold.ttf", "Arial.ttf",
           "Helvetica.ttc", "arialbd.ttf"),
}
_font_cache = {}


def _font(size, bold=False):
    key = (int(size), bool(bold))
    if key in _font_cache:
        return _font_cache[key]
    for name in _FONT_FILES[bool(bold)]:
        for d in _FONT_DIRS:
            try:
                f = ImageFont.truetype(os.path.join(d, name), size)
                _font_cache[key] = f
                return f
            except OSError:
                continue
    try:
        f = ImageFont.load_default(size)     # Pillow >= 10.1 honours the size
    except TypeError:
        f = ImageFont.load_default()
    _font_cache[key] = f
    return f


def horizon_from_dem(dem_path, lat, lon, azimuths, eye_h=1.6, canopy_m=0.0,
                     max_range_m=40000.0):
    """Horizon altitude, range and reach at each azimuth. Needs rasterio.

    reach is how far along each ray the DEM actually had ground; short reach
    means the horizon is only known that far out.
    """
    from pyproj import Transformer
    from .horizon import horizon_point
    from .search import read_dem
    z, tr, crs, cell = read_dem(dem_path)
    x, y = Transformer.from_crs("EPSG:4326", crs, always_xy=True).transform(lon, lat)
    col, row = (~tr) * (x, y)
    h = horizon_point(z, cell, row - 0.5, col - 0.5, azimuths,
                      max_range_m=max_range_m, eye_h=eye_h, canopy_m=canopy_m)
    alt = np.array([h[a].alt_deg for a in azimuths], float)
    rng = np.array([h[a].range_m for a in azimuths], float)
    reach = np.array([h[a].reach_m for a in azimuths], float)
    return alt, rng, reach


def render_view(lat, lon, when, hor_az, hor_alt, hor_rng, out_path,
                elev_m=0.0, tz=0.0, size=(2200, 1100), fov=None,
                title=None, subtitle=None):
    """Render the eclipse over a synthetic skyline.

    hor_az/hor_alt/hor_rng: horizon samples, degrees and metres. They are
    interpolated across the frame, so a coarse fan still renders sensibly.
    """
    c = circumstances(lat, lon, elev_m, when)
    if c.kind == "none":
        raise ValueError("no eclipse at this location on this date")

    hor_az = np.asarray(hor_az, float)
    order = np.argsort(hor_az)
    hor_az, hor_alt, hor_rng = hor_az[order], np.asarray(hor_alt, float)[order], \
        np.asarray(hor_rng, float)[order]

    if not np.isfinite(hor_alt).any():
        raise ValueError("horizon profile is entirely undefined")

    # sun track over the whole eclipse
    track = []
    t = c.c1
    while t <= c.c4:
        alt, az, _, _, _ = geometry(lat, lon, elev_m, t)
        track.append((t, alt, az))
        t += dt.timedelta(minutes=1)

    # Unwrap the track, then slide it by whole turns onto the same branch as
    # the horizon fan. Without this an eclipse straddling due north puts the
    # Sun at 359 deg and the fan at 1 deg, and the frame tears in half.
    fan_mid = float(hor_az[0] + hor_az[-1]) / 2
    azs = np.degrees(np.unwrap(np.radians([a for _, _, a in track])))
    azs = azs + 360.0 * round((fan_mid - float(azs.mean())) / 360.0)
    track = [(t, v, a) for (t, v, _), a in zip(track, azs)]

    span = fov or max(60.0, (azs.max() - azs.min()) + 40.0)
    az0 = (azs.min() + azs.max()) / 2 - span / 2
    frame_mid = az0 + span / 2

    def A(a):
        """Put a raw compass azimuth on the frame's branch."""
        return a - 360.0 * round((a - frame_mid) / 360.0)

    def hor_at(a):
        """Fan altitude at one azimuth, NaN where the fan does not reach."""
        return float(np.interp(A(a), hor_az, hor_alt,
                               left=np.nan, right=np.nan))

    alt_top = max(max(a for _, a, _ in track), float(np.nanmax(hor_alt))) + 6
    alt_bot = min(-3.0, float(np.nanmin(hor_alt)) - 2)

    W, H = size
    X = lambda a: (A(a) - az0) / span * W
    Y = lambda v: H - (v - alt_bot) / (alt_top - alt_bot) * H

    im = Image.new("RGB", (W, H), SKY_HI)
    d = ImageDraw.Draw(im)
    horiz_y = Y(0.0)
    for y in range(H):
        f = min(1.0, max(0.0, y / max(1.0, horiz_y)))
        d.line([(0, y), (W, y)],
               fill=tuple(int(SKY_HI[i] + (SKY_LO[i] - SKY_HI[i]) * f ** 1.6)
                          for i in range(3)))

    # Skyline, shaded by how far away it is. Columns the fan does not cover are
    # left as sky and flagged: np.interp would otherwise clamp to the end of
    # the fan and draw a confident flat ridge across ground nobody surveyed.
    cols = np.arange(W)
    az_of = az0 + cols / W * span
    alt_of = np.interp(az_of, hor_az, hor_alt, left=np.nan, right=np.nan)
    rng_of = np.interp(az_of, hor_az, hor_rng, left=np.nan, right=np.nan)
    known = np.isfinite(alt_of)
    haze = np.clip(np.log10(np.maximum(np.nan_to_num(rng_of, nan=50.0), 50.0)
                            / 50.0) / 2.6, 0, 1)
    for x in range(W):
        if known[x]:
            col = tuple(int(NEAR[i] + (FAR[i] - NEAR[i]) * haze[x])
                        for i in range(3))
            d.line([(x, Y(alt_of[x])), (x, H)], fill=col)
        elif (x // 14) % 2 == 0:
            d.line([(x, H - 30), (x, H)], fill=UNSURVEYED)

    # altitude grid. Blend by hand: this is an RGB image, so a fourth channel
    # in the fill tuple is silently dropped and the line comes out solid white.
    f_small = _font(26)
    for v in range(0, int(alt_top) + 1, 5):
        y = Y(v)
        base = im.getpixel((W // 2, max(0, min(H - 1, int(y)))))
        d.line([(0, y), (W, y)],
               fill=tuple(int(b + (255 - b) * 0.16) for b in base), width=1)
        d.text((10, y - 30), f"{v}\u00b0", font=f_small, fill=(255, 255, 255))
    d.line([(0, horiz_y), (W, horiz_y)], fill=(255, 255, 255), width=2)

    # azimuth ticks, labelled as compass bearings rather than frame offsets
    step = 5 if span <= 80 else 10
    for a in range(int(math.ceil(az0 / step) * step), int(az0 + span) + 1, step):
        d.line([(X(a), H - 46), (X(a), H - 16)], fill=(255, 255, 255), width=2)
        d.text((X(a), H - 78), f"{a % 360}\u00b0", font=f_small, anchor="ma",
               fill=(255, 255, 255))

    # track, solid where the Sun is above the skyline
    pts = [(X(a), Y(v)) for _, v, a in track]
    vis = [v > hor_at(a) for _, v, a in track]   # NaN compares False: unproven
    for i in range(len(pts) - 1):
        if vis[i] and vis[i + 1]:
            d.line([pts[i], pts[i + 1]], fill=AMBER, width=8)
        else:
            if i % 6 < 3:
                d.line([pts[i], pts[i + 1]], fill=(150, 150, 150), width=4)

    f_lab = _font(30, True)
    for tt, v, a in track:
        loc = tt + dt.timedelta(hours=tz)
        if loc.minute % 15 or loc.second > 30:
            continue
        x, y = X(a), Y(v)
        on = v > hor_at(a)
        d.ellipse([x - 13, y - 13, x + 13, y + 13],
                  fill=AMBER if on else (140, 140, 140), outline=(40, 30, 10), width=2)
        d.text((x, y - 34), loc.strftime("%H:%M"), font=f_lab, anchor="ms",
               fill=(255, 255, 255), stroke_width=4, stroke_fill=(30, 30, 40))

    def mark(t, label, colour, dy):
        if t is None:
            return
        alt, az, _, _, _ = geometry(lat, lon, elev_m, t)
        x, y = X(az), Y(alt)
        d.ellipse([x - 26, y - 26, x + 26, y + 26], outline=colour, width=5)
        d.text((x, y + dy), label, font=_font(34, True), anchor="ms",
               fill=colour, stroke_width=5, stroke_fill=(20, 20, 30))

    mark(c.c1, "partial begins", (150, 210, 255), -60)
    mark(c.c4, "partial ends", (150, 210, 255), 78)
    if c.c2:
        alt, az, _, _, _ = geometry(lat, lon, elev_m, c.cmax)
        x, y = X(az), Y(alt)
        for r, w in ((92, 4), (58, 3)):
            d.ellipse([x - r, y - r, x + r, y + r], outline=PEARL, width=w)
        d.ellipse([x - 30, y - 30, x + 30, y + 30], fill=(14, 16, 26),
                  outline=PEARL, width=6)
        loc = c.cmax + dt.timedelta(hours=tz)
        d.text((x, y - 118), f"TOTALITY {loc:%H:%M:%S}  \u00b7  {c.duration_s:.0f} s",
               font=_font(42, True), anchor="ms", fill=(255, 120, 145),
               stroke_width=6, stroke_fill=(20, 20, 30))

    # verdict panel
    hz = hor_at(c.sun_az_at_max)
    if not np.isfinite(hz):
        raise ValueError(
            f"the horizon fan ({hor_az[0]:.1f} to {hor_az[-1]:.1f} deg) does "
            f"not cover the Sun's azimuth at maximum ({c.sun_az_at_max:.1f} "
            f"deg), so there is no clearance to report")
    clear = c.sun_alt_at_max - hz
    d.rectangle([0, 0, W, 150], fill=UMBRA)
    d.text((26, 18), title or f"{lat:.5f}, {lon:.5f}   {when.isoformat()}",
           font=_font(42, True), fill=PEARL)
    verdict = ("VISIBLE, comfortable" if clear > 5 else
               "VISIBLE, tight" if clear > 0 else "BLOCKED by terrain")
    d.text((26, 76), subtitle or
           f"{c.kind} \u00b7 Sun {c.sun_alt_at_max:.1f}\u00b0 alt, "
           f"{c.sun_az_at_max:.1f}\u00b0 az \u00b7 skyline {hz:+.1f}\u00b0 \u00b7 "
           f"clearance {clear:+.1f}\u00b0 \u00b7 {verdict}",
           font=_font(32), fill=(150, 255, 190) if clear > 0 else (255, 130, 130))
    im.save(out_path, quality=92)
    return {"clearance_deg": round(clear, 2), "horizon_deg": round(hz, 2),
            "verdict": verdict, "circumstances": c}
