"""Solar eclipse circumstances for a point on the Earth.

Uses pyephem, which needs no downloaded ephemeris file. Positions are apparent
and topocentric, so lunar parallax is included -- essential for eclipses.

Accuracy note: pyephem uses truncated VSOP87 / ELP2000. Contact times are good
to a few seconds, which is far below the tolerance needed for site selection.
For published timings, cross-check against a national observatory.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, asdict
from datetime import datetime, timedelta

import ephem

__all__ = ["Circumstances", "geometry", "obscuration", "circumstances"]


@dataclass
class Circumstances:
    lat: float
    lon: float
    elev_m: float
    kind: str  # 'total' | 'annular' | 'partial' | 'none'
    c1: datetime | None = None  # first contact, partial begins
    c2: datetime | None = None  # second contact, totality begins
    cmax: datetime | None = None  # maximum eclipse
    c3: datetime | None = None  # third contact, totality ends
    c4: datetime | None = None  # fourth contact, partial ends
    duration_s: float = 0.0
    max_obscuration: float = 0.0
    sun_alt_at_max: float = 0.0
    sun_az_at_max: float = 0.0

    def as_dict(self):
        d = asdict(self)
        for k, v in d.items():
            if isinstance(v, datetime):
                d[k] = v.isoformat()
        return d


def _observer(lat, lon, elev_m):
    o = ephem.Observer()
    o.lat = str(lat)
    o.lon = str(lon)
    o.elevation = elev_m
    o.pressure = 0  # geometric positions; refraction is applied separately
    return o


def geometry(lat, lon, elev_m, when_utc, refracted=True):
    """Return (sun_alt, sun_az, separation, r_sun, r_moon), all degrees.

    Altitude and azimuth are refracted (what you would see). Separation and
    radii are unrefracted apparent topocentric values, which is what contact
    timing needs.

    refracted=False skips the second solar solution and returns the geometric
    altitude instead. Contact search only reads the separation and the radii,
    and it calls this function some hundreds of times per site, so the saving
    is worth having. Never use it for anything you intend to display.
    """
    o = _observer(lat, lon, elev_m)
    o.date = ephem.Date(when_utc)
    sun, moon = ephem.Sun(), ephem.Moon()

    sun.compute(o)
    moon.compute(o)
    sep = math.degrees(ephem.separation((sun.ra, sun.dec), (moon.ra, moon.dec)))
    r_sun = math.degrees(sun.radius)
    r_moon = math.degrees(moon.radius)
    if not refracted:
        return math.degrees(sun.alt), math.degrees(sun.az), sep, r_sun, r_moon

    o.pressure = 1013  # now ask for the refracted altitude
    o.date = ephem.Date(when_utc)
    sun.compute(o)
    return math.degrees(sun.alt), math.degrees(sun.az), sep, r_sun, r_moon


def obscuration(sep, r_sun, r_moon):
    """Fraction of the solar disc area covered by the Moon. All degrees."""
    if sep >= r_sun + r_moon:
        return 0.0
    if sep <= abs(r_moon - r_sun):
        return 1.0 if r_moon >= r_sun else (r_moon / r_sun) ** 2
    d, r, R = sep, r_sun, r_moon
    a1 = math.acos(max(-1, min(1, (d * d + r * r - R * R) / (2 * d * r))))
    a2 = math.acos(max(-1, min(1, (d * d + R * R - r * r) / (2 * d * R))))
    area = (r * r * (a1 - math.sin(2 * a1) / 2)
            + R * R * (a2 - math.sin(2 * a2) / 2))
    return area / (math.pi * r * r)


def _root(f, t_early, t_late, tol_s=0.05):
    """Bisect f on an ordered bracket [t_early, t_late] containing a sign change.

    Returns None if the bracket does not actually straddle a root.
    """
    if t_late < t_early:
        t_early, t_late = t_late, t_early
    f0, f1 = f(t_early), f(t_late)
    if (f0 < 0) == (f1 < 0):
        return None
    while (t_late - t_early).total_seconds() > tol_s:
        tm = t_early + (t_late - t_early) / 2
        fm = f(tm)
        if (f0 < 0) == (fm < 0):
            t_early, f0 = tm, fm
        else:
            t_late = tm
    return t_early + (t_late - t_early) / 2


def circumstances(lat, lon, elev_m=0.0, date=None, search_hours=4):
    """Find eclipse contacts for one location on one UTC day.

    date: a datetime.date. The search runs over that whole UTC day.
    Returns a Circumstances object; kind == 'none' if nothing happens.
    """
    if date is None:
        raise ValueError("date is required")
    day = datetime(date.year, date.month, date.day)

    def sep(t):
        return geometry(lat, lon, elev_m, t, refracted=False)[2]

    def sep_minus_sum(t):
        _, _, s, rs, rm = geometry(lat, lon, elev_m, t, refracted=False)
        return s - (rs + rm)

    def sep_minus_diff(t):
        _, _, s, rs, rm = geometry(lat, lon, elev_m, t, refracted=False)
        return s - abs(rm - rs)

    # coarse scan for the minimum separation
    step = timedelta(minutes=5)
    t, best = day, None
    end = day + timedelta(hours=24)
    while t <= end:
        s = sep(t)
        if best is None or s < best[1]:
            best = (t, s)
        t += step
    t_near = best[0]

    # Refine the moment of maximum by ternary search. Each round keeps two
    # thirds of the bracket, so run to a tolerance rather than a fixed count:
    # 10 minutes down to 0.05 s takes about 24 rounds, not 60.
    lo, hi = t_near - step, t_near + step
    while (hi - lo).total_seconds() > 0.05:
        a = lo + (hi - lo) / 3
        b = hi - (hi - lo) / 3
        if sep(a) < sep(b):
            hi = b
        else:
            lo = a
    cmax = lo + (hi - lo) / 2
    alt, az, sep, rs, rm = geometry(lat, lon, elev_m, cmax)

    res = Circumstances(lat=lat, lon=lon, elev_m=elev_m, kind="none",
                        sun_alt_at_max=alt, sun_az_at_max=az)
    if sep >= rs + rm:
        return res

    res.cmax = cmax
    res.max_obscuration = obscuration(sep, rs, rm)
    res.kind = "total" if (sep <= abs(rm - rs) and rm >= rs) else (
        "annular" if sep <= abs(rm - rs) else "partial")

    span = timedelta(hours=search_hours)
    res.c1 = _root(sep_minus_sum, cmax - span, cmax)
    res.c4 = _root(sep_minus_sum, cmax, cmax + span)
    if res.kind in ("total", "annular"):
        res.c2 = _root(sep_minus_diff, cmax - span, cmax)
        res.c3 = _root(sep_minus_diff, cmax, cmax + span)
        if res.c2 and res.c3:
            res.duration_s = (res.c3 - res.c2).total_seconds()
    return res
