"""Horizon angles from a digital elevation model.

Two engines:

  horizon_point   exact ray march with bilinear sampling. Accurate in the near
                  field, which is where obstructions usually decide a site.
                  Use for confirming a shortlist on a fine DEM.

  horizon_grid    shift-and-max over the whole raster. Fast, but the ray is
                  quantised to whole pixels, so near-field azimuths are
                  approximate. Use for regional search on a coarse DEM.

Both apply curvature and standard atmospheric refraction via an effective
Earth radius, R_eff = R / (1 - k), with k = 0.13 by default. Both take
canopy_m in METRES, added to sampled ground heights, never to an angle.

Neither engine reports an open view where it has no data. horizon_grid masks
such cells with NaN; horizon_point reports how far the ray actually had
terrain, so the caller can see that the answer is only good to that range.

Two coverage helpers answer the prior question, which is what radius a site can
honestly be judged on at all:

  edge_reach      free. Four subtractions on the raster shape, no data read.
                  Proves a radius does not fit inside the raster.
  coverage_point  exact. Marches the ray asking only whether there is ground,
                  so it also catches nodata inside the bounds.

The DEM must be in a projected CRS with metres for both axes.
"""
from __future__ import annotations

import math
from typing import NamedTuple

import numpy as np

__all__ = ["horizon_point", "horizon_grid", "geometric_ranges", "Horizon",
           "coverage_point", "Coverage", "edge_reach", "edge_reach_point",
           "R_EARTH"]

R_EARTH = 6371000.0


class Horizon(NamedTuple):
    """One azimuth's answer from horizon_point.

    alt_deg    horizon altitude, degrees
    range_m    distance at which that horizon forms
    reach_m    farthest distance along the ray that had terrain data. Less
               than max_range_m means the ray left the raster (or crossed
               nodata) and anything beyond reach_m is unknown, not empty.
    """
    alt_deg: float
    range_m: float
    reach_m: float


class Coverage(NamedTuple):
    """One azimuth's answer from coverage_point.

    contiguous_m  distance out to which every sample had ground. Beyond the
                  first gap the terrain is unknown however much data lies
                  further out, so this is the radius the ray can be trusted to.
    reach_m       farthest sample that had ground, the same loose measure
                  horizon_point reports. Exceeds contiguous_m when the ray
                  crosses a hole and finds ground again on the far side.
    """
    contiguous_m: float
    reach_m: float


def _r_eff(k):
    return R_EARTH / (1.0 - k)


def edge_reach_point(shape, cell_m, row, col):
    """Largest radius around one cell that stays inside the raster, metres.

    Four subtractions and a min, with no data read at all, so it costs nothing
    to call on every candidate site. The nearest edge is what binds: corners
    have more room and cannot decide a circle.

    This is a necessary condition and not a sufficient one. It proves a radius
    does not fit inside the raster, but says nothing about nodata within the
    bounds. coverage_point is what settles that.
    """
    ny, nx = shape
    return float(min(row, col, (ny - 1) - row, (nx - 1) - col)) * cell_m


def edge_reach(shape, cell_m):
    """edge_reach_point for every cell at once, as a 2D array of metres."""
    ny, nx = shape
    rows = np.arange(ny, dtype=np.float32)[:, None]
    cols = np.arange(nx, dtype=np.float32)[None, :]
    return np.float32(cell_m) * np.minimum(
        np.minimum(rows, np.float32(ny - 1) - rows),
        np.minimum(cols, np.float32(nx - 1) - cols))


def geometric_ranges(cell_m, max_range_m, growth=1.02):
    """Sample distances that start at one cell and grow geometrically."""
    out, r = [], float(cell_m)
    while r <= max_range_m:
        out.append(r)
        r *= growth
    return out


def _span(n, d):
    """Destination and source slices for out[i] = src[i + d] along one axis.

    Bounds are clamped to [0, n]; Python would otherwise treat a negative stop
    as an index from the end and silently wrap.
    """
    d_lo = min(max(-d, 0), n)
    d_hi = min(max(n - d, 0), n)
    s_lo = min(max(d, 0), n)
    s_hi = min(max(n + d, 0), n)
    return slice(d_lo, d_hi), slice(s_lo, s_hi)


def _shift(z, dy, dx):
    """Array of z[i+dy, j+dx], with NaN where that falls outside."""
    ny, nx = z.shape
    out = np.full(z.shape, np.nan, dtype=z.dtype)
    yd, ys = _span(ny, dy)
    xd, xs = _span(nx, dx)
    if yd.stop > yd.start and xd.stop > xd.start:
        out[yd, xd] = z[ys, xs]
    return out


def horizon_grid(z, cell_m, azimuths_deg, max_range_m=40000.0,
                 eye_h=1.6, k=0.13, growth=1.02, canopy_m=0.0, progress=None):
    """Horizon altitude in degrees for every cell, for each azimuth.

    z            2D elevation array, metres, north-up (row 0 is the north edge)
    cell_m       pixel size in metres
    azimuths_deg iterable of compass bearings, degrees true, 0 = north
    canopy_m     metres of cover added to sampled ground, as in horizon_point.
                 Added before the max, so it can correctly move the horizon
                 onto a nearer ridge. The observer stands on bare ground.
    returns      dict {azimuth: 2D float32 array of horizon altitude in degrees}
    """
    z = np.asarray(z, dtype=np.float32)
    eye = z + np.float32(eye_h)
    canopy = np.float32(canopy_m)
    reff = _r_eff(k)
    ranges = geometric_ranges(cell_m, max_range_m, growth)
    result = {}

    for ai, az in enumerate(azimuths_deg):
        a = math.radians(az)
        sx, sy = math.sin(a), -math.cos(a)  # +x east, +y south (row index)
        best = np.full(z.shape, -90.0, dtype=np.float32)
        seen = set()
        far = (0, 0)
        for r in ranges:
            dx = int(round(r * sx / cell_m))
            dy = int(round(r * sy / cell_m))
            if (dy, dx) == (0, 0) or (dy, dx) in seen:
                continue
            seen.add((dy, dx))
            far = (dy, dx)
            # true distance of the quantised sample, not the requested one
            d = math.hypot(dx * cell_m, dy * cell_m)
            zs = _shift(z, dy, dx)
            ang = np.degrees(np.arctan(
                (zs + canopy - eye - d * d / (2 * reff)) / d))
            np.fmax(best, ang, out=best)
        # A cell whose ray left the raster before max_range has an unknown
        # horizon: terrain beyond the edge could be anything. Mark it NaN
        # rather than reporting the sentinel as a wide-open view.
        reached = np.isfinite(_shift(z, *far)) if far != (0, 0) else np.zeros(
            z.shape, bool)
        # An observer standing on nodata has no elevation of its own, so every
        # angle along its ray is NaN and np.fmax quietly preserves the -90
        # sentinel. That is the same lie as an unreachable cell, and it ranked
        # a nodata strip at +100 degrees of clearance on the real Asturias DEM.
        reached &= np.isfinite(z)
        best[~reached] = np.nan
        result[az] = best
        if progress:
            progress(ai + 1, len(azimuths_deg))
    return result


def _bilinear(z, row, col):
    ny, nx = z.shape
    if row < 0 or col < 0 or row > ny - 1.001 or col > nx - 1.001:
        return np.nan
    r0, c0 = int(row), int(col)
    fr, fc = row - r0, col - c0
    return (z[r0, c0] * (1 - fr) * (1 - fc) + z[r0, c0 + 1] * (1 - fr) * fc
            + z[r0 + 1, c0] * fr * (1 - fc) + z[r0 + 1, c0 + 1] * fr * fc)


def horizon_point(z, cell_m, row, col, azimuths_deg, max_range_m=40000.0,
                  eye_h=1.6, k=0.13, step_m=None, canopy_m=0.0):
    """Horizon altitude in degrees at one cell, marching the exact ray.

    row, col may be fractional. canopy_m is added to sampled ground heights,
    a crude allowance for tree cover; leave at 0 for bare terrain.

    Returns dict {azimuth: Horizon(alt_deg, range_m, reach_m)}. Check reach_m
    against max_range_m before trusting a low horizon: a ray that ran out of
    raster at 3 km has not established that the next 37 km are clear, it has
    only failed to look. This is the point-engine counterpart of the NaN
    masking horizon_grid does for the same reason.
    """
    z = np.asarray(z, dtype=np.float32)
    step = step_m or cell_m / 2.0
    reff = _r_eff(k)
    z0 = _bilinear(z, row, col)
    if not np.isfinite(z0):
        raise ValueError("observer falls outside the DEM")
    eye = float(z0) + eye_h

    out = {}
    for az in azimuths_deg:
        a = math.radians(az)
        sx, sy = math.sin(a), -math.cos(a)
        best, best_r, reach = -90.0, 0.0, 0.0
        d = step
        while d <= max_range_m:
            zs = _bilinear(z, row + sy * d / cell_m, col + sx * d / cell_m)
            if np.isfinite(zs):
                reach = d
                ang = math.degrees(math.atan(
                    (float(zs) + canopy_m - eye - d * d / (2 * reff)) / d))
                if ang > best:
                    best, best_r = ang, d
            d += step
        out[az] = Horizon(best, best_r, reach)
    return out


def coverage_point(z, cell_m, row, col, azimuths_deg, max_range_m=40000.0,
                   step_m=None):
    """How far the DEM actually has ground along each ray from one cell.

    The cheap counterpart to horizon_point. It asks only whether there is data,
    never how high it is, so it needs no curvature, no eye height and no fine
    step. Run it before paying for a horizon to find out what radius the site
    can honestly be judged on.

    step_m defaults to one cell rather than the half cell horizon_point uses.
    The edge of valid data is a pixel-scale feature and half steps buy nothing;
    a hole smaller than a cell makes the bilinear sample NaN in any case.

    Returns dict {azimuth: Coverage(contiguous_m, reach_m)}. Compare
    contiguous_m against the range you mean to ask for. A site whose rays run
    out at 15 km has not been shown to be better than one surveyed to 40 km,
    it has only been looked at less.
    """
    z = np.asarray(z, dtype=np.float32)
    step = step_m or cell_m
    if not np.isfinite(_bilinear(z, row, col)):
        raise ValueError("observer falls outside the DEM")

    out = {}
    for az in azimuths_deg:
        a = math.radians(az)
        sx, sy = math.sin(a), -math.cos(a)
        contiguous, reach, broken = 0.0, 0.0, False
        d = step
        while d <= max_range_m:
            if np.isfinite(_bilinear(z, row + sy * d / cell_m,
                                     col + sx * d / cell_m)):
                reach = d
                if not broken:
                    contiguous = d
            else:
                broken = True
            d += step
        out[az] = Coverage(contiguous, reach)
    return out
