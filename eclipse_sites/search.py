"""Regional search: find every place with a clear view of totality.

Strategy is staged, because a fine DEM over a whole region is expensive and
unnecessary for the first pass:

  1. sweep a coarse DEM (25 m is plenty) with the shift-and-max engine
  2. rank the cells by clearance at maximum eclipse
  3. re-check the shortlist on a fine DEM (5 m) with the exact ray march

Eclipse circumstances vary across a large region, so the Sun's azimuth and
altitude are computed on a coarse lat/lon lattice and interpolated per cell
rather than assumed constant.
"""
from __future__ import annotations

import csv
import math
from datetime import date as _date

import numpy as np

from .ephemeris import circumstances
from .horizon import horizon_grid, horizon_point

__all__ = ["read_dem", "sun_fields", "sweep", "confirm", "write_geotiff"]


def read_dem(path, band=1):
    """Return (z, transform, crs, cell_m). Requires a projected CRS in metres."""
    import rasterio
    with rasterio.open(path) as src:
        z = src.read(band, masked=True).filled(np.nan).astype(np.float32)
        tr, crs = src.transform, src.crs
        if crs is None or crs.is_geographic:
            raise ValueError(
                "DEM must be in a projected CRS in metres. Reproject first, e.g.\n"
                "  gdalwarp -t_srs EPSG:25830 -r bilinear in.tif out.tif")
        cx, cy = abs(tr.a), abs(tr.e)
        if abs(cx - cy) > 1e-6:
            raise ValueError(f"non-square pixels ({cx} x {cy}); resample first")
    return z, tr, crs, cx


def _lattice(z, tr, crs, n=3):
    """n x n lattice of (row, col, lat, lon) spanning the raster."""
    from pyproj import Transformer
    to_wgs = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
    ny, nx = z.shape
    rows = np.linspace(0, ny - 1, n)
    cols = np.linspace(0, nx - 1, n)
    pts = []
    for r in rows:
        for c in cols:
            x, y = tr * (c + 0.5, r + 0.5)
            lon, lat = to_wgs.transform(x, y)
            pts.append((r, c, lat, lon))
    return rows, cols, pts


def sun_fields(z, tr, crs, when, lattice_n=3, elev_m=0.0, progress=None):
    """Per-cell Sun altitude and azimuth at maximum eclipse.

    when: a datetime.date. Returns (alt2d, az2d, sample_circumstances).
    """
    rows, cols, pts = _lattice(z, tr, crs, lattice_n)
    alt = np.zeros((lattice_n, lattice_n))
    az = np.zeros((lattice_n, lattice_n))
    sample = None
    for i, (r, c, lat, lon) in enumerate(pts):
        cc = circumstances(lat, lon, elev_m, when)
        alt.flat[i] = cc.sun_alt_at_max
        az.flat[i] = cc.sun_az_at_max
        if sample is None or cc.kind == "total":
            sample = cc
        if progress:
            progress(i + 1, len(pts))

    ny, nx = z.shape
    rr = np.arange(ny, dtype=np.float32)
    cc_ = np.arange(nx, dtype=np.float32)

    def bilerp(grid):
        a = np.empty((ny, len(cols)), np.float32)
        for j in range(len(cols)):
            a[:, j] = np.interp(rr, rows, grid[:, j])
        out = np.empty((ny, nx), np.float32)
        for i in range(ny):
            out[i] = np.interp(cc_, cols, a[i])
        return out

    return bilerp(alt), bilerp(az), sample


def _az_bins(az2d, az_step):
    """Azimuth bins spanning a field of Sun azimuths.

    Returns (az_linearised, lo, bin_azimuths). Two things to get right:

    * A region straddling due north holds azimuths near both 0 and 360.
      Binning those raw sweeps the whole compass and interpolates the horizon
      the long way round, so unwrap onto one branch first.
    * There must always be at least two bins. A field that lands on a single
      bin edge otherwise indexes one past the end of the horizon stack.
    """
    az_lin = np.where(np.nanmax(az2d) - np.nanmin(az2d) > 180.0,
                      np.where(az2d < 180.0, az2d + 360.0, az2d), az2d)
    lo = math.floor(np.nanmin(az_lin) / az_step) * az_step
    hi = math.ceil(np.nanmax(az_lin) / az_step) * az_step
    azs = list(np.arange(lo, hi + az_step / 2, az_step))
    if len(azs) < 2:
        azs.append(lo + az_step)
    return az_lin, lo, azs


def sweep(dem_path, when, eye_h=1.6, max_range_m=40000.0, az_step=2.5,
          canopy_m=0.0, lattice_n=3, progress=None):
    """Coarse regional pass. Returns (clearance2d, z, transform, crs, cell, circ).

    canopy_m is METRES, and is handed to the horizon engine to be added to
    ground heights, exactly as confirm() does on the fine DEM. It used to be
    added to the horizon ANGLE here, which made the coarse clearance column
    disagree with the fine one by tens of degrees for the same flag.
    """
    z, tr, crs, cell = read_dem(dem_path)
    alt2d, az2d, circ = sun_fields(z, tr, crs, when, lattice_n, progress=progress)

    az_lin, lo, azs = _az_bins(az2d, az_step)

    hz = horizon_grid(z, cell, azs, max_range_m=max_range_m, eye_h=eye_h,
                      canopy_m=canopy_m, progress=progress)
    stack = np.stack([hz[a] for a in azs])          # (n_az, ny, nx)
    idx = np.clip((az_lin - lo) / az_step, 0, len(azs) - 1.0001)
    i0 = idx.astype(np.int32)
    f = (idx - i0).astype(np.float32)
    horizon = (np.take_along_axis(stack, i0[None], 0)[0] * (1 - f)
               + np.take_along_axis(stack, i0[None] + 1, 0)[0] * f)
    return (alt2d - horizon).astype(np.float32), z, tr, crs, cell, circ


def _best_first(flat, want):
    """Indices of flat in descending order, cheapest way that yields `want`.

    A full argsort of a regional raster sorts tens of millions of cells to use
    the first few dozen. Take a generous slice with argpartition instead, and
    only fall back to the full sort if spacing rejects the whole slice.
    """
    k = min(flat.size, max(1024, want))
    part = np.argpartition(flat, -k)[-k:]
    yield from part[np.argsort(flat[part])[::-1]]
    if k < flat.size:
        rest = np.setdiff1d(np.argsort(flat)[::-1], part, assume_unique=False)
        yield from rest[np.argsort(flat[rest])[::-1]]


def top_sites(clearance, z, tr, crs, n=25, min_sep_cells=20):
    """Rank cells by clearance, keeping results spread out."""
    from pyproj import Transformer
    to_wgs = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
    flat = np.where(np.isfinite(clearance), clearance, -np.inf).ravel()
    ny, nx = clearance.shape
    # enough candidates that the spacing filter can still find n of them
    order = _best_first(flat, n * max(1, min_sep_cells) ** 2)
    picked, out = [], []
    for k in order:
        if len(out) >= n:
            break
        if not np.isfinite(flat[k]):
            break
        r, c = divmod(int(k), nx)
        if any(abs(r - pr) < min_sep_cells and abs(c - pc) < min_sep_cells
               for pr, pc in picked):
            continue
        picked.append((r, c))
        x, y = tr * (c + 0.5, r + 0.5)
        lon, lat = to_wgs.transform(x, y)
        out.append({"rank": len(out) + 1, "lat": round(lat, 6),
                    "lon": round(lon, 6), "elev_m": round(float(z[r, c]), 1),
                    "clearance_deg": round(float(clearance[r, c]), 2),
                    "row": r, "col": c})
    return out


def confirm(fine_dem_path, sites, when, eye_h=1.6, canopy_m=0.0,
            max_range_m=40000.0):
    """Re-check a shortlist on a fine DEM with the exact ray march.

    Returns (confirmed, dropped). A site the fine DEM does not cover cannot be
    confirmed, and used to disappear from the table without a word, leaving
    gaps in the rank column that read like a display fault.
    """
    from pyproj import Transformer
    z, tr, crs, cell = read_dem(fine_dem_path)
    to_dem = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    inv = ~tr
    out, dropped = [], []
    for s in sites:
        cc = circumstances(s["lat"], s["lon"], s.get("elev_m", 0.0), when)
        x, y = to_dem.transform(s["lon"], s["lat"])
        col, row = inv * (x, y)
        try:
            h = horizon_point(z, cell, row - 0.5, col - 0.5, [cc.sun_az_at_max],
                              max_range_m=max_range_m, eye_h=eye_h,
                              canopy_m=canopy_m)
        except ValueError:
            dropped.append(dict(s, why="outside the fine DEM"))
            continue
        hz, hr, reach = h[cc.sun_az_at_max]
        r = dict(s)
        r["coarse_rank"] = r.pop("rank", None)
        r.update({"fine_clearance_deg": round(cc.sun_alt_at_max - hz, 2),
                  "horizon_deg": round(hz, 2), "horizon_range_m": round(hr),
                  "reach_m": round(reach),
                  "ray_truncated": reach < max_range_m - cell,
                  "sun_alt_deg": round(cc.sun_alt_at_max, 2),
                  "sun_az_deg": round(cc.sun_az_at_max, 2),
                  "totality_s": round(cc.duration_s, 1)})
        out.append(r)
    out.sort(key=lambda d: -d["fine_clearance_deg"])
    for i, r in enumerate(out, 1):
        r["rank"] = i
    return out, dropped


def write_geotiff(path, arr, tr, crs, nodata=-9999.0):
    import rasterio
    a = np.where(np.isfinite(arr), arr, nodata).astype(np.float32)
    with rasterio.open(path, "w", driver="GTiff", height=a.shape[0],
                       width=a.shape[1], count=1, dtype="float32",
                       crs=crs, transform=tr, nodata=nodata,
                       compress="deflate") as dst:
        dst.write(a, 1)


def write_csv(path, rows):
    if not rows:
        return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
