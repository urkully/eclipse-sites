"""Turn a folder of downloaded tiles into rasters this package can use.

Point it at whatever came out of the download. Zips, mixed products, mixed UTM
zones, ESRI ASCII grids with no projection attached. It sorts by product, works
out the CRS, merges each group and writes one clean GeoTIFF per product into a
directory you name.

Two things it is careful about, because both silently produce wrong answers:

  * IGN ASCII grids usually carry no CRS. The zone is inferred from the file
    name (HU29, H30 and similar) and, failing that, you are asked rather than
    guessed at.
  * Normalised surface models hold heights above ground, not altitudes, and are
    resampled with MAX rather than bilinear. For a horizon you want the tallest
    thing in the cell; averaging a treeline is how you lose it.
"""
from __future__ import annotations

import json
import os
import re
import shutil
import zipfile
from dataclasses import dataclass, asdict, field

__all__ = ["scan", "prepare", "PRODUCTS"]

RASTER_EXT = {".tif", ".tiff", ".asc", ".txt", ".vrt"}

# name fragment -> (output name, nominal metres, is_normalised_height)
PRODUCTS = [
    (r"MDT02|DTM02", "dtm02", 2.0, False),
    (r"MDT05|DTM05", "dtm05", 5.0, False),
    (r"MDT25|DTM25", "dtm25", 25.0, False),
    (r"MDT200|DTM200", "dtm200", 200.0, False),
    (r"MDT50|DTM50", "dtm50cm", 0.5, False),
    (r"MDSNV|NDSM2?,?5?_?V|VEGET", "ndsm_vegetation", 2.5, True),
    (r"MDSNE|NDSM2?,?5?_?E|EDIFIC|BUILDING", "ndsm_building", 2.5, True),
    (r"MDS02|DSM02", "dsm02", 2.0, False),
    (r"MDS05|DSM05", "dsm05", 5.0, False),
]


@dataclass
class Tile:
    path: str
    product: str
    crs: str | None
    res: tuple | None
    bounds: tuple | None
    note: str = ""


@dataclass
class Report:
    out_dir: str
    target_crs: str = ""
    written: list = field(default_factory=list)
    skipped: list = field(default_factory=list)
    warnings: list = field(default_factory=list)


def _unzip(indir, workdir):
    """Extract any zips into workdir and return the extra search root."""
    found = False
    for root, _, files in os.walk(indir):
        for f in files:
            if f.lower().endswith(".zip"):
                try:
                    with zipfile.ZipFile(os.path.join(root, f)) as z:
                        z.extractall(workdir)
                    found = True
                except zipfile.BadZipFile:
                    pass
    return workdir if found else None


def _classify(name):
    up = name.upper().replace("-", "_")
    for pat, out, res, norm in PRODUCTS:
        if re.search(pat, up):
            return out, res, norm
    return None, None, False


def _zone_from_name(name):
    """IGN encodes the UTM zone as HU29, H30, HUSO29 and similar."""
    m = re.search(r"H(?:U|USO)?[_-]?(29|30|31|28)", name.upper())
    if m:
        z = int(m.group(1))
        return f"EPSG:258{z}"          # ETRS89 / UTM zone n N
    return None


def scan(indir, src_crs=None):
    """Inspect every raster under indir. Returns a list of Tile."""
    import rasterio
    tiles = []
    for root, _, files in os.walk(indir):
        for f in sorted(files):
            if os.path.splitext(f)[1].lower() not in RASTER_EXT:
                continue
            p = os.path.join(root, f)
            product, _, _ = _classify(f)
            crs = res = bounds = None
            note = ""
            try:
                with rasterio.open(p) as ds:
                    crs = str(ds.crs) if ds.crs else None
                    res = (abs(ds.transform.a), abs(ds.transform.e))
                    bounds = tuple(ds.bounds)
            except Exception as e:
                note = f"unreadable: {e}"
            if crs is None:
                crs = src_crs or _zone_from_name(f)
                if crs:
                    note = (note + " " if note else "") + f"CRS assumed {crs}"
            if product is None:
                product = f"res{int(res[0])}m" if res else "unknown"
                note = (note + " " if note else "") + "product guessed from pixel size"
            tiles.append(Tile(p, product, crs, res, bounds, note.strip()))
    return tiles


def _target_crs(tiles):
    """Most common source CRS wins; ties break toward the lower zone number."""
    counts = {}
    for t in tiles:
        if t.crs:
            counts[t.crs] = counts.get(t.crs, 0) + 1
    if not counts:
        return None
    return sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]


def prepare(indir, outdir, src_crs=None, target_crs=None, only=None,
            res=None, max_cells=4e8, progress=None):
    """Merge and reproject everything under indir into outdir.

    max_cells caps the mosaic. merge() builds the whole thing in memory, so a
    scattered download that spans a country is an out-of-memory kill rather
    than an error message. The estimate is made from the warped bounds, before
    anything is allocated.
    """
    import rasterio
    from rasterio.merge import merge
    from rasterio.vrt import WarpedVRT
    from rasterio.enums import Resampling
    from rasterio.warp import calculate_default_transform

    os.makedirs(outdir, exist_ok=True)
    work = os.path.join(outdir, "_unzipped")
    extra = _unzip(indir, work)

    tiles = scan(indir, src_crs)
    if extra:
        tiles += scan(extra, src_crs)
    rep = Report(out_dir=outdir)

    usable = [t for t in tiles if t.bounds and t.crs]
    keep = {id(t) for t in usable}          # identity, not O(n) field compare
    for t in tiles:
        if id(t) not in keep:
            rep.skipped.append({"path": t.path, "why": t.note or "no CRS"})
    if not usable:
        rep.warnings.append(
            "nothing usable. If these are ESRI ASCII grids with no projection, "
            "pass --src-crs, e.g. EPSG:25829 for UTM zone 29.")
        shutil.rmtree(work, ignore_errors=True)
        return rep

    tgt = target_crs or _target_crs(usable)
    rep.target_crs = tgt

    groups = {}
    for t in usable:
        groups.setdefault(t.product, []).append(t)

    for i, (product, members) in enumerate(sorted(groups.items())):
        if only and product not in only:
            continue
        # Classify on the basename, as scan() did. Matching the full path lets
        # a parent directory named for another product override the file's own
        # name, which would quietly turn MAX resampling off on a vegetation
        # model while still labelling the output ndsm_vegetation.
        _, nominal, normalised = _classify(os.path.basename(members[0].path))
        out_res = res or nominal or members[0].res[0]
        # MAX is a warp-only resampling: it aggregates a treeline correctly
        # instead of averaging it away. It is rejected by merge(), which reads.
        # So each tile is warped to the target grid first, then merged with no
        # further resampling.
        method = Resampling.max if normalised else Resampling.bilinear
        srcs = []
        try:
            for t in members:
                ds = rasterio.open(t.path)
                tr, w, h = calculate_default_transform(
                    ds.crs, tgt, ds.width, ds.height, *ds.bounds,
                    resolution=out_res)
                srcs.append(WarpedVRT(ds, crs=tgt, transform=tr, width=w,
                                      height=h, resampling=method))

            # Size the mosaic from the warped bounds before merge() allocates
            # it. This check used to run on the finished array, which is too
            # late to be useful when the array is what killed the process.
            xs = [b for s in srcs for b in (s.bounds[0], s.bounds[2])]
            ys = [b for s in srcs for b in (s.bounds[1], s.bounds[3])]
            est = ((max(xs) - min(xs)) / out_res) * ((max(ys) - min(ys)) / out_res)
            if est > max_cells:
                rep.warnings.append(
                    f"{product}: mosaic would be about "
                    f"{(max(xs)-min(xs))/out_res:.0f}x"
                    f"{(max(ys)-min(ys))/out_res:.0f} cells "
                    f"({est/1e6:.0f} million), over the {max_cells/1e6:.0f} "
                    "million limit. Tiles are probably scattered rather than "
                    "adjacent, so most of it would be empty. Download a "
                    "contiguous block, or raise --max-cells.")
                continue

            # Gaps between tiles MUST be declared nodata. If they default to
            # 0 the horizon engine reads them as ground at sea level, sees a
            # cliff-edge drop and reports a clearer view than exists.
            nod = srcs[0].nodata
            if nod is None:
                nod = -9999.0
            arr, transform = merge(srcs, res=(out_res, out_res),
                                   resampling=Resampling.nearest, nodata=nod)
            path = os.path.join(outdir, f"{product}.tif")
            with rasterio.open(
                    path, "w", driver="GTiff", height=arr.shape[1],
                    width=arr.shape[2], count=1, dtype="float32", crs=tgt,
                    transform=transform, nodata=nod,
                    compress="deflate", tiled=True) as dst:
                dst.write(arr[0].astype("float32"), 1)
            rep.written.append({
                "product": product, "path": path, "tiles": len(members),
                "res_m": out_res, "crs": tgt,
                "resampling": "max" if normalised else "bilinear",
                "size": [int(arr.shape[2]), int(arr.shape[1])],
                "normalised_heights": normalised})
            if normalised:
                rep.warnings.append(
                    f"{product}: heights are above ground, not altitudes. "
                    "Do not pass it as --dem.")
        except Exception as e:
            rep.warnings.append(f"{product}: failed, {e}")
        finally:
            for s in srcs:
                try:
                    s.close()
                except Exception:
                    pass
        if progress:
            progress(i + 1, len(groups))

    shutil.rmtree(work, ignore_errors=True)
    with open(os.path.join(outdir, "manifest.json"), "w") as fh:
        json.dump(asdict(rep), fh, indent=2)
    return rep
