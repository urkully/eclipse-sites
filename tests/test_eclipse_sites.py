"""Run with:  python -m pytest tests/ -v

No network and no DEM download needed. Terrain is synthetic with analytic
answers; the ephemeris case is a real eclipse with published timings.
"""
import datetime as dt
import math

import numpy as np
import pytest

from eclipse_sites.ephemeris import circumstances, obscuration, geometry
from eclipse_sites.horizon import (horizon_point, horizon_grid, _shift,
                                   coverage_point, edge_reach,
                                   edge_reach_point)

CELL = 10.0
N = 401
C = 200
R_EFF = 6371000.0 / 0.87


# --------------------------------------------------------------- ephemeris

def test_candamo_2026_contacts():
    """Real case, published locally as C1 ~19:30, C2 ~20:26, C4 ~21:20 CEST."""
    c = circumstances(43.4524508, -6.0705374, 117.6, dt.date(2026, 8, 12))
    assert c.kind == "total"
    cest = lambda t: t + dt.timedelta(hours=2)
    assert cest(c.c1).strftime("%H:%M") == "19:31"
    assert cest(c.c2).strftime("%H:%M") == "20:26"
    assert cest(c.c4).strftime("%H:%M") == "21:20"
    assert 100 < c.duration_s < 130
    assert c.max_obscuration == pytest.approx(1.0)
    assert c.sun_alt_at_max == pytest.approx(10.46, abs=0.15)
    assert c.sun_az_at_max == pytest.approx(280.6, abs=0.3)


def test_no_eclipse_returns_none_kind():
    c = circumstances(43.45, -6.07, 100.0, dt.date(2026, 8, 20))
    assert c.kind == "none"
    assert c.c1 is None


def test_ordering_of_contacts():
    c = circumstances(43.4524508, -6.0705374, 117.6, dt.date(2026, 8, 12))
    assert c.c1 < c.c2 < c.cmax < c.c3 < c.c4


def test_obscuration_bounds():
    assert obscuration(2.0, 0.26, 0.27) == 0.0        # discs apart
    assert obscuration(0.0, 0.26, 0.27) == 1.0        # total
    mid = obscuration(0.30, 0.26, 0.27)
    assert 0.0 < mid < 1.0
    # annular: Moon smaller than Sun, centred
    assert obscuration(0.0, 0.27, 0.26) == pytest.approx((0.26 / 0.27) ** 2)


def test_geometry_matches_known_sun_position():
    alt, az, sep, rs, rm = geometry(43.4524508, -6.0705374, 117.6,
                                    dt.datetime(2026, 8, 12, 18, 27, 25))
    assert alt == pytest.approx(10.5, abs=0.2)
    assert az == pytest.approx(280.5, abs=0.3)
    assert sep < abs(rm - rs)      # inside totality


# ------------------------------------------------------------------- shift

def test_shift_matches_bruteforce_including_oversized_offsets():
    z = np.arange(20, dtype=np.float32).reshape(4, 5)
    for dy in range(-7, 8):
        for dx in range(-8, 9):
            got = _shift(z, dy, dx)
            ref = np.full_like(z, np.nan)
            for i in range(4):
                for j in range(5):
                    si, sj = i + dy, j + dx
                    if 0 <= si < 4 and 0 <= sj < 5:
                        ref[i, j] = z[si, sj]
            assert np.array_equal(np.nan_to_num(got, nan=-1),
                                  np.nan_to_num(ref, nan=-1)), (dy, dx)


# ----------------------------------------------------------------- horizon

def test_flat_plane_is_slightly_below_zero():
    z = np.zeros((N, N), np.float32)
    r = horizon_point(z, CELL, C, C, [0, 90, 180, 270], max_range_m=1500)
    for az, (ang, _, _) in r.items():
        assert -0.2 < ang < 0.0     # curvature only


def test_wall_east_matches_analytic():
    z = np.zeros((N, N), np.float32)
    z[:, C + 50:] = 100.0
    d = 500.0
    expect = math.degrees(math.atan((100 - 1.6 - d * d / (2 * R_EFF)) / d))
    ang, rng, _ = horizon_point(z, CELL, C, C, [90], max_range_m=1500)[90]
    assert ang == pytest.approx(expect, abs=1e-3)
    assert rng == pytest.approx(d, abs=CELL)


def test_azimuth_convention_zero_is_north():
    z = np.zeros((N, N), np.float32)
    z[:C - 50, :] = 100.0            # low row index = north; edge at 510 m
    r = horizon_point(z, CELL, C, C, [0, 90, 180, 270], max_range_m=1500)
    d = 510.0
    expect = math.degrees(math.atan((100 - 1.6 - d * d / (2 * R_EFF)) / d))
    assert r[0][0] == pytest.approx(expect, abs=1e-3)
    assert r[180][0] < 0
    assert r[90][0] < 0


def test_grid_agrees_with_point():
    z = np.zeros((N, N), np.float32)
    z[:, C + 50:] = 100.0
    g = horizon_grid(z, CELL, [90], max_range_m=1500)[90]
    p = horizon_point(z, CELL, C, C, [90], max_range_m=1500)[90][0]
    assert g[C, C] == pytest.approx(p, abs=1e-3)


def test_horizon_falls_as_you_retreat_from_a_cone():
    z = np.zeros((N, N), np.float32)
    yy, xx = np.mgrid[0:N, 0:N]
    d = np.hypot(yy - C, xx - (C + 50)) * CELL
    z = np.maximum(0, 150 - d * 0.5).astype(np.float32)
    near = horizon_point(z, CELL, C, C, [90], max_range_m=2000)[90][0]
    far = horizon_point(z, CELL, C, C - 30, [90], max_range_m=2000)[90][0]
    assert far < near


def test_canopy_raises_the_horizon():
    z = np.zeros((N, N), np.float32)
    z[:, C + 50:] = 100.0
    bare = horizon_point(z, CELL, C, C, [90], max_range_m=1500)[90][0]
    tree = horizon_point(z, CELL, C, C, [90], max_range_m=1500,
                         canopy_m=20.0)[90][0]
    assert tree > bare + 1.0


def test_unreachable_cells_are_masked_not_wide_open():
    """A ray that leaves the raster must give NaN, never the -90 sentinel."""
    # 50 cells at 10 m = 500 m wide. Marching 200 m west needs 20 cells of
    # room, so the western strip is unknowable and the eastern side is fine.
    z = np.zeros((50, 50), np.float32)
    g = horizon_grid(z, CELL, [270], max_range_m=200)[270]
    assert np.isnan(g[:, :15]).all()     # ray exits the raster
    assert np.isfinite(g[:, 25:]).all()  # room to march
    assert np.nanmax(g) < 0.0            # flat ground, nothing above eye
    assert not (g == -90.0).any()        # sentinel must never survive


def test_observer_on_nodata_is_masked_not_wide_open():
    """A cell with no elevation of its own cannot report a horizon.

    Every angle along its ray is NaN, and np.fmax ignores NaN, so the -90
    sentinel survived and read as a wide open view. On the real Asturias DEM
    this ranked a strip of nodata at +100 degrees of clearance, top of the
    whole region.
    """
    z = np.zeros((80, 80), np.float32)
    z[30:35, 30:35] = np.nan            # a hole the observer stands in
    g = horizon_grid(z, CELL, [90], max_range_m=200)[90]
    assert np.isnan(g[30:35, 30:35]).all()
    assert not (g == -90.0).any()
    assert np.nanmax(g) < 0.0           # flat ground, nothing above eye


def test_observer_outside_dem_raises():
    z = np.zeros((20, 20), np.float32)
    with pytest.raises(ValueError):
        horizon_point(z, CELL, 500, 500, [0])


def test_point_engine_reports_where_it_ran_out_of_dem():
    """A short ray must not read as a confirmed clear view to max_range.

    horizon_grid masks these cells NaN. The point engine is what produces the
    final answer, so it has to say the same thing in its own way: reach_m is
    how far it actually looked.
    """
    z = np.zeros((60, 60), np.float32)          # 600 m across at CELL = 10
    h = horizon_point(z, CELL, 30, 30, [270], max_range_m=20000)[270]
    assert h.reach_m < 400                       # ran out after ~300 m
    assert h.reach_m < 20000
    # and the grid engine refuses the same query outright
    g = horizon_grid(z, CELL, [270], max_range_m=20000)[270]
    assert np.isnan(g).all()


def test_reach_equals_max_range_when_the_dem_is_big_enough():
    z = np.zeros((N, N), np.float32)
    h = horizon_point(z, CELL, C, C, [90], max_range_m=1500)[90]
    assert h.reach_m == pytest.approx(1500, abs=CELL)


# ---------------------------------------------------------------- coverage

def test_edge_reach_is_the_distance_to_the_nearest_edge():
    """The nearest edge binds a circle. Corners have more room and cannot."""
    assert edge_reach_point((N, N), CELL, C, C) == pytest.approx(2000.0)
    assert edge_reach_point((N, N), CELL, 50, C) == pytest.approx(500.0)
    assert edge_reach_point((N, N), CELL, C, N - 1 - 50) == pytest.approx(500.0)


def test_edge_reach_grid_matches_the_point_version():
    g = edge_reach((N, N), CELL)
    for r, c in ((C, C), (0, 0), (50, 137), (N - 1, N - 1), (3, N - 4)):
        assert g[r, c] == pytest.approx(edge_reach_point((N, N), CELL, r, c))


def test_coverage_reports_the_raster_edge_when_asked_beyond_it():
    z = np.zeros((N, N), np.float32)
    c = coverage_point(z, CELL, C, C, [90], max_range_m=4000)[90]
    assert c.contiguous_m == pytest.approx(2000.0, abs=2 * CELL)


def test_coverage_stops_at_the_first_gap_not_the_last_data():
    """Ground beyond a hole was seen. The hole itself was not.

    A ray that resumes at 500 m has still not shown what sits at 350 m, so the
    trustworthy radius ends at the gap even though data continues past it.
    This is the distinction reach_m alone cannot make.
    """
    z = np.zeros((N, N), np.float32)
    z[:, C + 30:C + 35] = np.nan          # nodata band 300 to 350 m east
    c = coverage_point(z, CELL, C, C, [90], max_range_m=1500)[90]
    assert c.contiguous_m == pytest.approx(280.0, abs=2 * CELL)
    assert c.reach_m == pytest.approx(1500.0, abs=CELL)
    assert c.reach_m > c.contiguous_m


def test_coverage_reach_agrees_with_horizon_point():
    """The loose measure must mean the same thing in both engines."""
    z = np.zeros((N, N), np.float32)
    z[:, C + 30:C + 35] = np.nan
    azs = [0, 45, 90, 200, 315]
    cov = coverage_point(z, CELL, C, C, azs, max_range_m=1500, step_m=CELL)
    hor = horizon_point(z, CELL, C, C, azs, max_range_m=1500, step_m=CELL)
    for az in azs:
        assert cov[az].reach_m == pytest.approx(hor[az].reach_m, abs=1e-6)


def test_edge_reach_is_necessary_but_not_sufficient():
    """The bounds test is free, and free buys only an upper bound.

    It proves a radius does not fit inside the raster. It knows nothing about
    nodata within the bounds, so it can only ever be optimistic, which is why
    coverage_point has to exist alongside it.
    """
    z = np.zeros((N, N), np.float32)
    z[:, C + 30:C + 35] = np.nan
    assert edge_reach_point(z.shape, CELL, C, C) == pytest.approx(2000.0)
    c = coverage_point(z, CELL, C, C, [90], max_range_m=1500)[90]
    assert c.contiguous_m < 400.0


def test_coverage_refuses_a_nodata_observer():
    """Same precondition as horizon_point, and for the same reason."""
    z = np.zeros((N, N), np.float32)
    z[C, C] = np.nan
    with pytest.raises(ValueError, match="outside the DEM"):
        coverage_point(z, CELL, C, C, [90])


# ----------------------------------------------------------------- horizon

def test_grid_canopy_matches_point_canopy():
    """Both engines take canopy in METRES, and agree on what that means.

    The observer sits on a pinnacle above the local treetops, because a global
    canopy offset otherwise puts 20 m of trees in the very next cell and the
    answer becomes whichever engine happens to sample nearest. That is real
    behaviour, not a bug, but it hides the thing under test here.
    """
    z = np.zeros((N, N), np.float32)
    z[:, C + 50:] = 100.0          # wall 500 m east
    z[C, C] = 60.0                 # eye at 61.6 m, above 20 m of canopy
    p = horizon_point(z, CELL, C, C, [90], max_range_m=1500, canopy_m=20.0)[90]
    g = horizon_grid(z, CELL, [90], max_range_m=1500, canopy_m=20.0)[90]
    assert g[C, C] == pytest.approx(p.alt_deg, abs=1e-3)
    # the wall is what forms the horizon, raised by exactly its 20 m of trees
    d = 500.0
    expect = math.degrees(math.atan((120 - 61.6 - d * d / (2 * R_EFF)) / d))
    assert p.alt_deg == pytest.approx(expect, abs=1e-3)
    assert p.range_m == pytest.approx(d, abs=CELL)


# ------------------------------------------------------------------ render

def test_render_view_reads_a_fixed_fan_at_the_sun_bearing(tmp_path):
    """Render maths against a fixed horizon fan. Not evidence about the site.

    The fan is EU-DEM derived and is now known to be wrong: it stopped looking
    at 14 km and missed the 753 m massif at 18.7 km that forms the real
    horizon, so it puts the clearance 1.2 deg too high. The live answer from
    the PNOA 5 m DTM is +8.6 deg. See FINDINGS-anchor.md and invariant 13.

    It is kept as a fixture because what is under test here is render_view's
    own arithmetic: given this fan, does it interpolate to the Sun's bearing
    and report this clearance. That is still worth pinning. It was once
    described as corroborated by a photogrammetric route agreeing to 0.08 deg,
    which was never independent: fit_camera solves a camera against the fan it
    is given and cannot disagree with its own input.
    """
    from eclipse_sites.render import render_view
    az = [230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330]
    alt = [5.84, 6.68, 1.85, 1.65, 2.67, 0.64, 1.97, 3.79, 8.45, 12.86, 16.86]
    rng = [3500, 3500, 3500, 7000, 3500, 3500, 14000, 900, 900, 80, 80]
    out = tmp_path / "view.png"
    r = render_view(43.4524508, -6.0705374, dt.date(2026, 8, 12),
                    az, alt, rng, str(out), elev_m=117.6, tz=2.0)
    assert out.exists() and out.stat().st_size > 10000
    assert r["clearance_deg"] == pytest.approx(9.8, abs=0.3)
    assert r["verdict"].startswith("VISIBLE")


def test_labels_use_a_real_scalable_font():
    """A hardcoded Linux font path fell through to load_default(), which
    ignores the size, so every label on macOS came out at 10 px."""
    from PIL import Image, ImageDraw
    from eclipse_sites.render import _font
    d = ImageDraw.Draw(Image.new("RGB", (10, 10)))
    small = d.textbbox((0, 0), "TOTALITY", font=_font(26))
    big = d.textbbox((0, 0), "TOTALITY", font=_font(42, True))
    assert big[2] > small[2] * 1.3        # 42 px must be visibly wider than 26
    assert big[2] > 120                   # and not a 10 px default


def test_render_refuses_when_the_fan_misses_the_sun(tmp_path):
    """np.interp clamps outside the fan, which used to invent a flat skyline
    across ground nobody surveyed and report a clearance from it."""
    from eclipse_sites.render import render_view
    with pytest.raises(ValueError, match="does not cover"):
        render_view(43.4524508, -6.0705374, dt.date(2026, 8, 12),
                    [100, 110, 120], [2.0, 2.0, 2.0], [3500, 3500, 3500],
                    str(tmp_path / "x.png"), elev_m=117.6, tz=2.0)


def test_render_is_indifferent_to_which_branch_the_fan_is_on(tmp_path):
    """A fan is a set of directions, not a set of numbers. 250 deg and -110 deg
    are the same direction and must render the same, which is what lets a view
    across due north work at all."""
    from eclipse_sites.render import render_view
    alt = [5.84, 6.68, 1.85, 1.65, 2.67, 0.64, 1.97]
    rng = [3500, 3500, 3500, 7000, 3500, 3500, 14000]
    plain = list(np.arange(250.0, 311.0, 10.0))
    shifted = [a - 360.0 for a in plain]          # same directions, other branch
    args = dict(elev_m=117.6, tz=2.0)
    a = render_view(43.4524508, -6.0705374, dt.date(2026, 8, 12), plain, alt,
                    rng, str(tmp_path / "a.png"), **args)
    b = render_view(43.4524508, -6.0705374, dt.date(2026, 8, 12), shifted, alt,
                    rng, str(tmp_path / "b.png"), **args)
    assert a["clearance_deg"] == pytest.approx(b["clearance_deg"], abs=1e-9)
    assert a["horizon_deg"] == pytest.approx(b["horizon_deg"], abs=1e-9)
    assert (tmp_path / "b.png").stat().st_size > 10000


def test_render_refuses_when_no_eclipse(tmp_path):
    from eclipse_sites.render import render_view
    with pytest.raises(ValueError):
        render_view(43.45, -6.07, dt.date(2026, 8, 20), [270, 280, 290],
                    [1, 1, 1], [1000, 1000, 1000], str(tmp_path / "x.png"))


# --------------------------------------------------------------- calibrate

def _synth_skyline(intr, az, alt, bearing, pitch):
    # Close the ring exactly as fit_camera does, so a full-circle scene means
    # the same thing to the fixture and to the code under test. Without this a
    # frame that spans north is synthesised from a scene clamped at the seam
    # and no solver can recover the camera that made it.
    from eclipse_sites.calibrate import pixel_to_altaz, _ring
    az, alt, _ = _ring(az, alt, np.zeros_like(np.asarray(az, float)))
    rows = np.arange(0, intr.height, 1.0)
    cols = np.arange(0, intr.width, 4)
    ys = np.full(cols.shape, np.nan)
    for i, x in enumerate(cols):
        a, z = pixel_to_altaz(np.full(rows.shape, float(x)), rows, intr,
                              bearing, pitch)
        s = np.where(np.diff(np.sign(
            a - np.interp(z, az, alt, left=np.nan, right=np.nan))) != 0)[0]
        if s.size:
            ys[i] = rows[s[0]]
    out = np.full(intr.width, np.nan)
    out[cols] = ys
    return out


def _scene():
    az = np.arange(200.0, 340.0, 1.0)
    alt = 3.0 + 5.0 * np.sin(np.radians((az - 200) * 3.0)) + 0.02 * (az - 260)
    return az, alt


def test_fit_recovers_known_camera():
    from eclipse_sites.calibrate import Intrinsics, fit_camera
    intr = Intrinsics(4000, 1848, 1500.0, 2000.0, 924.0)
    az, alt = _scene()
    sky = _synth_skyline(intr, az, alt, 283.0, 14.7)
    f = fit_camera(sky, intr, az, alt, np.full_like(az, 5000.0),
                   bearing_prior=270.0, bearing_window=60.0)
    assert f.bearing_deg == pytest.approx(283.0, abs=0.6)
    assert f.pitch_deg == pytest.approx(14.7, abs=0.4)


def test_canopy_does_not_corrupt_a_correctly_anchored_fit():
    """Trees on near ground must land in the residual, not in the camera angles."""
    from eclipse_sites.calibrate import Intrinsics, fit_camera
    intr = Intrinsics(4000, 1848, 1500.0, 2000.0, 924.0)
    az, alt = _scene()
    rng = np.full_like(az, 5000.0)
    rng[az > 300] = 80.0
    treed = alt.copy()
    treed[az > 300] += 8.0
    sky = _synth_skyline(intr, az, treed, 283.0, 14.7)

    good = fit_camera(sky, intr, az, alt, rng, bearing_prior=270.0,
                      bearing_window=60.0, anchor_min_range_m=2000.0)
    bad = fit_camera(sky, intr, az, alt, rng, bearing_prior=270.0,
                     bearing_window=60.0, anchor_min_range_m=0.0)
    assert good.bearing_deg == pytest.approx(283.0, abs=0.8)
    assert abs(bad.bearing_deg - 283.0) > 5.0    # the documented failure mode


def test_canopy_height_recovered_in_metres():
    from eclipse_sites.calibrate import (Intrinsics, fit_camera,
                                         canopy_from_residual)
    intr = Intrinsics(4000, 1848, 1500.0, 2000.0, 924.0)
    az, alt = _scene()
    rng = np.full_like(az, 5000.0)
    rng[az > 300] = 80.0
    treed = alt.copy()
    treed[az > 300] += 8.0
    sky = _synth_skyline(intr, az, treed, 283.0, 14.7)
    f = fit_camera(sky, intr, az, alt, rng, bearing_prior=270.0,
                   bearing_window=60.0)
    rows = canopy_from_residual(sky, intr, f, az, alt, rng, eye_elev_m=117.6,
                                az_bin=5.0)
    near = [r["extra_height_m"] for r in rows if r["az_deg"] > 302]
    expect = 80 * math.tan(math.radians(8))      # 11.2 m
    assert np.median(near) == pytest.approx(expect, abs=1.5)


def test_canopy_profile_round_trip_reaches_the_horizon():
    """The whole point of the feature: measured trees change a clearance.

    Same synthesised camera as the metres test, carried through the profile and
    into horizon_point, which must move the horizon by the measured height.
    """
    from eclipse_sites.calibrate import (Intrinsics, fit_camera,
                                         canopy_from_residual, canopy_profile,
                                         profile_to_azimuths)
    intr = Intrinsics(4000, 1848, 1500.0, 2000.0, 924.0)
    az, alt = _scene()
    rng = np.full_like(az, 5000.0)
    rng[az > 300] = 80.0
    treed = alt.copy()
    treed[az > 300] += 8.0
    sky = _synth_skyline(intr, az, treed, 283.0, 14.7)
    f = fit_camera(sky, intr, az, alt, rng, bearing_prior=270.0,
                   bearing_window=60.0)
    rows = canopy_from_residual(sky, intr, f, az, alt, rng, eye_elev_m=117.6,
                                az_bin=5.0)
    prof = canopy_profile(rows, f, max_rms_deg=5.0)
    got = profile_to_azimuths([305.0], prof, max_gap_deg=8.0)[0]
    assert got == pytest.approx(80 * math.tan(math.radians(8)), abs=2.0)

    # the angle is what gets applied, and it must exceed bare earth by the
    # 8 deg of trees that were injected
    ang = profile_to_azimuths([305.0], prof, max_gap_deg=8.0,
                              key="skyline_alt_deg")[0]
    bare = float(np.interp(305.0, az, alt))
    assert ang - bare == pytest.approx(8.0, abs=1.5)


def test_photo_is_applied_as_an_angle_not_as_metres_of_canopy():
    """Measured metres must never be fed to the engine as canopy_m.

    canopy_m is added to every ground sample along the ray, including the one
    half a cell from the eye. Twenty metres of real, correctly measured tree
    lands there as an 83 degree wall and the site reads as blocked. The
    photograph measured an angle; the round trip through metres is what breaks
    it, so the profile carries the angle and site applies it as a floor.
    """
    from eclipse_sites.calibrate import profile_to_azimuths
    z = np.zeros((N, N), np.float32)
    z[:, C + 50:] = 100.0
    bare = horizon_point(z, CELL, C, C, [90], max_range_m=1500)[90].alt_deg
    as_metres = horizon_point(z, CELL, C, C, [90], max_range_m=1500,
                              canopy_m=20.0)[90].alt_deg
    assert as_metres > 70.0              # the wall, which is why we do not

    prof = {"canopy": [{"az_deg": a, "skyline_alt_deg": 15.0, "canopy_m": 20.0}
                       for a in (88.0, 90.0, 92.0)]}
    floor = profile_to_azimuths([90.0], prof, fallback=-90.0,
                                key="skyline_alt_deg")[0]
    assert bare == pytest.approx(11.1, abs=0.2)      # the wall at 500 m
    assert max(bare, floor) == pytest.approx(15.0, abs=1e-6)


def test_canopy_profile_drops_what_it_cannot_convert():
    """Range, sign and plausibility filters, all at write time."""
    from eclipse_sites.calibrate import canopy_profile
    fit = type("F", (), {"bearing_deg": 280.0, "pitch_deg": 10.0,
                         "rms_deg": 0.4, "n_anchor": 200})()
    rows = [
        {"az_deg": 10.0, "photo_alt_deg": 1.0, "dem_alt_deg": 0.5,
         "residual_deg": +0.2, "horizon_range_m": 18700,
         "extra_height_m": 65.0},                       # horizon far away
        {"az_deg": 12.0, "photo_alt_deg": -1.0, "dem_alt_deg": 0.5,
         "residual_deg": -1.2, "horizon_range_m": 500,
         "extra_height_m": -9.0},                       # trace or fit is wrong
        {"az_deg": 14.0, "photo_alt_deg": 5.0, "dem_alt_deg": 0.5,
         "residual_deg": +4.0, "horizon_range_m": 900,
         "extra_height_m": 300.0},                      # no tree is 300 m
        {"az_deg": 16.0, "photo_alt_deg": 1.0, "dem_alt_deg": 0.5,
         "residual_deg": -0.05, "horizon_range_m": 900,
         "extra_height_m": -0.4},                       # noise, clamps to zero
        {"az_deg": 18.0, "photo_alt_deg": 2.0, "dem_alt_deg": 0.5,
         "residual_deg": +1.0, "horizon_range_m": 900,
         "extra_height_m": 15.0},                       # a real measurement
    ]
    p = canopy_profile(rows, fit)
    assert [e["az_deg"] for e in p["canopy"]] == [16.0, 18.0]
    assert p["canopy"][0]["canopy_m"] == 0.0            # clamped, not negative
    d = p["meta"]["bins_dropped"]
    assert (d["horizon_too_far"], d["residual_negative"],
            d["height_implausible"]) == (1, 1, 1)


def test_canopy_profile_refuses_a_fit_it_does_not_believe():
    from eclipse_sites.calibrate import canopy_profile
    bad = type("F", (), {"bearing_deg": 0.0, "pitch_deg": 0.0,
                         "rms_deg": 3.5, "n_anchor": 40})()
    with pytest.raises(ValueError, match="rms"):
        canopy_profile([], bad)


def test_a_dropped_bin_does_not_come_back_as_interpolation():
    """The filters are worthless if np.interp rebuilds what they discarded.

    A bin dropped because its horizon forms 18 km away must not reappear as the
    average of its neighbours. Only the width of the hole can decide this: an
    azimuth 3 deg from a surviving bin is still inside an 8 deg gap.
    """
    from eclipse_sites.calibrate import profile_to_azimuths
    prof = {"canopy": [{"az_deg": 277.0, "canopy_m": 20.0},
                       {"az_deg": 285.0, "canopy_m": 20.0}]}
    assert profile_to_azimuths([277.0], prof, max_gap_deg=4.0)[0] == 20.0
    assert profile_to_azimuths([280.6], prof, max_gap_deg=4.0)[0] == 0.0
    # widen the tolerance past the hole and interpolation is allowed again
    assert profile_to_azimuths([280.6], prof, max_gap_deg=10.0)[0] == 20.0


def test_canopy_profile_falls_back_outside_the_photographed_arc():
    """A profile is an arc, not a ring, and says nothing behind the camera."""
    from eclipse_sites.calibrate import profile_to_azimuths
    prof = {"canopy": [{"az_deg": 350.0, "canopy_m": 10.0},
                       {"az_deg": 353.0, "canopy_m": 10.0},
                       {"az_deg": 356.0, "canopy_m": 10.0},
                       {"az_deg": 2.0, "canopy_m": 10.0}]}
    inside = profile_to_azimuths([352.0, 355.0], prof, fallback=99.0)
    assert list(inside) == [10.0, 10.0]                 # arc crossing north
    assert profile_to_azimuths([180.0], prof, fallback=99.0)[0] == 99.0


def test_coverage_penalty_blocks_the_sliver_overfit():
    """Without it the solver slides the frame off the anchors and 'wins'."""
    from eclipse_sites.calibrate import Intrinsics, fit_camera
    intr = Intrinsics(4000, 1848, 1500.0, 2000.0, 924.0)
    az, alt = _scene()
    sky = _synth_skyline(intr, az, alt, 283.0, 14.7)
    f = fit_camera(sky, intr, az, alt, np.full_like(az, 5000.0),
                   bearing_prior=270.0, bearing_window=60.0)
    usable = int(np.isfinite(sky[::8]).sum())
    assert f.n_anchor > 0.25 * usable


def test_intrinsics_reads_exif_sub_ifd(tmp_path):
    """Phones store FocalLengthIn35mmFilm in the Exif IFD, not the top level."""
    from PIL import Image
    from eclipse_sites.calibrate import intrinsics_from_exif
    p = tmp_path / "x.jpg"
    im = Image.new("RGB", (4000, 1848), (128, 128, 128))
    ex = im.getexif()
    ex.get_ifd(0x8769)[0xA405] = 13
    im.save(p, exif=ex)
    intr = intrinsics_from_exif(str(p))
    assert 1400 < intr.f_px < 1600
    assert intr.source.startswith("exif")


def test_intrinsics_overrides(tmp_path):
    """--f-px and --hfov must bypass EXIF entirely, including when it is absent."""
    from PIL import Image
    from eclipse_sites.calibrate import intrinsics_from_exif
    p = tmp_path / "bare.jpg"
    Image.new("RGB", (4000, 1848), (128, 128, 128)).save(p)   # no EXIF at all

    with pytest.raises(ValueError):
        intrinsics_from_exif(str(p))

    i = intrinsics_from_exif(str(p), f_px=1234.0)
    assert i.f_px == 1234.0 and i.cx == 2000 and i.cy == 924
    assert "f_px" in i.source

    i = intrinsics_from_exif(str(p), hfov_deg=90.0)
    assert i.f_px == pytest.approx(2000.0)      # half width / tan(45 deg)
    assert "hfov" in i.source


def test_reported_residual_is_a_true_rms_not_the_search_cost():
    """The CLI prints rms_deg as degrees and warns above 2.0, so it must be
    degrees. It used to carry the sqrt(n_cols/n) coverage penalty baked in."""
    from eclipse_sites.calibrate import Intrinsics, fit_camera, pixel_to_altaz
    intr = Intrinsics(4000, 1848, 1500.0, 2000.0, 924.0)
    az, alt = _scene()
    sky = _synth_skyline(intr, az, alt, 283.0, 14.7)
    f = fit_camera(sky, intr, az, alt, np.full_like(az, 5000.0),
                   bearing_prior=270.0, bearing_window=60.0)

    cols = np.arange(0, intr.width, 8)
    ys = sky[cols]
    ok = np.isfinite(ys)
    a, z = pixel_to_altaz(cols[ok], ys[ok], intr, f.bearing_deg, f.pitch_deg,
                          f.roll_deg, f.f_scale)
    m = (z >= az[0]) & (z <= az[-1])
    r = np.interp(z[m], az, alt) - a[m]
    assert f.rms_deg == pytest.approx(float(np.sqrt(np.mean(r * r))), abs=1e-6)
    assert f.cost >= f.rms_deg          # the penalty is kept, just separately


def _ring_scene():
    """A full-circle skyline with distinct landmarks, so the fit is unique.

    A periodic scene like sin(3*az) fits equally well at three bearings 120
    apart, which is a real ambiguity and not something a solver should be
    asked to resolve.
    """
    az = np.arange(0.0, 360.0, 1.0)
    d = lambda c: (az - c + 180.0) % 360.0 - 180.0
    alt = (2.0 + 9.0 * np.exp(-d(20.0) ** 2 / 120.0)
           + 6.0 * np.exp(-d(350.0) ** 2 / 60.0)
           + 4.0 * np.exp(-d(120.0) ** 2 / 300.0))
    return az, alt


def test_anchor_sectors_may_wrap_through_north():
    """hi < lo means the sector runs through north, e.g. (340, 40)."""
    from eclipse_sites.calibrate import Intrinsics, fit_camera
    intr = Intrinsics(4000, 1848, 1500.0, 2000.0, 924.0)
    az, alt = _ring_scene()
    sky = _synth_skyline(intr, az, alt, 10.0, 14.7)
    f = fit_camera(sky, intr, az, alt, np.full_like(az, 5000.0),
                   bearing_prior=10.0, bearing_window=25.0,
                   anchor_sectors=[(340.0, 40.0)])
    assert f.n_anchor > 0
    assert f.bearing_deg == pytest.approx(10.0, abs=1.5)


def test_full_circle_fan_has_no_seam_at_north():
    """A 0 to 359.5 fan has a hole in it as far as np.interp is concerned, and
    anything landing in that hole clamps to an endpoint instead of wrapping."""
    from eclipse_sites.calibrate import _ring
    az = np.arange(0.0, 360.0, 0.5)
    alt = np.where(az < 180.0, 5.0, 9.0)
    raz, ralt, _ = _ring(az, alt, np.full_like(az, 5000.0))
    assert raz[0] < 0.0 and raz[-1] > 359.5      # seam closed on both sides
    # Crossing the seam must interpolate from the 9 side round to the 5 side.
    # Without the ring closed, both of these clamp to an endpoint instead.
    assert 5.0 < np.interp(359.75, raz, ralt) < 9.0
    assert 5.0 < np.interp(-0.25, raz, ralt) < 9.0
    assert np.interp(359.75, az, alt) == 9.0     # the clamped, wrong answer


def test_extract_skyline_finds_a_hard_edge(tmp_path):
    from PIL import Image
    from eclipse_sites.calibrate import extract_skyline
    a = np.zeros((200, 300, 3), np.uint8)
    a[:80] = 220          # bright grey sky
    a[80:] = 30           # dark ground
    p = tmp_path / "s.png"
    Image.fromarray(a).save(p)
    prof = extract_skyline(str(p), smooth=1)
    assert np.nanmedian(prof) == pytest.approx(80, abs=2)


def test_masked_columns_are_ignored(tmp_path):
    from PIL import Image
    from eclipse_sites.calibrate import extract_skyline
    a = np.zeros((200, 300, 3), np.uint8)
    a[:80] = 220
    a[80:] = 30
    a[:, :50] = 30        # a pole or parasol filling the left edge
    p = tmp_path / "m.png"
    Image.fromarray(a).save(p)
    prof = extract_skyline(str(p), masks=[(0, 0, 50, 200)], smooth=1)
    assert np.isnan(prof[:50]).all()
    assert np.nanmedian(prof[60:]) == pytest.approx(80, abs=2)


# ------------------------------------------------------------------ search

def _proj_dem(path, z, cell, x0, y0, crs="EPSG:25829"):
    """Write a projected GeoTIFF. x0, y0 is the top left corner in metres."""
    import rasterio
    from rasterio.transform import from_origin
    with rasterio.open(path, "w", driver="GTiff", height=z.shape[0],
                       width=z.shape[1], count=1, dtype="float32", crs=crs,
                       transform=from_origin(x0, y0, cell, cell),
                       nodata=-9999.0) as d:
        d.write(z.astype("float32"), 1)


# UTM 29N metres near Candamo, so the 2026 eclipse is total over the raster
X0, Y0, DCELL, DN = 690000.0, 4818000.0, 25.0, 240


RANGE = 1000.0          # 40 cells at 25 m, so columns below 40 are unknowable


def _region(ridge_h=0.0):
    """Flat ground, optionally with a north-south ridge at columns 60 to 70.

    The ridge sits well inside the western edge so that cells east of it still
    have room to march a full RANGE without leaving the raster.
    """
    z = np.zeros((DN, DN), np.float32)
    if ridge_h:
        z[:, 60:70] = ridge_h          # low column index = west
    return z


def test_az_bins_always_gives_at_least_two():
    from eclipse_sites.search import _az_bins
    single = np.full((4, 4), 280.0)          # lands exactly on a bin edge
    _, lo, azs = _az_bins(single, 2.5)
    assert len(azs) >= 2
    assert lo == 280.0


def test_az_bins_unwraps_a_region_straddling_north():
    from eclipse_sites.search import _az_bins
    az = np.array([[359.0, 1.0], [358.5, 0.5]])
    az_lin, lo, azs = _az_bins(az, 2.5)
    assert az_lin.max() - az_lin.min() < 5.0      # not 358
    assert len(azs) < 10                          # not the whole compass
    assert all(a > 350.0 for a in azs)


def test_sweep_runs_and_finds_the_ridge(tmp_path):
    """End to end on a synthetic region. search.py had no test at all."""
    from eclipse_sites.search import sweep, top_sites
    p = tmp_path / "coarse.tif"
    _proj_dem(str(p), _region(ridge_h=300.0), DCELL, X0, Y0)
    clear, z, tr, crs, cell, circ = sweep(
        str(p), dt.date(2026, 8, 12), max_range_m=RANGE, az_step=2.5)
    assert clear.shape == (DN, DN)
    assert cell == DCELL
    assert circ.kind == "total"
    # Sun is at ~280 deg, just north of west. Cells just east of the ridge look
    # straight into it; cells far east cannot reach it within RANGE.
    blocked = np.nanmean(clear[:, 80:100])
    open_ = np.nanmean(clear[:, 180:220])
    assert blocked < -10.0 < 0.0 < open_
    sites = top_sites(clear, z, tr, crs, n=5, min_sep_cells=10)
    assert len(sites) == 5
    assert all(-7 < s["lon"] < -5 and 43 < s["lat"] < 44 for s in sites)
    assert sites[0]["clearance_deg"] >= sites[-1]["clearance_deg"]


def test_sweep_canopy_is_metres_not_degrees(tmp_path):
    """The regression that made the coarse and fine columns incomparable.

    sweep used to add canopy_m straight onto the horizon ANGLE, so --canopy 20
    meant 20 metres of trees on the fine pass and a flat 20 degrees on the
    coarse one. The fix is to hand the metres to the engine, so the drop in
    clearance must be exactly the rise the engine gives, not the raw number.
    """
    from eclipse_sites.search import sweep
    p = tmp_path / "c.tif"
    _proj_dem(str(p), _region(ridge_h=300.0), DCELL, X0, Y0)
    kw = dict(when=dt.date(2026, 8, 12), max_range_m=RANGE, az_step=2.5)
    bare, *_ = sweep(str(p), **kw)
    treed, *_ = sweep(str(p), canopy_m=20.0, **kw)
    drop = (bare - treed)[120, 150]

    # what the engine itself says 20 m of canopy is worth on that ray
    z = _region(ridge_h=300.0)
    az = 280.0
    h0 = horizon_grid(z, DCELL, [az], max_range_m=RANGE)[az][120, 150]
    h1 = horizon_grid(z, DCELL, [az], max_range_m=RANGE,
                      canopy_m=20.0)[az][120, 150]
    assert drop == pytest.approx(h1 - h0, abs=0.5)
    assert drop != pytest.approx(20.0, abs=1.0)   # the old, wrong answer


def test_sweep_edge_cells_are_masked_not_wide_open(tmp_path):
    """Invariant 1, at the level where it actually bit."""
    from eclipse_sites.search import sweep
    p = tmp_path / "e.tif"
    _proj_dem(str(p), _region(), DCELL, X0, Y0)
    clear, *_ = sweep(str(p), dt.date(2026, 8, 12), max_range_m=RANGE)
    # The Sun is at ~280 deg, so rays run west and slightly north. Cells within
    # RANGE of the western edge cannot see where they are looking.
    assert np.isnan(clear[20:, :30]).all()       # ray exits to the west
    assert np.isfinite(clear[20:, 60:]).all()    # room to march
    assert np.nanmax(clear) < 90.0               # no +96 deg nonsense
    assert np.nanmax(clear) == pytest.approx(10.5, abs=1.0)   # flat, Sun alt


def test_confirm_reports_sites_the_fine_dem_misses(tmp_path):
    """They used to vanish silently, leaving holes in the rank column."""
    from eclipse_sites.search import confirm
    fine = tmp_path / "fine.tif"
    _proj_dem(str(fine), _region(), 5.0, X0, Y0)     # only 1.2 km across
    from pyproj import Transformer
    to_wgs = Transformer.from_crs("EPSG:25829", "EPSG:4326", always_xy=True)
    inside = to_wgs.transform(X0 + 600, Y0 - 600)
    outside = to_wgs.transform(X0 + 50000, Y0 - 600)
    sites = [{"rank": 1, "lat": inside[1], "lon": inside[0], "elev_m": 0.0},
             {"rank": 2, "lat": outside[1], "lon": outside[0], "elev_m": 0.0}]
    ok, dropped = confirm(str(fine), sites, dt.date(2026, 8, 12),
                          max_range_m=1000.0)
    assert len(ok) == 1 and len(dropped) == 1
    assert "outside" in dropped[0]["why"]
    assert ok[0]["rank"] == 1                     # renumbered after re-sorting
    assert ok[0]["ray_truncated"] is True         # 1.2 km raster, 1 km ray


# ----------------------------------------------------------------- prepare

def _tile(path, crs, x0, y0, cell, n=60, base=100.0):
    import rasterio
    from rasterio.transform import from_origin
    z = (base + 20 * np.sin(np.arange(n * n).reshape(n, n) / 37.0)).astype("float32")
    with rasterio.open(path, "w", driver="GTiff", height=n, width=n, count=1,
                       dtype="float32", crs=crs,
                       transform=from_origin(x0, y0, cell, cell)) as d:
        d.write(z, 1)


def _dump(tmp_path):
    import zipfile
    d = tmp_path / "dump"
    (d / "sub").mkdir(parents=True)
    _tile(str(d / "PNOA_MDT05_ETRS89_HU29_0026_LID.tif"), "EPSG:25829", 700000, 4820000, 5)
    _tile(str(d / "sub" / "PNOA_MDT05_ETRS89_HU29_0027_LID.tif"), "EPSG:25829", 700300, 4820000, 5)
    _tile(str(d / "PNOA_MDT25_ETRS89_HU29_0026.tif"), "EPSG:25829", 700000, 4820000, 25)
    _tile(str(d / "PNOA_MDSNV25_ETRS89_HU29_0026.tif"), "EPSG:25829", 700000, 4820000, 2.5, base=12.0)
    (d / "readme.pdf").write_text("not a raster")
    z = tmp_path / "t.tif"
    _tile(str(z), "EPSG:25829", 700600, 4820000, 5)
    with zipfile.ZipFile(d / "extra.zip", "w") as zf:
        zf.write(z, "PNOA_MDT05_HU29_0030_LID.tif")
    return d


def test_prepare_sorts_merges_and_unzips(tmp_path):
    from eclipse_sites.prepare import prepare
    rep = prepare(str(_dump(tmp_path)), str(tmp_path / "out"))
    got = {w["product"]: w for w in rep.written}
    assert set(got) == {"dtm05", "dtm25", "ndsm_vegetation"}
    assert got["dtm05"]["tiles"] == 3          # 2 loose + 1 from the zip
    assert rep.target_crs == "EPSG:25829"


def test_prepare_ignores_non_rasters(tmp_path):
    from eclipse_sites.prepare import scan
    tiles = scan(str(_dump(tmp_path)))
    assert all(not t.path.endswith(".pdf") for t in tiles)


def test_normalised_products_use_max_resampling(tmp_path):
    """Averaging a treeline loses it; the horizon needs the tallest thing."""
    from eclipse_sites.prepare import prepare
    rep = prepare(str(_dump(tmp_path)), str(tmp_path / "out"))
    veg = [w for w in rep.written if w["product"] == "ndsm_vegetation"][0]
    assert veg["resampling"] == "max"
    assert veg["normalised_heights"] is True
    assert any("not altitudes" in w for w in rep.warnings)


def test_gaps_are_declared_nodata_not_zero(tmp_path):
    """A gap read as 0 m would look like a cliff and fake a clear view."""
    import rasterio
    from eclipse_sites.prepare import prepare
    d = tmp_path / "gappy"
    d.mkdir()
    # deliberately NOT adjacent: 60 cells at 5 m spans 300 m, so leaving 2 km
    # between the two tiles guarantees empty ground in the mosaic
    _tile(str(d / "PNOA_MDT05_HU29_a.tif"), "EPSG:25829", 700000, 4820000, 5)
    _tile(str(d / "PNOA_MDT05_HU29_b.tif"), "EPSG:25829", 702000, 4820000, 5)
    out = tmp_path / "out"
    prepare(str(d), str(out))
    with rasterio.open(out / "dtm05.tif") as f:
        assert f.nodata is not None
        raw = f.read(1)
        assert (raw == f.nodata).sum() > 0             # the gap exists
        assert not ((raw > -1) & (raw < 1)).any()      # and is not 0 m ground
        arr = f.read(1, masked=True).filled(np.nan)
        assert np.isnan(arr).sum() > 0                 # reads as NaN downstream


def test_product_is_classified_by_filename_not_by_parent_directory(tmp_path):
    """A folder named for one product must not reclassify the files inside it.

    scan() matched the basename but prepare() matched the whole path, so a
    vegetation model sitting in a folder called MDT25_download was labelled
    ndsm_vegetation and then resampled with bilinear, quietly averaging away
    the treeline that MAX exists to preserve.
    """
    from eclipse_sites.prepare import prepare
    d = tmp_path / "MDT25_download"
    d.mkdir()
    _tile(str(d / "PNOA_MDSNV25_ETRS89_HU29_0026.tif"), "EPSG:25829",
          700000, 4820000, 2.5, base=12.0)
    rep = prepare(str(d), str(tmp_path / "out"))
    veg = [w for w in rep.written if w["product"] == "ndsm_vegetation"]
    assert len(veg) == 1
    assert veg[0]["resampling"] == "max"
    assert veg[0]["normalised_heights"] is True


def test_oversized_mosaic_is_refused_before_it_is_allocated(tmp_path):
    """merge() builds the whole array in memory, so scattered tiles used to be
    an OOM kill. The estimate now comes from the warped bounds."""
    from eclipse_sites.prepare import prepare
    d = tmp_path / "scattered"
    d.mkdir()
    # Two 5 m tiles 400 km apart on BOTH axes. Separating them on one axis
    # only gives a long thin mosaic that is quite cheap; it is the area that
    # kills you.
    _tile(str(d / "PNOA_MDT05_HU29_a.tif"), "EPSG:25829", 700000, 4820000, 5)
    _tile(str(d / "PNOA_MDT05_HU29_b.tif"), "EPSG:25829", 1100000, 4420000, 5)
    rep = prepare(str(d), str(tmp_path / "out"), max_cells=4e8)
    assert not rep.written
    assert any("over the" in w and "limit" in w for w in rep.warnings)

    # The limit is a limit, not a refusal to work: adjacent tiles of the same
    # resolution pass comfortably. Raising it on the scattered pair instead
    # would make this test allocate the 25 GB the check exists to prevent.
    close = tmp_path / "adjacent"
    close.mkdir()
    _tile(str(close / "PNOA_MDT05_HU29_a.tif"), "EPSG:25829", 700000, 4820000, 5)
    _tile(str(close / "PNOA_MDT05_HU29_b.tif"), "EPSG:25829", 700300, 4820000, 5)
    rep2 = prepare(str(close), str(tmp_path / "out2"), max_cells=4e8)
    assert [w["product"] for w in rep2.written] == ["dtm05"]
    rep3 = prepare(str(close), str(tmp_path / "out3"), max_cells=100)
    assert not rep3.written


def test_zone_inferred_from_filename_when_crs_missing(tmp_path):
    from eclipse_sites.prepare import _zone_from_name
    assert _zone_from_name("PNOA_MDT05_ETRS89_HU29_0026_LID.asc") == "EPSG:25829"
    assert _zone_from_name("PNOA_MDT05_ETRS89_HU30_0028_LID.asc") == "EPSG:25830"
    assert _zone_from_name("random_file.tif") is None


def test_src_crs_override_rescues_projectionless_input(tmp_path):
    from eclipse_sites.prepare import scan
    d = tmp_path / "d"
    d.mkdir()
    _tile(str(d / "mystery_MDT05.tif"), None, 700000, 4820000, 5)
    assert scan(str(d))[0].crs is None
    assert scan(str(d), src_crs="EPSG:25829")[0].crs == "EPSG:25829"
