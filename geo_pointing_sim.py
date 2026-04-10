#!/usr/bin/env python3
"""
geo_pointing_sim.py
===================
Simulates an optimal sidereal-tracking pointing strategy for a multi-camera
array imaging the geostationary (GEO) satellite belt.

Cameras are arranged in a 1×n row along the ecliptic (numbered 1→n, West→East).
The mount tracks at sidereal rate for --dwell minutes, then slews East by one
full array-width to the next pointing.  The simulation runs for four key dates
(vernal equinox, summer solstice, autumnal equinox, winter solstice) from
Elginfield Observatory (43.1925°N, 81.3158°W, 325 m).

Usage examples
--------------
  python geo_pointing_sim.py                              # quick summary, all 4 dates
  python geo_pointing_sim.py --plot                      # + save PNG figures
  python geo_pointing_sim.py --n-cameras 6 --fov 3.5    # 6 cameras, 3.5° FOV each
  python geo_pointing_sim.py --refresh                   # re-download GEO TLEs first

Requirements: skyfield, matplotlib, numpy  (see requirements.txt)
"""

import argparse
import sys
import os
import itertools
from collections import defaultdict
from datetime import datetime, timedelta, date, timezone

import numpy as np

try:
    from skyfield.api import load, wgs84, EarthSatellite
    from skyfield import almanac
    from skyfield.framelib import ecliptic_frame
except ImportError:
    sys.exit("Error: 'skyfield' is required.  Install with:  pip install skyfield")

# ---------------------------------------------------------------------------
# Constants – Elginfield Observatory
# ---------------------------------------------------------------------------
ELGINFIELD_LAT  =  43.1925   # degrees North
ELGINFIELD_LON  = -81.3158   # degrees East (negative = West)
ELGINFIELD_ELEV =  325.0     # metres

# Four simulation dates (month, day) – always use current --year
SIM_DATES = [
    ("Vernal Equinox",    3, 20),
    ("Summer Solstice",   6, 21),
    ("Autumnal Equinox",  9, 23),
    ("Winter Solstice",  12, 21),
]

CELESTRAK_GEO_URL = "https://celestrak.org/NORAD/elements/gp.php?GROUP=geo&FORMAT=TLE"
DEFAULT_TLE_FILE  = "tle_group_geo.txt"


# ---------------------------------------------------------------------------
# TLE loading
# ---------------------------------------------------------------------------

def load_tle_file(tle_path, ts):
    """Load EarthSatellite objects from a local TLE file."""
    if not os.path.exists(tle_path):
        return []
    sats = []
    with open(tle_path, "r") as f:
        lines = [l.rstrip() for l in f if l.strip()]
    i = 0
    while i + 2 < len(lines):
        name = lines[i]
        line1 = lines[i + 1]
        line2 = lines[i + 2]
        if line1.startswith("1 ") and line2.startswith("2 "):
            try:
                sats.append(EarthSatellite(line1, line2, name, ts))
                i += 3
            except Exception:
                i += 1
        else:
            i += 1
    return sats


def refresh_tles(tle_path):
    """Download fresh GEO TLEs from Celestrak."""
    try:
        import requests
        print(f"Downloading GEO TLEs from Celestrak…")
        r = requests.get(CELESTRAK_GEO_URL, timeout=30)
        r.raise_for_status()
        with open(tle_path, "w") as f:
            f.write(r.text)
        print(f"  Saved {tle_path}")
    except Exception as e:
        print(f"  Warning: could not refresh TLEs – {e}")


# ---------------------------------------------------------------------------
# Twilight / dark-night window
# ---------------------------------------------------------------------------

def get_night_window(observer, ts, year, month, day, min_alt_deg=10.0):
    """
    Return (t_start, t_end) skyfield Time objects bracketing astronomical night
    (sun below −18° altitude) for the given date at Elginfield.
    Falls back to −12° (nautical) then −6° (civil) if needed.
    Also clips to mount's minimum altitude constraint on the GEO arc.
    """
    eph = load("de421.bsp")

    # Search window: noon local → noon next day  (UTC offset ≈ −5 h)
    utc_offset = 5  # hours behind UTC
    noon_utc = ts.utc(year, month, day, 17, 0, 0)   # ~noon local
    noon_utc2 = ts.utc(year, month, day + 1, 17, 0, 0)

    f = almanac.dark_twilight_day(eph, observer)

    try:
        times, events = almanac.find_discrete(noon_utc, noon_utc2, f)
    except Exception:
        # fallback: use a 10-hour window centred on local midnight
        midnight = ts.utc(year, month, day + 1, 5, 0, 0)
        t_start = ts.utc(year, month, day + 1, 0, 0, 0)
        t_end   = ts.utc(year, month, day + 1, 10, 0, 0)
        return t_start, t_end

    # almanac dark_twilight_day states:
    # 0 = night, 1 = astronomical twilight, 2 = nautical, 3 = civil, 4 = day
    # We want transitions INTO night (event=0) and OUT OF night (event!=0 following 0)
    # --- simplified: find the longest contiguous stretch where f <= 1 (astro night)
    night_start = None
    night_end   = None
    for t, ev in zip(times, events):
        if ev == 0 and night_start is None:
            night_start = t
        if ev != 0 and night_start is not None and night_end is None:
            night_end = t

    if night_start is None:
        # No astronomical night (high summer); use nautical twilight
        for t, ev in zip(times, events):
            if ev <= 1 and night_start is None:
                night_start = t
            if ev > 1 and night_start is not None and night_end is None:
                night_end = t

    if night_start is None:
        # Ultimate fallback: 10 PM – 5 AM UTC (local ~6 PM – midnight)
        night_start = ts.utc(year, month, day + 1, 1, 0, 0)
        night_end   = ts.utc(year, month, day + 1, 10, 0, 0)

    if night_end is None:
        night_end = ts.utc(year, month, day + 1, 11, 0, 0)

    return night_start, night_end


# ---------------------------------------------------------------------------
# GEO belt / ecliptic pointing geometry
# ---------------------------------------------------------------------------

def ecliptic_to_radec(ecl_lon_deg, ecl_lat_deg, ts_ref):
    """Convert ecliptic (lon, lat) in degrees to (RA, Dec) degrees at J2000."""
    from skyfield.positionlib import Astrometric
    import math
    # Obliquity of ecliptic J2000
    eps = math.radians(23.4393)
    l = math.radians(ecl_lon_deg)
    b = math.radians(ecl_lat_deg)
    # Rotation
    ra = math.atan2(math.sin(l) * math.cos(eps) - math.tan(b) * math.sin(eps),
                    math.cos(l))
    dec = math.asin(math.sin(b) * math.cos(eps) + math.cos(b) * math.sin(eps) * math.sin(l))
    return math.degrees(ra) % 360, math.degrees(dec)


def build_pointing_sequence(n_cameras, fov_deg, gap_deg, t_start, t_end,
                             observer, ts, eph, min_alt_deg=10.0):
    """
    Build the list of mount pointing centres (RA, Dec) in degrees.

    Strategy:
    - The GEO belt runs along the celestial equator / ecliptic.
    - We sweep ecliptic longitude from roughly 270° (rising) to 90° (setting)
      where the belt is above min_alt as seen from Elginfield.
    - Step size = n_cameras × (fov_deg + gap_deg) in RA.
    - Cameras are numbered 1→n from West→East relative to the mount centre.
    """
    array_width_deg = n_cameras * (fov_deg + gap_deg) - gap_deg  # total width

    # Sample the night in 1-min steps and find the range of LST / RA where
    # the ecliptic is above min_alt. We do this by checking a point on the
    # celestial equator (dec=0) across RA sweeps.
    sun = eph["sun"]
    earth = eph["earth"]

    # Find the range of RA (centre) that is above min_alt over the night
    t_mid_jd = (t_start.tt + t_end.tt) / 2.0
    t_mid = ts.tt_jd(t_mid_jd)

    from skyfield.api import Star
    visible_ra = []
    for ra_test in np.arange(0, 360, 2.0):
        star = Star(ra_hours=ra_test / 15.0, dec_degrees=0.0)
        astrometric = (earth + observer).at(t_mid).observe(star)
        alt, az, _ = astrometric.apparent().altaz()
        if alt.degrees >= min_alt_deg:
            visible_ra.append(ra_test)

    if not visible_ra:
        # Fallback: full range
        visible_ra = list(np.arange(0, 360, 2.0))

    ra_min = min(visible_ra)
    ra_max = max(visible_ra)

    # Handle wrap-around (e.g. 300°→60° crossing 0°)
    if ra_max - ra_min > 180:
        # Wrap case: scan up from ra_min → 360, then 0 → ra_max
        centres = []
        ra_cur = ra_min
        while ra_cur <= 360.0:
            centres.append(ra_cur % 360)
            ra_cur += array_width_deg
        ra_cur = 0.0
        while ra_cur <= ra_max:
            centres.append(ra_cur)
            ra_cur += array_width_deg
    else:
        ra_cur = ra_min
        centres = []
        while ra_cur <= ra_max + array_width_deg:
            centres.append(ra_cur % 360)
            ra_cur += array_width_deg

    # Return list of (ra_center, dec_center=0) tuples
    return [(ra, 0.0) for ra in centres]


# ---------------------------------------------------------------------------
# Camera FOV intersection
# ---------------------------------------------------------------------------

def camera_fovs(mount_ra, mount_dec, n_cameras, fov_deg, gap_deg):
    """
    Return list of (ra_center, dec_center, half_fov) for each camera.
    Cameras tile East-West (in RA), numbered 1→n West→East.
    mount_ra / mount_dec = pointing of camera 1 (westernmost).
    The mount slews so that the *array centre* is at mount_ra.
    """
    # Array centre at mount_ra; camera 1 is at the western edge
    step = fov_deg + gap_deg
    total = n_cameras * step - gap_deg
    # Camera 1 centre = mount_ra - total/2 + fov_deg/2
    cam1_ra = mount_ra - total / 2.0 + fov_deg / 2.0
    cams = []
    for i in range(n_cameras):
        cam_ra = (cam1_ra + i * step) % 360
        cams.append((cam_ra, mount_dec, fov_deg / 2.0))
    return cams


def sat_in_fov(sat_ra, sat_dec, fov_list):
    """
    Check if satellite at (sat_ra, sat_dec) falls within any camera FOV.
    FOV is a square box (half_fov×half_fov) centred on (cam_ra, cam_dec).
    Returns camera number (1-based) or 0 if not in any FOV.
    """
    for idx, (cam_ra, cam_dec, half) in enumerate(fov_list):
        dra = abs(((sat_ra - cam_ra + 180) % 360) - 180)  # handle wrap
        ddec = abs(sat_dec - cam_dec)
        if dra <= half and ddec <= half:
            return idx + 1
    return 0


# ---------------------------------------------------------------------------
# Main simulation for one date
# ---------------------------------------------------------------------------

def simulate_night(label, year, month, day, satellites, ts, eph,
                   n_cameras, fov_deg, gap_deg, dwell_min, min_alt_deg,
                   verbose=True):
    """
    Simulate one night and return a results dict with:
      pointings         - list of (ra, dec) unique mount pointing centres
      dwell_records     - list of {pointing_pos, revisit, t_start_utc, n_unique, cam_hits}
      total_dark_min    - total dark-night duration in minutes
      unique_sats_total - # unique GEO sats seen over the whole night
      pct_time_sats     - % of dark time with >=1 sat in any FOV

    The mount CYCLES through all pointing positions repeatedly for the entire
    night (dwell --> slew --> dwell --> slew --> ... until dawn).
    Each field will therefore be revisited multiple times per night.
    """
    observer = wgs84.latlon(ELGINFIELD_LAT, ELGINFIELD_LON,
                             elevation_m=ELGINFIELD_ELEV)

    # --- Night window ---
    t_night_start, t_night_end = get_night_window(observer, ts, year, month, day,
                                                   min_alt_deg)
    total_dark_min = (t_night_end.tt - t_night_start.tt) * 24 * 60
    if total_dark_min <= 0:
        total_dark_min = 480.0   # fallback 8 hours

    if verbose:
        print(f"\n{'='*70}")
        print(f"  {label}  ({year}-{month:02d}-{day:02d})")
        print(f"{'='*70}")
        ts_str = t_night_start.utc_strftime("%H:%M UTC")
        te_str = t_night_end.utc_strftime("%H:%M UTC")
        print(f"  Night window: {ts_str} -> {te_str}  ({total_dark_min:.0f} min)")

    # --- Pointing sequence ---
    array_width = n_cameras * (fov_deg + gap_deg) - gap_deg

    pointings = build_pointing_sequence(
        n_cameras, fov_deg, gap_deg,
        t_night_start, t_night_end,
        observer, ts, eph, min_alt_deg
    )

    n_pointings = len(pointings)
    # Estimate how many times we cycle through all pointings in the night
    cycle_time_min = n_pointings * dwell_min
    n_cycles_est   = total_dark_min / cycle_time_min if cycle_time_min > 0 else 1

    if verbose:
        print(f"  Unique pointings in one sweep  : {n_pointings}")
        print(f"  Cycle time (one full sweep)    : {cycle_time_min:.0f} min")
        print(f"  Estimated revisits per field   : {n_cycles_est:.1f}x")
        print(f"  Array width                    : {array_width:.1f}deg  "
              f"({n_cameras} cameras x {fov_deg}deg FOV)")

    # --- Per-dwell satellite check ---
    dwell_records  = []
    total_dwell_with_sats = 0.0
    all_detected_ids = set()

    t_cursor_tt = t_night_start.tt
    t_end_tt    = t_night_end.tt

    # Pre-compute RA/Dec for all satellites at 1-min resolution for the whole night
    n_steps = max(1, int(round(total_dark_min)))
    t_grid_tt = np.linspace(t_night_start.tt, t_night_end.tt, n_steps)
    t_grid = ts.tt_jd(t_grid_tt)

    if verbose:
        print(f"  Propagating {len(satellites)} GEO satellites x {n_steps} timesteps...")

    # Vectorised propagation:
    #   sat_radec[s] = (ra_deg_array, dec_deg_array)  -- geocentric, for FOV matching
    #   sat_alt_deg[s] = alt_degrees_array             -- topocentric, for visibility window
    sat_radec   = []
    sat_alt_deg = []
    for sat in satellites:
        try:
            astr = sat.at(t_grid)
            ra, dec, _ = astr.radec()
            sat_radec.append((ra._degrees, dec.degrees))
            # Topocentric elevation from Elginfield
            topo = (sat - observer).at(t_grid)
            alt_topo, _, _ = topo.altaz()
            sat_alt_deg.append(alt_topo.degrees)
        except Exception:
            sat_radec.append((None, None))
            sat_alt_deg.append(None)

    def tt_to_idx(tt_val):
        return int(np.clip(
            round((tt_val - t_night_start.tt) / (t_night_end.tt - t_night_start.tt)
                  * (n_steps - 1)),
            0, n_steps - 1
        ))

    # Track how many times each pointing position has been visited
    revisit_count  = [0] * n_pointings   # indexed by position in pointings list
    dwell_seq_idx  = 0                   # global sequential dwell counter
    # Per-satellite: total minutes actually imaged during the night
    sat_imaged_min = defaultdict(float)  # s_idx -> minutes in FOV

    # CYCLE through all pointings repeatedly until the night ends
    for pos_idx, (mount_ra, mount_dec) in enumerate(itertools.cycle(pointings)):
        if t_cursor_tt >= t_end_tt:
            break

        t_dwell_end_tt = t_cursor_tt + dwell_min / 1440.0
        g_idx_start = tt_to_idx(t_cursor_tt)
        g_idx_end   = tt_to_idx(min(t_dwell_end_tt, t_end_tt))

        local_pos = pos_idx % n_pointings   # position index within one cycle
        revisit   = revisit_count[local_pos]
        revisit_count[local_pos] += 1

        fov_list = camera_fovs(mount_ra, mount_dec, n_cameras, fov_deg, gap_deg)

        unique_this_dwell = set()
        cam_hits = [0] * n_cameras

        for s_idx in range(len(satellites)):
            ra_arr, dec_arr = sat_radec[s_idx]
            if ra_arr is None:
                continue
            for g in range(g_idx_start, g_idx_end + 1):
                cam_num = sat_in_fov(float(ra_arr[g]), float(dec_arr[g]), fov_list)
                if cam_num > 0:
                    unique_this_dwell.add(s_idx)
                    cam_hits[cam_num - 1] += 1
                    break

        actual_dwell_min = (min(t_dwell_end_tt, t_end_tt) - t_cursor_tt) * 1440.0
        if unique_this_dwell:
            total_dwell_with_sats += actual_dwell_min
        all_detected_ids.update(unique_this_dwell)
        # Accumulate per-satellite imaged time
        for s_idx in unique_this_dwell:
            sat_imaged_min[s_idx] += actual_dwell_min

        t_start_obj = ts.tt_jd(t_cursor_tt)
        dwell_records.append({
            "dwell_seq":    dwell_seq_idx,
            "pointing_pos": local_pos,          # 0-based index in pointings list
            "revisit":      revisit,             # 0 = first visit, 1 = second, ...
            "mount_ra":     mount_ra,
            "mount_dec":    mount_dec,
            "t_start_utc":  t_start_obj.utc_strftime("%H:%M"),
            "dwell_min":    actual_dwell_min,
            "n_unique":     len(unique_this_dwell),
            "sat_ids":      unique_this_dwell,
            "cam_hits":     cam_hits,
            "fov_list":     fov_list,
        })

        t_cursor_tt = t_dwell_end_tt
        dwell_seq_idx += 1

    pct_time = 100.0 * total_dwell_with_sats / total_dark_min if total_dark_min > 0 else 0.0
    actual_revisits = max(revisit_count) if revisit_count else 0

    # --- Per-satellite average coverage efficiency ---
    # For each detected satellite:
    #   max_vis_min = minutes it was above min_alt during the night
    #   ratio       = actual_imaged_min / max_vis_min
    # Then we average ratio over all detected satellites.
    min_per_step = total_dark_min / n_steps if n_steps > 1 else 1.0
    coverage_ratios = []
    for s_idx, imaged_min in sat_imaged_min.items():
        alt_arr = sat_alt_deg[s_idx]
        if alt_arr is None:
            continue
        max_vis_steps = int(np.sum(alt_arr >= min_alt_deg))
        max_vis_min   = max_vis_steps * min_per_step
        if max_vis_min > 0:
            coverage_ratios.append(min(1.0, imaged_min / max_vis_min))
    avg_sat_coverage_pct = 100.0 * np.mean(coverage_ratios) if coverage_ratios else 0.0

    if verbose:
        print(f"\n  Results:")
        print(f"    Unique GEO sats detected       : {len(all_detected_ids)}")
        print(f"    Total dwell intervals          : {dwell_seq_idx}")
        print(f"    Actual revisits per field      : {actual_revisits}x")
        print(f"    Avg sat coverage efficiency    : {avg_sat_coverage_pct:.1f}%")
        print(f"      (mean of: sat imaged time / sat visibility window, over all detected sats)")
        cam_totals = [0] * n_cameras
        for rec in dwell_records:
            for ci, h in enumerate(rec["cam_hits"]):
                cam_totals[ci] += h
        print(f"    Satellite detections per camera:")
        for ci, tot in enumerate(cam_totals):
            print(f"      Camera {ci+1}: {tot} detection-timesteps")

    return {
        "label":                label,
        "date":                 f"{year}-{month:02d}-{day:02d}",
        "total_dark_min":       total_dark_min,
        "unique_sats_total":    len(all_detected_ids),
        "pct_time_sats":        pct_time,
        "avg_sat_coverage_pct": avg_sat_coverage_pct,
        "actual_revisits":      actual_revisits,
        "n_pointings":          n_pointings,
        "cycle_time_min":       cycle_time_min,
        "dwell_records":        dwell_records,
        "n_cameras":            n_cameras,
        "fov_deg":              fov_deg,
        "gap_deg":              gap_deg,
        "dwell_min":            dwell_min,
    }


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_results(results, save=True):
    """Generate a multi-panel figure per simulation date."""
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        import matplotlib.cm as cm
    except ImportError:
        print("matplotlib not available – skipping plots.")
        return

    CAMERA_COLORS = ["#E63946", "#457B9D", "#2A9D8F", "#E9C46A",
                     "#F4A261", "#264653", "#A8DADC", "#8338EC"]

    for res in results:
        label  = res["label"]
        date_s = res["date"]
        recs   = res["dwell_records"]
        n_cam  = res["n_cameras"]
        fov    = res["fov_deg"]
        gap    = res["gap_deg"]

        fig = plt.figure(figsize=(18, 12))
        fig.patch.set_facecolor("#0d1117")
        title = (f"GEO Belt Pointing Simulation – {label} ({date_s})\n"
                 f"Elginfield Observatory  |  {n_cam} cameras  |  "
                 f"FOV {fov}°×{fov}°/cam  |  gap {gap}°")
        fig.suptitle(title, color="white", fontsize=13, y=0.98)

        ax_style = dict(facecolor="#161b22", labelcolor="white",
                        tick_params=dict(colors="white"))

        # ── Panel 1: Camera array layout (RA/Dec plane) ─────────────────────
        ax1 = fig.add_subplot(2, 3, 1)
        ax1.set_facecolor("#161b22")
        ax1.tick_params(colors="white")
        ax1.xaxis.label.set_color("white"); ax1.yaxis.label.set_color("white")
        for sp in ax1.spines.values(): sp.set_edgecolor("#444")

        # Draw a representative pointing mid-night
        mid_idx = len(recs) // 2 if recs else 0
        if recs:
            ex = recs[mid_idx]
            for ci, (cam_ra, cam_dec, half) in enumerate(ex["fov_list"]):
                col = CAMERA_COLORS[ci % len(CAMERA_COLORS)]
                rect = mpatches.FancyArrowPatch
                sq = mpatches.Rectangle(
                    (cam_ra - half, cam_dec - half), 2 * half, 2 * half,
                    linewidth=2, edgecolor=col, facecolor=col + "33",
                    label=f"Cam {ci+1}"
                )
                ax1.add_patch(sq)
                ax1.text(cam_ra, cam_dec, str(ci + 1), color="white",
                         ha="center", va="center", fontsize=11, fontweight="bold")
            # Horizon reference (dec=0 = equator / GEO belt)
            x0 = ex["fov_list"][0][0] - half - 1
            x1 = ex["fov_list"][-1][0] + half + 1
            ax1.axhline(0, color="#aaaaaa", linestyle="--", linewidth=1,
                        label="Celestial equator (GEO belt)")
            ax1.set_xlim(x0, x1)
            ax1.set_ylim(-fov * 2, fov * 2)
        ax1.set_xlabel("Right Ascension (°)", fontsize=9)
        ax1.set_ylabel("Declination (°)", fontsize=9)
        ax1.set_title("Camera Array Layout\n(representative mid-night pointing)",
                      color="white", fontsize=10)
        ax1.legend(loc="upper right", fontsize=7,
                   facecolor="#161b22", labelcolor="white", framealpha=0.8)

        # ── Panel 2: Pointing sequence arc in (RA, alt-elevation proxy) ─────
        ax2 = fig.add_subplot(2, 3, 2)
        ax2.set_facecolor("#161b22")
        ax2.tick_params(colors="white")
        ax2.xaxis.label.set_color("white"); ax2.yaxis.label.set_color("white")
        for sp in ax2.spines.values(): sp.set_edgecolor("#444")

        p_ras  = [r["mount_ra"] for r in recs]
        p_nsats = [r["n_unique"] for r in recs]
        if p_ras:
            scatter = ax2.scatter(p_ras, p_nsats,
                                  c=p_nsats, cmap="plasma", s=60,
                                  vmin=0, vmax=max(p_nsats) or 1, zorder=3)
            ax2.plot(p_ras, p_nsats, color="#444", linewidth=0.8, zorder=2)
            cb = plt.colorbar(scatter, ax=ax2)
            cb.set_label("# unique GEO sats", color="white", fontsize=8)
            cb.ax.yaxis.set_tick_params(color="white")
            plt.setp(plt.getp(cb.ax.axes, "yticklabels"), color="white")
        ax2.set_xlabel("Mount RA (°)", fontsize=9)
        ax2.set_ylabel("Unique sats in FOV", fontsize=9)
        ax2.set_title("Satellites per Pointing\n(along ecliptic track)",
                      color="white", fontsize=10)

        # ── Panel 3: Bar chart – sats per dwell over time ───────────────────
        ax3 = fig.add_subplot(2, 3, 3)
        ax3.set_facecolor("#161b22")
        ax3.tick_params(colors="white")
        ax3.xaxis.label.set_color("white"); ax3.yaxis.label.set_color("white")
        for sp in ax3.spines.values(): sp.set_edgecolor("#444")

        times_str = [r["t_start_utc"] for r in recs]
        nsats_arr = [r["n_unique"] for r in recs]
        x_pos = range(len(times_str))
        bar_cols = ["#2A9D8F" if n > 0 else "#444" for n in nsats_arr]
        ax3.bar(x_pos, nsats_arr, color=bar_cols, width=0.8)
        # Tick every ~15 positions
        tick_step = max(1, len(times_str) // 10)
        ax3.set_xticks(list(x_pos)[::tick_step])
        ax3.set_xticklabels(times_str[::tick_step], rotation=45, fontsize=7)
        ax3.set_xlabel("Dwell start (UTC)", fontsize=9)
        ax3.set_ylabel("# unique GEO sats", fontsize=9)
        ax3.set_title("GEO Satellites per Dwell\n(green = ≥1 detection)",
                      color="white", fontsize=10)

        # ── Panel 4: Per-camera detection histogram ──────────────────────────
        ax4 = fig.add_subplot(2, 3, 4)
        ax4.set_facecolor("#161b22")
        ax4.tick_params(colors="white")
        ax4.xaxis.label.set_color("white"); ax4.yaxis.label.set_color("white")
        for sp in ax4.spines.values(): sp.set_edgecolor("#444")

        cam_totals = [0] * n_cam
        for rec in recs:
            for ci, h in enumerate(rec["cam_hits"]):
                cam_totals[ci] += h
        cam_labels = [f"Cam {i+1}" for i in range(n_cam)]
        cam_cols   = [CAMERA_COLORS[i % len(CAMERA_COLORS)] for i in range(n_cam)]
        ax4.bar(cam_labels, cam_totals, color=cam_cols)
        ax4.set_xlabel("Camera", fontsize=9)
        ax4.set_ylabel("Total detection-timesteps", fontsize=9)
        ax4.set_title("Detections per Camera\n(timesteps satellite in FOV)",
                      color="white", fontsize=10)
        for tick in ax4.get_xticklabels():
            tick.set_color("white")

        # ── Panel 5: Histogram of sats per dwell ────────────────────────────
        ax5 = fig.add_subplot(2, 3, 5)
        ax5.set_facecolor("#161b22")
        ax5.tick_params(colors="white")
        ax5.xaxis.label.set_color("white"); ax5.yaxis.label.set_color("white")
        for sp in ax5.spines.values(): sp.set_edgecolor("#444")

        if nsats_arr and max(nsats_arr) > 0:
            ax5.hist(nsats_arr, bins=min(30, max(nsats_arr) + 1),
                     color="#457B9D", edgecolor="#aaa", linewidth=0.5)
        ax5.set_xlabel("# unique GEO sats per dwell", fontsize=9)
        ax5.set_ylabel("# of dwell intervals", fontsize=9)
        ax5.set_title("Distribution of Detections\nper Dwell", color="white", fontsize=10)

        # ── Panel 6: Summary stats text ──────────────────────────────────────
        ax6 = fig.add_subplot(2, 3, 6)
        ax6.set_facecolor("#161b22")
        ax6.axis("off")
        summary_lines = [
            f"Date:            {date_s}",
            f"Location:        Elginfield, ON",
            f"# Cameras:       {n_cam}",
            f"FOV/cam:         {fov}° × {fov}°",
            f"Gap:             {gap}°",
            f"Array width:     {n_cam*(fov+gap)-gap:.1f}°",
            f"Dwell time:      {res['dwell_min']} min",
            f"Dark night:      {res['total_dark_min']:.0f} min",
            f"# Pointings:     {len(recs)}",
            f"Unique GEO sats: {res['unique_sats_total']}",
            f"% time w/ sats:  {res['pct_time_sats']:.1f}%",
        ]
        ax6.text(0.05, 0.95, "\n".join(summary_lines),
                 transform=ax6.transAxes, color="white",
                 fontsize=11, va="top", fontfamily="monospace",
                 bbox=dict(facecolor="#1e2531", edgecolor="#444", boxstyle="round,pad=0.5"))
        ax6.set_title("Summary", color="white", fontsize=10)

        plt.tight_layout(rect=[0, 0, 1, 0.96])

        if save:
            safe_label = label.replace(" ", "_").lower()
            fname = f"geo_pointing_plot_{safe_label}.png"
            plt.savefig(fname, dpi=150, facecolor="#0d1117")
            print(f"  Saved: {fname}")
        try:
            plt.show(block=False)
            plt.pause(0.5)
        except Exception:
            pass
        plt.close()


# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

def print_summary_table(results):
    w = 105
    print("\n" + "=" * w)
    print(f"{'Date':<20} {'Night (min)':<13} {'Unique GEO sats':<17} "
          f"{'Revisits/field':<16} {'Avg sat coverage %'}")
    print(f"                                                                   "
          f"                 (imaged / visibility window)")
    print("-" * w)
    for r in results:
        revisits = r.get("actual_revisits", "?")
        print(f"{r['label']:<20} {r['total_dark_min']:<13.0f} "
              f"{r['unique_sats_total']:<17} "
              f"{str(revisits) + 'x':<16} "
              f"{r['avg_sat_coverage_pct']:.1f}%")
    print("=" * w)



# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Simulate GEO belt pointing pattern for a multi-camera tracking mount.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-n", "--n-cameras", type=int,   default=4,
                        help="Number of cameras on the mount.")
    parser.add_argument("-f", "--fov",       type=float, default=5.0,
                        help="Per-camera square FOV width (degrees).")
    parser.add_argument("-g", "--gap",       type=float, default=0.0,
                        help="Gap between adjacent cameras (degrees, 0=adjacent).")
    parser.add_argument("-d", "--dwell",     type=float, default=10.0,
                        help="Sidereal tracking dwell time (minutes).")
    parser.add_argument("--min-alt",   type=float, default=10.0, dest="min_alt",
                        help="Minimum pointing altitude (degrees).")
    parser.add_argument("--year",      type=int,   default=2026,
                        help="Year for simulation.")
    parser.add_argument("--tle-file",  type=str,   default=DEFAULT_TLE_FILE,
                        help="Local GEO TLE file path.")
    parser.add_argument("--refresh",   action="store_true",
                        help="Re-download GEO TLEs from Celestrak before running.")
    parser.add_argument("--plot",      action="store_true",
                        help="Save PNG plots for each date.")
    parser.add_argument("--dates",     type=str,   default="all",
                        help=("Which dates to simulate: 'all', 'equinox', 'solstice', "
                              "or comma-separated YYYY-MM-DD values."))
    parser.add_argument("--no-summary", action="store_true",
                        help="Suppress per-night verbose output.")
    return parser.parse_args()


def select_dates(dates_arg, year):
    """Parse --dates argument into list of (label, month, day) tuples."""
    arg_lower = dates_arg.lower().strip()
    if arg_lower == "all":
        return [(lbl, m, d) for lbl, m, d in SIM_DATES]
    if arg_lower == "equinox":
        return [(lbl, m, d) for lbl, m, d in SIM_DATES if "Equinox" in lbl]
    if arg_lower == "solstice":
        return [(lbl, m, d) for lbl, m, d in SIM_DATES if "Solstice" in lbl]
    # Check if the string matches one of our standard labels (case-insensitive)
    label_matches = [(lbl, m, d) for lbl, m, d in SIM_DATES
                     if lbl.lower() == arg_lower]
    if label_matches:
        return label_matches
    # Try comma-separated list of YYYY-MM-DD or label names
    chosen = []
    for part in dates_arg.split(","):
        part = part.strip()
        # Check label names first
        lbl_match = [(lbl, m, d) for lbl, m, d in SIM_DATES
                     if lbl.lower() == part.lower()]
        if lbl_match:
            chosen.extend(lbl_match)
            continue
        # Try YYYY-MM-DD
        try:
            dt = datetime.strptime(part, "%Y-%m-%d")
            chosen.append((part, dt.month, dt.day))
        except ValueError:
            print(f"  Warning: could not parse date '{part}' – skipping.")
    return chosen


def main():
    args = parse_args()

    # ── Load ephemeris ───────────────────────────────────────────────────────
    ts = load.timescale()
    try:
        eph = load("de421.bsp")
    except Exception as e:
        sys.exit(f"Error loading DE421 ephemeris: {e}\n"
                 "Download with:  python -c \"from skyfield.api import load; load('de421.bsp')\"")

    # ── TLE loading ──────────────────────────────────────────────────────────
    if args.refresh:
        refresh_tles(args.tle_file)

    satellites = load_tle_file(args.tle_file, ts)
    if not satellites:
        print(f"No TLEs in '{args.tle_file}'. Downloading from Celestrak…")
        refresh_tles(args.tle_file)
        satellites = load_tle_file(args.tle_file, ts)

    if not satellites:
        sys.exit("Error: Could not load any GEO TLEs. Check network or TLE file path.")

    print(f"\nLoaded {len(satellites)} GEO satellites from '{args.tle_file}'")
    print(f"Configuration: {args.n_cameras} cameras  |  FOV {args.fov}°×{args.fov}°/cam  "
          f"|  gap {args.gap}°  |  dwell {args.dwell} min\n")

    # ── Dates ────────────────────────────────────────────────────────────────
    chosen_dates = select_dates(args.dates, args.year)
    if not chosen_dates:
        sys.exit("No valid simulation dates selected.")

    # ── Run simulations ──────────────────────────────────────────────────────
    all_results = []
    for label, month, day in chosen_dates:
        res = simulate_night(
            label=label, year=args.year, month=month, day=day,
            satellites=satellites, ts=ts, eph=eph,
            n_cameras=args.n_cameras, fov_deg=args.fov, gap_deg=args.gap,
            dwell_min=args.dwell, min_alt_deg=args.min_alt,
            verbose=not args.no_summary
        )
        all_results.append(res)

    # ── Summary table ────────────────────────────────────────────────────────
    print_summary_table(all_results)

    # ── Plots ────────────────────────────────────────────────────────────────
    if args.plot:
        print("\nGenerating plots…")
        plot_results(all_results, save=True)

    print("\nDone.")


if __name__ == "__main__":
    main()
