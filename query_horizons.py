import argparse
import sys
import datetime
import os
import requests
import math
from astropy.coordinates import SkyCoord
import astropy.units as u

# Default Elginfield location
DEFAULT_LAT = 43.1925
DEFAULT_LON = -81.3158
DEFAULT_ELEV_KM = 0.325

def parse_time(t_str):
    """Parse YYYYMMDD_HHmmss to datetime, or return None."""
    try:
        return datetime.datetime.strptime(t_str, "%Y%m%d_%H%M%S")
    except ValueError:
        raise argparse.ArgumentTypeError(f"Time must be in YYYYMMDD_HHmmss format, got {t_str}")

def get_moon_limb_sep(target_ra, target_dec, moon_ra, moon_dec, moon_delta_km):
    """Calculate the separation from the Moon's physical limb in degrees."""
    try:
        sc = SkyCoord(target_ra, target_dec, unit=(u.hourangle, u.deg))
        mc = SkyCoord(moon_ra, moon_dec, unit=(u.hourangle, u.deg))
        sep = sc.separation(mc).degree
        # Moon angular radius (1737.4 km physical radius)
        moon_ang_rad = math.degrees(math.asin(1737.4 / float(moon_delta_km)))
        return sep - moon_ang_rad
    except Exception as e:
        print(f"Exception parsing SkyCoord: {e}")
        return None

def estimate_magnitude(delta_km, phase_angle_deg, area_m2=20.0, albedo=0.2):
    """Estimate visual apparent magnitude using a diffuse sphere model."""
    import math
    v_sun = -26.74
    alpha = math.radians(phase_angle_deg)
    
    # Phase function
    phi_alpha = (math.sin(alpha) + (math.pi - alpha) * math.cos(alpha)) / math.pi
    if phi_alpha <= 0:
        return None
        
    d_m = delta_km * 1000.0
    val = (area_m2 * albedo * phi_alpha) / (d_m**2)
    
    if val <= 0:
        return None
        
    return v_sun - 2.5 * math.log10(val)

def apply_airmass_extinction(unc_mag, elev_deg, k=0.2):
    """Calculate airmass correction using the Kasten-Young (1989) formula."""
    import math
    if elev_deg <= 0:
        return None
    # Kasten-Young formula for geometric airmass
    h_rad = math.radians(elev_deg)
    # Protection against zero/negative fractional power bases
    denominator = math.sin(h_rad) + 0.50572 * math.pow(6.07995 + elev_deg, -1.6364)
    if denominator <= 0:
        return None
        
    airmass = 1.0 / denominator
    return unc_mag + (k * airmass)

def main():
    now = datetime.datetime.utcnow()
    default_start = now.strftime("%Y%m%d_%H%M%S")
    default_stop = (now + datetime.timedelta(days=1)).strftime("%Y%m%d_%H%M%S")

    parser = argparse.ArgumentParser(description="Query JPL Horizons for ephemerides and export to CSV.")
    parser.add_argument("--obj", type=str, default="-1024", help="Object ID (default: -1024)")
    parser.add_argument("--loc", type=str, default="Elginfield", help="Location name (default: Elginfield)")
    parser.add_argument("--latlon", type=float, nargs=2, metavar=('LAT', 'LON'), help="Override location by providing latitude and longitude in degrees")
    parser.add_argument("--start", type=str, default=default_start, help=f"Start time in YYYYMMDD_HHmmss format (default: {default_start})")
    parser.add_argument("--stop", type=str, default=default_stop, help=f"Stop time in YYYYMMDD_HHmmss format (default: {default_stop})")
    parser.add_argument("--step", type=str, default="30", help="Time step interval in minutes (default: 30). Optional: provide '1h' or '1d' for hours/days.")
    parser.add_argument("--out", type=str, default=None, help="Output CSV filename (default: auto-generated based on target and time)")
    parser.add_argument("--tle", type=str, default=None, metavar="TLE_FILE_OR_LINES",
                        help="Path to a TLE file (name line + line 1 + line 2, or just line 1 + line 2) "
                             "OR two raw TLE lines separated by a pipe '|'. "
                             "If provided, skyfield is used instead of JPL Horizons. "
                             "'_TLE' is appended to the auto-generated CSV filename.")
    
    args = parser.parse_args()

    if args.tle:
        run_tle_mode(args)
        return
    
    # Process step interval to default to minutes if numeric alone is passed
    step_val = args.step
    if step_val.isdigit():
        step_val = f"{step_val}m"

    start_dt = parse_time(args.start)
    stop_dt = parse_time(args.stop)

    # East Longitude is positive, West is negative.
    if args.latlon:
        lat, lon = args.latlon
        site_coord = f"{lon},{lat},0.0" # Default altitude 0 km if latlon provided
        center = "coord@399"
    elif args.loc.lower() == "elginfield":
        site_coord = f"{DEFAULT_LON},{DEFAULT_LAT},{DEFAULT_ELEV_KM}"
        center = "coord@399"
    else:
        center = args.loc
        site_coord = None

    params = {
        "format": "text",
        "COMMAND": f"'{args.obj}'",
        "OBJ_DATA": "NO",
        "MAKE_EPHEM": "YES",
        "EPHEM_TYPE": "OBSERVER",
        "CENTER": f"'{center}'",
        "START_TIME": f"'{start_dt.strftime('%Y-%m-%d %H:%M:%S')}'",
        "STOP_TIME": f"'{stop_dt.strftime('%Y-%m-%d %H:%M:%S')}'",
        "STEP_SIZE": f"'{step_val}'",
        "CSV_FORMAT": "YES",
        # Request commonly used quantities: 1 (Astrometric RA/Dec), 2 (Apparent RA/Dec),
        # 3 (Rates of RA/Dec), 4 (Apparent Az/El), 9 (Visual Mag), 20 (Observer range/rate), 
        # 24 (Sun-Target-Obs angle)
        "QUANTITIES": "'1,2,3,4,9,20,24'",
        "RANGE_UNITS": "'km'",
    }
    if site_coord:
        params["SITE_COORD"] = f"'{site_coord}'"

    print(f"Querying JPL Horizons API for object '{args.obj}'...")
    url = "https://ssd.jpl.nasa.gov/api/horizons.api"
    try:
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        sys.exit(f"Error querying Horizons: {e}")

    text = response.text

    if "$$SOE" not in text:
        print("Error: Horizons response did not contain ephemeris data ($$SOE not found).")
        print("Full response below:")
        print(text)
        sys.exit(1)

    print("Querying JPL Horizons API for Moon (301) to calculate limb separation...")
    moon_params = params.copy()
    moon_params["COMMAND"] = "'301'"
    try:
        moon_response = requests.get(url, params=moon_params, timeout=30)
        moon_response.raise_for_status()
        moon_text = moon_response.text
        if "$$SOE" not in moon_text:
            moon_text = None
    except Exception as e:
        print(f"Warning: Could not fetch Moon ephemeris: {e}")
        moon_text = None

    def extract_soe_block(res_text):
        if not res_text: return [], ""
        lns = res_text.splitlines()
        s_idx, e_idx = -1, -1
        for i, l in enumerate(lns):
            if l.strip() == "$$SOE": s_idx = i
            elif l.strip() == "$$EOE": e_idx = i
        if s_idx != -1 and e_idx != -1:
            header = lns[s_idx - 2].strip() if s_idx >= 2 else ""
            return lns[s_idx+1:e_idx], header
        return [], ""

    target_lines, headers_line = extract_soe_block(text)
    moon_lines, moon_headers = extract_soe_block(moon_text)
    
    # Parse Moon Headers
    moon_ra_idx, moon_dec_idx, moon_delta_idx = -1, -1, -1
    if moon_headers:
        m_cols = [c.strip() for c in moon_headers.split(',')]
        if 'S-brt' in m_cols:
            m_cols.remove('S-brt')
        moon_ra_idx = next((i for i, c in enumerate(m_cols) if 'R.A.' in c and 'a-app' in c), -1)
        moon_dec_idx = next((i for i, c in enumerate(m_cols) if 'DEC' in c and 'a-app' in c), -1)
        moon_delta_idx = next((i for i, c in enumerate(m_cols) if 'delta' in c), -1)

    csv_data = []
    sbrt_idx = -1
    apmag_idx = -1
    delta_idx = -1
    sto_idx = -1
    elev_idx = -1
    
    if headers_line:
        cols = [c.strip() for c in headers_line.split(',')]
        if 'S-brt' in cols:
            sbrt_idx = cols.index('S-brt')
            cols.pop(sbrt_idx)
            
        # Add units and fill empty status header names
        for j, col in enumerate(cols):
            if col == '' and j == 1:
                cols[j] = 'Solar_Presence'
            elif col == '' and j == 2:
                cols[j] = 'Lunar_Presence'
            elif col == 'delta':
                cols[j] = 'delta_(km)'
            elif col == 'deldot':
                cols[j] = 'deldot_(km/s)'
            elif col == 'S-T-O':
                cols[j] = 'S-T-O_(deg)'
                
        apmag_idx = cols.index('APmag') if 'APmag' in cols else -1
        delta_idx = cols.index('delta_(km)') if 'delta_(km)' in cols else -1
        sto_idx = cols.index('S-T-O_(deg)') if 'S-T-O_(deg)' in cols else -1
        elev_idx = cols.index('Elev_(a-app)') if 'Elev_(a-app)' in cols else -1
        target_ra_idx = next((i for i, c in enumerate(cols) if 'R.A.' in c and 'a-app' in c), -1)
        target_dec_idx = next((i for i, c in enumerate(cols) if 'DEC' in c and 'a-app' in c), -1)
        
        # Process RA/Dec rate headers to rename them
        ra_rate_idx = cols.index('dRA*cosD') if 'dRA*cosD' in cols else -1
        dec_rate_idx = cols.index('d(DEC)/dt') if 'd(DEC)/dt' in cols else -1
        
        if ra_rate_idx != -1:
            cols[ra_rate_idx] = 'dRA*cosD_(arcsec/s)'
        if dec_rate_idx != -1:
            cols[dec_rate_idx] = 'd(DEC)/dt_(arcsec/s)'
        
        # Inject the new column Header right after APmag if it exists, otherwise at the end
        if apmag_idx != -1:
            cols.insert(apmag_idx + 1, 'APmag_corr')
        else:
            cols.append('APmag_corr')

        if cols and cols[-1] == '':
            cols[-1] = 'Moon_Limb_Sep_(deg)'
            cols.append('')
        else:
            cols.append('Moon_Limb_Sep_(deg)')
            
        headers_line = ",".join(cols)
        csv_data.append(headers_line)
        
    for i in range(len(target_lines)):
        line_str = target_lines[i].strip()
        if not line_str:
            continue
            
        parts = [p.strip() for p in line_str.split(',')]
        
        if sbrt_idx != -1 and len(parts) > sbrt_idx:
            parts.pop(sbrt_idx)
            
        # Estimate magnitude for Artemis II specifically
        if args.obj == '-1024' and apmag_idx != -1 and delta_idx != -1 and sto_idx != -1:
            if len(parts) > max(apmag_idx, delta_idx, sto_idx):
                if parts[apmag_idx] == 'n.a.':
                    try:
                        delta = float(parts[delta_idx])
                        sto = float(parts[sto_idx])
                        est_mag = estimate_magnitude(delta, sto)
                        if est_mag is not None:
                            parts[apmag_idx] = f"{est_mag:.2f}"
                    except ValueError:
                        pass
        
        # Calculate Airmass correction globally
        corr_val = 'n.a.'
        if apmag_idx != -1 and elev_idx != -1 and len(parts) > max(apmag_idx, elev_idx):
            try:
                base_mag = float(parts[apmag_idx])
                elev = float(parts[elev_idx])
                corr_mag = apply_airmass_extinction(base_mag, elev)
                if corr_mag is not None:
                    corr_val = f"{corr_mag:.2f}"
            except ValueError:
                pass
                
        # Convert Apparent Rates of Motion from arcsec/hr to arcsec/s
        if ra_rate_idx != -1 and len(parts) > ra_rate_idx and parts[ra_rate_idx] != 'n.a.':
            try:
                parts[ra_rate_idx] = f"{float(parts[ra_rate_idx]) / 3600.0:.6f}"
            except ValueError:
                pass
        if dec_rate_idx != -1 and len(parts) > dec_rate_idx and parts[dec_rate_idx] != 'n.a.':
            try:
                parts[dec_rate_idx] = f"{float(parts[dec_rate_idx]) / 3600.0:.6f}"
            except ValueError:
                pass
                
        # Insert correction column
        if apmag_idx != -1:
            parts.insert(apmag_idx + 1, corr_val)
        else:
            parts.append(corr_val)
            
        # Calculate Moon Limb Separation
        limb_val = 'n.a.'
        if moon_lines and i < len(moon_lines):
            # Same treatment for trailing commas
            m_line_str = moon_lines[i].strip()
            if m_line_str:
                m_parts = [p.strip() for p in m_line_str.split(',')]
                # S-brt was removed from m_cols natively, so we just remove from m_parts without relying on m_cols shifting
                # We need to find S-brt in the RAW header since m_parts contains the RAW columns still
                raw_m_cols = [c.strip() for c in moon_headers.split(',')]
                m_s_idx = raw_m_cols.index('S-brt') if 'S-brt' in raw_m_cols else -1
                if m_s_idx != -1 and len(m_parts) > m_s_idx:
                    m_parts.pop(m_s_idx)

                if target_ra_idx != -1 and target_dec_idx != -1 and moon_ra_idx != -1 and moon_dec_idx != -1 and moon_delta_idx != -1:
                    if len(parts) > max(target_ra_idx, target_dec_idx) and len(m_parts) > max(moon_ra_idx, moon_dec_idx, moon_delta_idx):
                        tra = parts[target_ra_idx]
                        tdec = parts[target_dec_idx]
                        mra = m_parts[moon_ra_idx]
                        mdec = m_parts[moon_dec_idx]
                        mdelt = m_parts[moon_delta_idx]
                        sep = get_moon_limb_sep(tra, tdec, mra, mdec, mdelt)
                        if sep is not None:
                            limb_val = f"{sep:.3f}"
        
        if parts and parts[-1] == '':
            parts[-1] = limb_val
            parts.append('')
        else:
            parts.append(limb_val)
            
        csv_data.append(",".join(parts))

    csv_output = "\n".join(csv_data)

    print("\n--- Ephemerides ---")
    print(csv_output)
    print("-------------------\n")

    if args.out:
        out_file = args.out
    else:
        safe_obj = args.obj.replace(' ', '_').replace('/', '_').replace('=', '_')
        out_file = f"ephemeris_{safe_obj}_{args.start}-{args.stop}.csv"

    try:
        with open(out_file, "w", encoding="utf-8") as f:
            f.write(csv_output + "\n")
        print(f"Successfully saved {len(csv_data)-1} ephemeris records to {out_file}")
    except OSError as e:
        print(f"Error saving to file {out_file}: {e}")


# ---------------------------------------------------------------------------
# TLE mode helpers
# ---------------------------------------------------------------------------

def _parse_step_to_seconds(step_str):
    """Convert step string (e.g. '30', '30m', '1h', '1d') to integer seconds."""
    s = step_str.strip()
    if s.endswith('d'):
        return int(s[:-1]) * 86400
    elif s.endswith('h'):
        return int(s[:-1]) * 3600
    elif s.endswith('m'):
        return int(s[:-1]) * 60
    else:
        return int(s) * 60  # bare number treated as minutes


def _load_tle_lines(tle_arg):
    """Return (name, line1, line2) from a file path or a pipe-separated pair of TLE lines."""
    if '|' in tle_arg:
        # Raw lines passed on command line separated by |
        parts = [p.strip() for p in tle_arg.split('|')]
        if len(parts) == 2:
            return ('OBJECT', parts[0], parts[1])
        elif len(parts) == 3:
            return (parts[0], parts[1], parts[2])
        else:
            sys.exit("Error: --tle pipe syntax expects 'LINE1|LINE2' or 'NAME|LINE1|LINE2'.")
    else:
        path = tle_arg.strip()
        if not os.path.isfile(path):
            sys.exit(f"Error: TLE file not found: {path}")
        with open(path, 'r') as f:
            raw = [l.rstrip() for l in f if l.strip()]
        if len(raw) == 2:
            return ('OBJECT', raw[0], raw[1])
        elif len(raw) >= 3:
            return (raw[0].strip(), raw[1], raw[2])
        else:
            sys.exit("Error: TLE file must contain at least two non-empty lines (line 1 and line 2).")


def run_tle_mode(args):
    """Compute ephemeris from a TLE using skyfield and write the same CSV format."""
    try:
        from skyfield.api import load, EarthSatellite, wgs84
        from skyfield.positionlib import Geocentric
    except ImportError:
        sys.exit("Error: skyfield is required for TLE mode. Install it with: pip install skyfield")

    tle_name, line1, line2 = _load_tle_lines(args.tle)
    print(f"Using TLE for: {tle_name}")
    print(f"  Line 1: {line1}")
    print(f"  Line 2: {line2}")

    # Observer location
    if args.latlon:
        obs_lat, obs_lon = args.latlon
        obs_elev_km = 0.0
    else:
        obs_lat = DEFAULT_LAT
        obs_lon = DEFAULT_LON
        obs_elev_km = DEFAULT_ELEV_KM

    ts = load.timescale()
    satellite = EarthSatellite(line1, line2, tle_name, ts)
    observer = wgs84.latlon(obs_lat, obs_lon, elevation_m=obs_elev_km * 1000.0)

    start_dt = parse_time(args.start)
    stop_dt  = parse_time(args.stop)
    step_str = args.step
    if step_str.isdigit():
        step_str = f"{step_str}m"
    step_sec = _parse_step_to_seconds(step_str)

    # Build list of UTC datetimes
    times_dt = []
    current = start_dt
    while current <= stop_dt:
        times_dt.append(current)
        current += datetime.timedelta(seconds=step_sec)

    if not times_dt:
        sys.exit("Error: No time steps generated. Check --start, --stop, and --step.")

    # Also load Moon ephemeris for limb separation
    try:
        eph = load('de421.bsp')
        earth = eph['earth']
        moon_body = eph['moon']
        have_moon = True
    except Exception as e:
        print(f"Warning: Could not load DE421 for Moon limb separation: {e}")
        have_moon = False

    # CSV header (mirrors Horizons output as closely as possible)
    header = ("Date__(UT)__HR:MN,Solar_Presence,Lunar_Presence,"
              "R.A._(ICRF),DEC_(ICRF),"
              "R.A._(a-app),DEC_(a-app),"
              "dRA*cosD_(arcsec/s),d(DEC)/dt_(arcsec/s),"
              "Azi_(a-app),Elev_(a-app),"
              "APmag,APmag_corr,"
              "delta_(km),deldot_(km/s),"
              "S-T-O_(deg),"
              "Moon_Limb_Sep_(deg)")

    csv_data = [header]

    DELTA_T = 1.0  # seconds, used for numerical rate derivatives

    for dt_utc in times_dt:
        t  = ts.from_datetime(dt_utc.replace(tzinfo=datetime.timezone.utc))
        dt_prev = dt_utc - datetime.timedelta(seconds=DELTA_T)
        dt_next = dt_utc + datetime.timedelta(seconds=DELTA_T)
        t_prev = ts.from_datetime(dt_prev.replace(tzinfo=datetime.timezone.utc))
        t_next = ts.from_datetime(dt_next.replace(tzinfo=datetime.timezone.utc))

        diff       = satellite - observer
        topocentric = diff.at(t)
        ra, dec, distance = topocentric.radec()          # ICRF / astrometric
        ra_app, dec_app, _ = topocentric.radec('date')   # apparent (of-date)
        alt, az, dist_km_obj = topocentric.altaz()

        delta_km = distance.km

        # Numerical range-rate (km/s)
        dist_prev = (satellite - observer).at(t_prev).radec()[2].km
        dist_next = (satellite - observer).at(t_next).radec()[2].km
        deldot = (dist_next - dist_prev) / (2.0 * DELTA_T)

        # Numerical RA/Dec rates (arcsec/s) – central difference
        ra_prev, dec_prev, _ = (satellite - observer).at(t_prev).radec('date')
        ra_next, dec_next, _ = (satellite - observer).at(t_next).radec('date')
        dec_mid_rad = math.radians(dec_app.degrees)
        dra_arcsec_s  = ((ra_next.hours - ra_prev.hours) * 15.0 * 3600.0
                         * math.cos(dec_mid_rad)) / (2.0 * DELTA_T)
        ddec_arcsec_s = ((dec_next.degrees - dec_prev.degrees) * 3600.0) / (2.0 * DELTA_T)

        elev_deg = alt.degrees

        # Airmass-corrected magnitude (APmag unavailable from TLE)
        apmag_str = 'n.a.'
        apmag_corr_str = 'n.a.'

        # S-T-O angle (Sun-Target-Observer) — unavailable without planetary ephemeris
        sto_str = 'n.a.'

        # Moon limb separation
        limb_val = 'n.a.'
        if have_moon:
            try:
                moon_topo = (earth + observer).at(t).observe(moon_body).apparent()
                mra, mdec, mdist = moon_topo.radec('date')
                moon_delta_km = mdist.km
                sat_ra_str  = _ra_hours_to_hms(ra_app.hours)
                sat_dec_str = f"{dec_app.degrees:.4f}"
                mra_str     = _ra_hours_to_hms(mra.hours)
                mdec_str    = f"{mdec.degrees:.4f}"
                limb_val_f  = get_moon_limb_sep(sat_ra_str, sat_dec_str, mra_str, mdec_str, moon_delta_km)
                if limb_val_f is not None:
                    limb_val = f"{limb_val_f:.3f}"
            except Exception:
                pass  # leave limb_val as n.a.

        date_str = dt_utc.strftime("%Y-%b-%d %H:%M")

        # Format RA as HH MM SS.ff
        ra_icrf_str  = _ra_hours_to_hms(ra.hours)
        dec_icrf_str = f"{dec.degrees:+.4f}"
        ra_app_str   = _ra_hours_to_hms(ra_app.hours)
        dec_app_str  = f"{dec_app.degrees:+.4f}"

        row = (
            f"{date_str},,"
            f"{ra_icrf_str},{dec_icrf_str},"
            f"{ra_app_str},{dec_app_str},"
            f"{dra_arcsec_s:.6f},{ddec_arcsec_s:.6f},"
            f"{az.degrees:.4f},{elev_deg:.4f},"
            f"{apmag_str},{apmag_corr_str},"
            f"{delta_km:.3f},{deldot:.6f},"
            f"{sto_str},"
            f"{limb_val}"
        )
        csv_data.append(row)

    csv_output = "\n".join(csv_data)

    print("\n--- Ephemerides (TLE mode) ---")
    print(csv_output)
    print("-----------------------------\n")

    if args.out:
        out_file = args.out
    else:
        safe_obj = args.obj.replace(' ', '_').replace('/', '_').replace('=', '_')
        out_file = f"ephemeris_{safe_obj}_{args.start}-{args.stop}_TLE.csv"

    try:
        with open(out_file, "w", encoding="utf-8") as f:
            f.write(csv_output + "\n")
        print(f"Successfully saved {len(csv_data)-1} ephemeris records to {out_file}")
    except OSError as e:
        print(f"Error saving to file {out_file}: {e}")


def _ra_hours_to_hms(hours):
    """Convert decimal hours to 'HH MM SS.ff' string matching Horizons format."""
    h = int(hours)
    remaining = (hours - h) * 60.0
    m = int(remaining)
    s = (remaining - m) * 60.0
    return f"{h:02d} {m:02d} {s:05.2f}"


if __name__ == "__main__":
    main()
