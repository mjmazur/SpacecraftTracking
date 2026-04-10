import requests
import argparse
from datetime import datetime, timedelta
import sys
import re

# Common satellites dictionary map (Name to Horizons ID)
# DES= indicates the NORAD ID. 
# Artemis P1 / P2 (THEMIS B/C) have Horizons IDs -117 and -118 respectively.
COMMON_SATELLITES = {
    "ISS": "DES=25544",
    "HST": "DES=20580",
    "CSS": "DES=39392",
    "TIANGONG": "DES=39392",
    "ARTEMIS-P1": "-117", # THEMIS-B / ARTEMIS P1
    "ARTEMIS-P2": "-118", # THEMIS-C / ARTEMIS P2
}

def resolve_target(target_query):
    """
    Resolves a user target query to a Horizons target parameter.
    If it's in the common list, return that. Otherwise, try to search for the Artemis string
    or assume the user provided a valid ID/Name.
    """
    key = target_query.upper().strip()
    if key in COMMON_SATELLITES:
        return COMMON_SATELLITES[key]
    return target_query

def get_horizons_ephemeris(target, location, start_time, stop_time, step_size):
    """
    Connects to JPL Horizons API and retrieves ephemeris data.
    """
    url = "https://ssd.jpl.nasa.gov/api/horizons.api"
    
    if location.lower() == "elginfield":
        location_param = "coord@399"
        site_coord = "-81.3158,43.1925,0.325"
    else:
        location_param = location

    # We request QUANTITIES='4,9'
    # 4: Apparent AZ & EL
    # 9: Visual magnitude & Surface Brightness (often required to check illumination/eclipse for satellites)
    params = {
        "format": "text",
        "COMMAND": target,
        "OBJ_DATA": "NO",
        "MAKE_EPHEM": "YES",
        "EPHEM_TYPE": "OBSERVER",
        "CENTER": location_param,
        "START_TIME": start_time,
        "STOP_TIME": stop_time,
        "STEP_SIZE": step_size,
        "QUANTITIES": "4,9",
        "CSV_FORMAT": "YES"
    }
    
    if location.lower() == "elginfield":
        params["SITE_COORD"] = site_coord

    try:
        response = requests.get(url, params=params)
        response.raise_for_status()
        return response.text
    except Exception as e:
        return f"Error connecting to JPL Horizons: {e}"

def search_artemis_satellites():
    """Dynamically returns known Artemis spacecraft from Horizons."""
    url = "https://ssd.jpl.nasa.gov/api/horizons_lookup.api"
    artemis_targets = []
    
    try:
        response = requests.get(url, params={"sstr": "Artemis", "format": "json"})
        response.raise_for_status()
        data = response.json()
        
        if "result" in data:
            for item in data["result"]:
                # Only include spacecraft
                if item.get("type", "").lower() == "spacecraft":
                    artemis_targets.append((item.get("name", "Artemis Satellite"), item.get("spkid")))
                    
    except Exception as e:
        print(f"Error fetching Artemis satellites dynamically: {e}")
        
    # Always include the THEMIS B/C ARTEMIS probes manually if not returned
    artemis_targets.extend([("ARTEMIS-P1 (THEMIS-B)", "-117"), ("ARTEMIS-P2 (THEMIS-C)", "-118")])
    return list(set(artemis_targets))  # remove potential duplicates

def get_all_spacecraft():
    """Fetches the list of all spacecraft from Horizons."""
    url = "https://ssd.jpl.nasa.gov/api/horizons.api"
    params = {
        "format": "text",
        "COMMAND": "SPACECRAFT"
    }
    targets = []
    try:
        response = requests.get(url, params=params)
        response.raise_for_status()
        text = response.text
        
        # Parse the block between the dashes and the end summary
        in_table = False
        for line in text.splitlines():
            if line.startswith("  -------"):
                in_table = True
                continue
            if in_table and "Number of matches" in line:
                break
            if in_table and line.strip():
                # Line format: "      -2  Mariner 2 (spacecraft)             1962-041A                         "
                parts = line.split(maxsplit=2)
                if len(parts) >= 2:
                    tid = parts[0].strip()
                    # The name is everything after the ID up to the designation. We can just use parts[1] + parts[2] mostly
                    # Using regex or simpler split:
                    full_rest = parts[1] + ((" " + parts[2]) if len(parts) > 2 else "")
                    # The name is usually before the two spaces
                    name_part = full_rest.split("  ")[0].strip()
                    targets.append((name_part, tid))
                    
    except Exception as e:
        print(f"Error fetching all spacecraft from Horizons: {e}")
        
    return targets

def parse_and_find_windows(text, target_name):
    """
    Parses CSV ephemeris text to find windows where Altitude > 0.
    Finds start time, stop time, peak altitude, and determines if it was sunlit.
    """
    if "Error" in text:
        return None
        
    lines = text.splitlines()
    in_data = False
    
    # Store data rows: (time, flags, az, el, mag)
    data_points = []

    for line in lines:
        if "$$SOE" in line:
            in_data = True
            continue
        if "$$EOE" in line:
            in_data = False
            continue
            
        if in_data:
            parts = line.split(',')
            if len(parts) >= 6:
                timestamp = parts[0].strip()
                flags = parts[1].strip()  # Usually daylight/eclipse flags
                try:
                    az = float(parts[2].strip())
                    el = float(parts[3].strip())
                    # Mag might be empty or 'n.a.' if eclipsed
                    mag_str = parts[4].strip()
                    mag = float(mag_str) if mag_str and mag_str != 'n.a.' else None
                    data_points.append((timestamp, flags, az, el, mag))
                except ValueError:
                    continue
    
    windows = []
    current_window = None
    
    for pt in data_points:
        timestamp, flags, az, el, mag = pt
        is_sunlit = mag is not None and 'u' not in flags.lower() # completely eclipsed usually marked by 'N.A.' or 'u'
        
        if el > 0:
            if current_window is None:
                current_window = {
                    "start": timestamp,
                    "end": timestamp,
                    "peak_el": el,
                    "sunlit": is_sunlit,
                    "track": [(az, el)]
                }
            else:
                current_window["end"] = timestamp
                current_window["track"].append((az, el))
                if el > current_window["peak_el"]:
                    current_window["peak_el"] = el
                if is_sunlit:
                    current_window["sunlit"] = True
        else:
            if current_window is not None:
                windows.append(current_window)
                current_window = None
                
    if current_window is not None:
        windows.append(current_window)
        
    return windows

def convert_time_format(time_str):
    """Converts YYYYMMDD_HHmmss to a format Horizons accepts: YYYY-MM-DD HH:MM:SS."""
    try:
        dt = datetime.strptime(time_str, "%Y%m%d_%H%M%S")
        return dt.strftime("%Y-%m-%d %H:%M:%S")
    except ValueError:
        # Fallback to passing the string directly if it doesn't match the expected format
        return time_str

def check_target(target_name, target_id, args, h_start, h_stop):
    print(f"\nChecking {target_name} ({target_id})...")
    result = get_horizons_ephemeris(target_id, args.location, h_start, h_stop, args.step)
    
    if "Multiple primary targets" in result or "No matches found" in result:
        print(f"  -> Could not resolve exact target. Horizons responded with an error or ambiguity.")
        return

    windows = parse_and_find_windows(result, target_name)
    
    if windows is None:
        print(f"  -> Error parsing data.")
        return
        
    if not windows:
        print(f"  -> Not visible above the horizon during this time interval.")
        return []
        
    for w in windows:
        w["sat_id"] = target_id
        
    print(f"  -> VISIBLE! Found {len(windows)} pass(es):")
    for i, w in enumerate(windows, 1):
        sunlit_str = "Yes" if w['sunlit'] else "No (Eclipsed)"
        print(f"     Pass {i}: Start: {w['start']} | End: {w['end']} | Peak Alt: {w['peak_el']:.2f} deg | Sunlit: {sunlit_str}")

    return windows

def check_target_celestrak(target_name, target_id, args, start_str, stop_str):
    print(f"\nChecking {target_name} ({target_id}) via Celestrak...")
    try:
        from skyfield.api import load, wgs84
    except ImportError:
        print("  -> Error: 'skyfield' library is required for Celestrak source.")
        return

    norad_id = None
    if target_id.startswith("DES="):
        norad_id = target_id.replace("DES=", "")
    else:
        # Handle Artemis specifically as Horizons uses negative IDs 
        if "ARTEMIS-P1" in target_name.upper() or target_id == "-117":
            norad_id = "31113"
        elif "ARTEMIS-P2" in target_name.upper() or target_id == "-118":
            norad_id = "31114"
        else:
            norad_id = target_id

    if not str(norad_id).isdigit():
        print(f"  -> Could not determine a valid NORAD ID from '{target_id}'. Celestrak requires NORAD ID.")
        return

    url = f"https://celestrak.org/NORAD/elements/gp.php?CATNR={norad_id}&FORMAT=TLE"
    
    ts = load.timescale()
    try:
        eph = load('de421.bsp')
    except Exception as e:
        print(f"  -> Error loading DE421 ephemeris: {e}")
        return

    try:
        # Load TLE from Celestrak
        satellites = load.tle_file(url, filename=f"tle_{norad_id}.txt")
        if not satellites:
            print(f"  -> No TLE found for NORAD ID {norad_id} on Celestrak.")
            return
        sat = satellites[0]
    except Exception as e:
        print(f"  -> Error fetching TLE from Celestrak: {e}")
        return

    # Parse location
    if args.location.lower() == "elginfield":
        lat, lon, elev = 43.1925, -81.3158, 325.0
    else:
        print("  -> Note: Custom locations not natively supported by CELESTRAK mode yet. Assuming Elginfield.")
        lat, lon, elev = 43.1925, -81.3158, 325.0

    observer = wgs84.latlon(lat, lon, elevation_m=elev)

    try:
        dt_start = datetime.strptime(start_str, "%Y-%m-%d %H:%M:%S")
        dt_stop = datetime.strptime(stop_str, "%Y-%m-%d %H:%M:%S")
    except ValueError:
        print("  -> Error parsing time format.")
        return

    step_minutes = 1
    if args.step.endswith('m'):
        step_minutes = float(args.step[:-1])
    elif args.step.endswith('s'):
        step_minutes = float(args.step[:-1]) / 60.0
        
    if getattr(args, 'plot', False):
        step_minutes /= 4.0
        
    current_dt = dt_start
    windows = []
    current_window = None

    while current_dt <= dt_stop:
        t = ts.utc(current_dt.year, current_dt.month, current_dt.day, current_dt.hour, current_dt.minute, current_dt.second)
        
        difference = sat - observer
        topocentric = difference.at(t)
        alt, az, distance = topocentric.altaz()
        
        el = alt.degrees
        is_sunlit = sat.at(t).is_sunlit(eph)
        timestamp_str = current_dt.strftime("%Y-%b-%d %H:%M:%S")
        
        if el > 0:
            if current_window is None:
                current_window = {
                    "start": timestamp_str,
                    "end": timestamp_str,
                    "peak_el": el,
                    "sunlit": is_sunlit,
                    "track": [(az.degrees, el)],
                    "sat_id": norad_id
                }
            else:
                current_window["end"] = timestamp_str
                current_window["track"].append((az.degrees, el))
                if el > current_window["peak_el"]:
                    current_window["peak_el"] = el
                if is_sunlit:
                    current_window["sunlit"] = True
        else:
            if current_window is not None:
                windows.append(current_window)
                current_window = None
                
        current_dt += timedelta(minutes=step_minutes)

    if current_window is not None:
        windows.append(current_window)

    if not windows:
        print(f"  -> Not visible above the horizon during this time interval.")
        return []
        
    print(f"  -> VISIBLE! Found {len(windows)} pass(es):")
    for i, w in enumerate(windows, 1):
        sunlit_str = "Yes" if w['sunlit'] else "No (Eclipsed)"
        print(f"     Pass {i}: Start: {w['start']} | End: {w['end']} | Peak Alt: {w['peak_el']:.2f} deg | Sunlit: {sunlit_str}")

    return windows

def check_celestrak_group(group_name, args, start_str, stop_str):
    print(f"\nChecking Celestrak group: {group_name.upper()}...")
    try:
        from skyfield.api import load, wgs84
    except ImportError:
        print("  -> Error: 'skyfield' library is required for Celestrak source.")
        return

    CELESTRAK_GROUPS = {
        "active": "active",
        "stations": "stations",
        "brightest": "visual",
        "starlink": "starlink",
        "oneweb": "oneweb",
        "geo": "geo",
        "leo": "leo",
        "last": "last-30-days",
        "gps": "gps-ops"
    }
    
    group_param = CELESTRAK_GROUPS.get(group_name.lower(), group_name.lower())
    url = f"https://celestrak.org/NORAD/elements/gp.php?GROUP={group_param}&FORMAT=TLE"
    
    ts = load.timescale()
    try:
        eph = load('de421.bsp')
    except Exception as e:
        print(f"  -> Error loading DE421 ephemeris: {e}")
        return

    try:
        # Load TLE from Celestrak
        satellites = load.tle_file(url, filename=f"tle_group_{group_param}.txt")
        if not satellites:
            print(f"  -> No TLEs found for group {group_name} on Celestrak.")
            return
        print(f"  -> Loaded {len(satellites)} satellites from {group_name}.")
    except Exception as e:
        print(f"  -> Error fetching TLE group from Celestrak: {e}")
        return

    # Parse location
    if args.location.lower() == "elginfield":
        lat, lon, elev = 43.1925, -81.3158, 325.0
    else:
        print("  -> Note: Custom locations not natively supported by CELESTRAK mode yet. Assuming Elginfield.")
        lat, lon, elev = 43.1925, -81.3158, 325.0

    observer = wgs84.latlon(lat, lon, elevation_m=elev)

    try:
        dt_start = datetime.strptime(start_str, "%Y-%m-%d %H:%M:%S")
        dt_stop = datetime.strptime(stop_str, "%Y-%m-%d %H:%M:%S")
    except ValueError:
        print("  -> Error parsing time format.")
        return

    step_minutes = 1
    if args.step.endswith('m'):
        step_minutes = float(args.step[:-1])
    elif args.step.endswith('s'):
        step_minutes = float(args.step[:-1]) / 60.0
        
    if getattr(args, 'plot', False):
        step_minutes /= 4.0
        
    times_dt = []
    current_dt = dt_start
    while current_dt <= dt_stop:
        times_dt.append(current_dt)
        current_dt += timedelta(minutes=step_minutes)
        
    years = [dt.year for dt in times_dt]
    months = [dt.month for dt in times_dt]
    days = [dt.day for dt in times_dt]
    hours = [dt.hour for dt in times_dt]
    minutes = [dt.minute for dt in times_dt]
    seconds = [dt.second for dt in times_dt]
    t_array = ts.utc(years, months, days, hours, minutes, seconds)
    
    visible_count = 0
    print("  -> Calculating visibility windows...")
    group_passes = []
    
    for sat in satellites:
        difference = sat - observer
        topocentric = difference.at(t_array)
        alt, az, distance = topocentric.altaz()
        
        el_array = alt.degrees
        
        if (el_array > 0).any():
            sunlit_array = sat.at(t_array).is_sunlit(eph)
            windows = []
            current_window = None
            az_array = az.degrees
            
            for idx, el in enumerate(el_array):
                if el > 0:
                    timestamp_str = times_dt[idx].strftime("%Y-%b-%d %H:%M:%S")
                    is_sunlit = bool(sunlit_array[idx])
                    az_val = az_array[idx]
                    
                    if current_window is None:
                        current_window = {
                            "start": timestamp_str,
                            "end": timestamp_str,
                            "peak_el": el,
                            "sunlit": is_sunlit,
                            "track": [(az_val, el)],
                            "sat_id": sat.model.satnum
                        }
                    else:
                        current_window["end"] = timestamp_str
                        current_window["track"].append((az_val, el))
                        if el > current_window["peak_el"]:
                            current_window["peak_el"] = el
                        if is_sunlit:
                            current_window["sunlit"] = True
                else:
                    if current_window is not None:
                        windows.append(current_window)
                        current_window = None
                        
            if current_window is not None:
                windows.append(current_window)
                
            if windows:
                visible_count += 1
                group_passes.append((sat.name, windows))
                print(f"\n  [VISIBLE] {sat.name}")
                for i, w in enumerate(windows, 1):
                    sunlit_str = "Yes" if w['sunlit'] else "No (Eclipsed)"
                    print(f"     Pass {i}: Start: {w['start']} | End: {w['end']} | Peak Alt: {w['peak_el']:.2f} deg | Sunlit: {sunlit_str}")

    print(f"\nFinished processing {group_name}. Found {visible_count} visible satellites.")
    return group_passes

def plot_skymap(all_visible_passes):
    """
    Plots all visibility tracks on a polar projection representing a skymap.
    """
    if not all_visible_passes:
        return
        
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("\nError: matplotlib/numpy is required for plotting. Install with 'pip install matplotlib numpy'.")
        return
        
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, polar=True)
    
    # Polar coordinates in matplotlib:
    # theta (angle) = azimuth in radians
    # r (radius) = zenith angle = 90 - elevation
    # Center is zenith (r=0, alt=90), Edge is horizon (r=90, alt=0)
    
    ax.set_theta_zero_location("N") # North at top
    ax.set_theta_direction(-1)      # Clockwise angles (N=0, E=90, S=180, W=270)
    
    # Set radial limits and ticks (90 degrees to 0 degrees mapped to 0 to 90 radius)
    ax.set_ylim(0, 90)
    ax.set_yticks([0, 30, 60, 90])
    ax.set_yticklabels(['90°\n(Zenith)', '60°', '30°', '0°\n(Horizon)'])
    
    ax.set_xticks(np.radians([0, 45, 90, 135, 180, 225, 270, 315]))
    ax.set_xticklabels(['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'])
    
    # Color map
    colors = plt.cm.tab10(np.linspace(0, 1, 10))
    color_idx = 0
    
    for sat_name, windows in all_visible_passes:
        color = colors[color_idx % 10]
        color_idx += 1
        
        for w_idx, w in enumerate(windows):
            track = w.get("track", [])
            if not track:
                continue
                
            az_vals = [pt[0] for pt in track]
            el_vals = [pt[1] for pt in track]
            
            theta = np.radians(az_vals)
            r = 90 - np.array(el_vals)
            
            label = f"{sat_name}" if w_idx == 0 else ""
            
            ax.plot(theta, r, color=color, linewidth=2, label=label)
            
            if len(theta) > 0:
                ax.plot(theta[0], r[0], 'go', markersize=5) # Start (green)
                ax.plot(theta[-1], r[-1], 'ro', markersize=5) # End (red)
                
    ax.set_title("Visible Satellite Passes", va='bottom', fontsize=14)
    # Legend max 20 items to avoid crowding
    handles, labels = ax.get_legend_handles_labels()
    if len(handles) > 0:
        ax.legend(handles[:20], labels[:20], loc='upper right', bbox_to_anchor=(1.3, 1.1))
    
    plt.tight_layout()
    try:
        plt.savefig("skymap.png", dpi=300)
        print("\nSaved skymap plot to 'skymap.png'.")
    except Exception as e:
        print(f"\nCould not save plot: {e}")
        
    try:
        plt.show(block=False)
        plt.pause(5) # show briefly if running interactively
    except Exception:
        pass

def main():
    now = datetime.utcnow()
    default_start = now.strftime("%Y%m%d_%H%M%S")
    default_stop = (now + timedelta(minutes=10)).strftime("%Y%m%d_%H%M%S")

    parser = argparse.ArgumentParser(description="Find visible satellites from JPL Horizons.")
    parser.add_argument("--targets", type=str, default="", help="Comma-separated list of target names or IDs (e.g. ISS,HST,25544).")
    parser.add_argument("--location", default="elginfield", help="Observer location (default: elginfield)")
    parser.add_argument("--start", default=default_start, help=f"Start time in UTC (default: {default_start})")
    parser.add_argument("--stop", default=default_stop, help=f"Stop time in UTC (default: {default_stop})")
    parser.add_argument("--step", default="1m", help="Step size (default: 1m)")
    parser.add_argument("--source", choices=["horizons", "celestrak"], default="horizons", help="Source for ephemeris (default: horizons)")
    parser.add_argument("--plot", action="store_true", help="Generate a polar skymap plot of visible passes")

    args = parser.parse_args()
    
    h_start = convert_time_format(args.start)
    h_stop = convert_time_format(args.stop)

    print(f"--- Satellite Visibility Checker ---")
    print(f"Location: {args.location}")
    print(f"Source:   {args.source.upper()}")
    print(f"Window:   {args.start} to {args.stop}")
    print(f"Step:     {args.step}")
    print(f"------------------------------------")

    # Determine targets
    targets_to_check = []
    
    if args.targets.lower() == "all":
        if args.source == "horizons":
            print("Fetching list of all spacecraft from JPL Horizons...")
            targets_to_check = get_all_spacecraft()
            print(f"Found {len(targets_to_check)} spacecraft to check.")
        else:
            print("Cannot fetch 'all' spacecraft natively using Celestrak. Please specify individual targets or a group (e.g., 'active', 'brightest').")
            sys.exit(1)
    elif args.targets:
        items = [i.strip() for i in args.targets.split(",")]
        for item in items:
            targets_to_check.append((item, resolve_target(item)))
    else:
        # Default behavior: Check common satellites and Artemis
        print("No targets specified. Checking default common satellites...")
        targets_to_check.append(("ISS", COMMON_SATELLITES["ISS"]))
        targets_to_check.append(("Hubble Space Telescope (HST)", COMMON_SATELLITES["HST"]))
        targets_to_check.append(("China Space Station (CSS)", COMMON_SATELLITES["CSS"]))
        targets_to_check.extend(search_artemis_satellites())

    all_visible_passes = []

    for name, tid in targets_to_check:
        # User requested Artemis explicitly included
        # the list is populated dynamically above.
        windows = None
        if args.source == "horizons":
            windows = check_target(name, tid, args, h_start, h_stop)
        elif args.source == "celestrak":
            if name.lower() in ["active", "stations", "brightest", "starlink", "oneweb", "geo", "leo", "last", "gps"]:
                group_passes = check_celestrak_group(name.lower(), args, h_start, h_stop)
                if group_passes:
                    all_visible_passes.extend(group_passes)
                continue
            else:
                windows = check_target_celestrak(name, tid, args, h_start, h_stop)
                
        if windows:
            all_visible_passes.append((name, windows))
            
    if all_visible_passes:
        print("\n" + "="*125)
        print(f"{'Satellite Name':<25} | {'ID':<8} | {'Start Time':<20} | {'End Time':<20} | {'Start Az':<8} | {'Start El':<8} | {'End Az':<8} | {'End El':<8}")
        print("-" * 125)
        for name, windows in all_visible_passes:
            for w in windows:
                sat_id = w.get("sat_id", "N/A")
                start_time = w.get("start", "N/A")
                end_time = w.get("end", "N/A")
                track = w.get("track", [])
                if track:
                    start_az, start_el = track[0]
                    end_az, end_el = track[-1]
                    print(f"{name[:25]:<25} | {str(sat_id):<8} | {start_time:<20} | {end_time:<20} | {start_az:6.2f}° | {start_el:6.2f}° | {end_az:6.2f}° | {end_el:6.2f}°")
                else:
                    print(f"{name[:25]:<25} | {str(sat_id):<8} | {start_time:<20} | {end_time:<20} | {'N/A    '} | {'N/A    '} | {'N/A    '} | {'N/A    '}")
        print("="*125 + "\n")
            
    if getattr(args, 'plot', False) and all_visible_passes:
        plot_skymap(all_visible_passes)

if __name__ == "__main__":
    main()
