import requests
import argparse
from datetime import datetime
import re

def get_horizons_elements(target, epoch_date):
    """
    Query JPL Horizons for orbital elements of a target object.
    """
    url = "https://ssd.jpl.nasa.gov/api/horizons.api"
    # Calculate a stop time 1 day after the epoch date to satisfy Horizons API
    try:
        import datetime as dt
        try:
            start_dt = datetime.strptime(epoch_date, "%Y-%m-%d %H:%M:%S")
            stop_date = (start_dt + dt.timedelta(days=1)).strftime("%Y-%m-%d %H:%M:%S")
        except ValueError:
            start_dt = datetime.strptime(epoch_date, "%Y-%m-%d")
            stop_date = (start_dt + dt.timedelta(days=1)).strftime("%Y-%m-%d")
    except ValueError:
        stop_date = epoch_date # fallback if it's some other format
        
    params = {
        "format": "text",
        "COMMAND": f"'{target}'",
        "OBJ_DATA": "YES",
        "MAKE_EPHEM": "YES",
        "EPHEM_TYPE": "VECTORS",
        "CENTER": "'500@399'",  # Earth barycenter for TLEs
        "REF_PLANE": "'FRAME'",
        "VEC_TABLE": "'2'",
        "OUT_UNITS": "'KM-S'",
        "START_TIME": f"'{epoch_date}'",
        "STOP_TIME": f"'{stop_date}'",
        "STEP_SIZE": "'1d'",
    }
    
    try:
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()
        text = response.text
        if "$$SOE" not in text:
            raise Exception("Horizons response did not contain ephemeris data ($$SOE not found):\n" + text[:500])
        return text
    except requests.exceptions.RequestException as e:
        raise Exception(f"Error querying JPL Horizons: {e}")

import math
from sgp4.ext import rv2coe
from sgp4.api import WGS84

def parse_elements(text):
    """
    Parse the Horizons VECTORS output to extract orbital parameters using SGP4 rv2coe.
    """
    elements = {}
    
    # Look for Target body name
    name_match = re.search(r'Target body name:\s*(.*?)\s+(?:\{|\n)', text)
    if name_match:
        elements['name'] = name_match.group(1).strip()
    else:
        elements['name'] = "UNKNOWN TARGET"
        
    # Use regex to find X, Y, Z, VX, VY, VZ
    matches = re.finditer(r'([A-Za-z]+)\s*=\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)', text)
    data = {m.group(1).upper(): float(m.group(2)) for m in matches}
    
    # Look for Epoch time
    if 'JDTDB' in text:
        # Example format: 2461133.500000000 = A.D. 2026-Apr-03 00:00:00.0000 TDB
        jd_match = re.search(r'([0-9]+\.[0-9]+)\s*=\s*A\.D\.', text)
        if jd_match:
            jd = float(jd_match.group(1))
            from astropy.time import Time
            t = Time(jd, format='jd', scale='tdb')
            elements['epoch'] = t.utc.datetime
        else:
            elements['epoch'] = datetime.utcnow()
    else:
        elements['epoch'] = datetime.utcnow()

    r = [data.get('X', 0.0), data.get('Y', 0.0), data.get('Z', 0.0)]
    v = [data.get('VX', 0.0), data.get('VY', 0.0), data.get('VZ', 0.0)]
    
    # Earth's standard gravitational parameter (WGS-84) in km^3/s^2
    mu = 398600.4418
    
    # Calculate Classical Orbital Elements via sgp4
    p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper = rv2coe(r, v, mu)
    
    elements['ecc'] = ecc
    # Convert radians to degrees for TLE format
    elements['inc'] = math.degrees(incl)
    elements['raan'] = math.degrees(omega) % 360.0
    elements['aop'] = math.degrees(argp) % 360.0
    elements['ma'] = math.degrees(m) % 360.0
    
    # Calculate Mean Motion in revs per day
    if a > 0:
        mm_rad_sec = math.sqrt(mu / (a**3))
        elements['mm'] = mm_rad_sec * (86400.0 / (2 * math.pi))
    else:
        # Hyperbolic or parabolic orbits can't have a standard mean motion (period is infinite).
        elements['mm'] = 0.0
    
    return elements

def format_tle(satnum, elements, epoch):
    """
    Format the elements into a TLE (Two-Line Element set).
    """
    # Title Line (Line 0)
    # Most TLE parsers expect a title line up to 24 chars long, but we'll include the full name extracted
    title_line = elements.get('name', f"TARGET {satnum}")
    
    # Line 1
    epoch_frac = epoch.hour / 24.0 + epoch.minute / 1440.0 + epoch.second / 86400.0 + epoch.microsecond / 86400000000.0
    epoch_frac_str = f"{epoch_frac:.8f}"[1:] if epoch_frac > 0 else ".00000000"
    epoch_day = epoch.timetuple().tm_yday
    
    line1 = f"1 {satnum:05d}U 00000A   {epoch.year % 100:02d}{epoch_day:03d}{epoch_frac_str}  .00000000  00000-0  00000-0 0  0001"
    
    # Line 2
    inc = elements['inc']
    raan = elements['raan']
    ecc = elements['ecc']
    aop = elements['aop']
    ma = elements['ma']
    mm = elements['mm']
    
    # Format ecc as 7 digits after decimal, but in TLE it's 0000000 for 0.0000000
    ecc_str = f"{int(ecc * 10000000):07d}"
    line2 = f"2 {satnum:05d} {inc:8.4f} {raan:8.4f} {ecc_str} {aop:8.4f} {ma:8.4f} {mm:11.8f}    01"
    
    # Calculate checksums
    line1 += str(calculate_checksum(line1))
    line2 += str(calculate_checksum(line2))
    
    return title_line + '\n' + line1 + '\n' + line2

def calculate_checksum(line):
    """
    Calculate the TLE checksum for a line.
    """
    checksum = 0
    for char in line:
        if char.isdigit():
            checksum += int(char)
        elif char == '-':
            checksum += 1
    return checksum % 10

def main():
    parser = argparse.ArgumentParser(description="Query JPL Horizons for an object and generate a TLE file.")
    parser.add_argument("--target", default="-1024", help="Target object ID (e.g., '-1024' or 'DES=25544')")
    parser.add_argument("--epoch", default=datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"), help="Epoch date for elements (YYYY-MM-DD HH:MM:SS)")
    parser.add_argument("--output", default=None, help="Output TLE filename")
    
    args = parser.parse_args()
    
    # Extract satellite number from target
    if args.target.startswith('DES='):
        satnum = int(args.target.split('=')[1])
    else:
        try:
            satnum = abs(int(args.target))
            if satnum > 99999:
                satnum = 99999
        except ValueError:
            # For other objects, use a dummy number
            satnum = 99999
    
    try:
        text = get_horizons_elements(args.target, args.epoch)
        elements = parse_elements(text)
        
        if not elements:
            print("Error: Could not parse orbital elements from Horizons response.")
            print("Response:")
            print(text)
            return
        
        epoch = elements['epoch']
        tle = format_tle(satnum, elements, epoch)
        
        if args.output:
            filename = args.output
        else:
            filename = f"{args.target.replace('=', '_')}.tle"
        
        with open(filename, 'w') as f:
            f.write(tle)
        
        print(f"TLE saved to {filename}")
        print(tle)
        
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()