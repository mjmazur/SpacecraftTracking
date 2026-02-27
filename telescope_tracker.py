import time
import argparse
import sys
from datetime import datetime, timedelta
from alpyca.alpaca import Telescope
from alpyca.alpaca import DeviceAddress
import requests

def get_current_ephemeris(target, location):
    """
    Fetches the current ephemeris for a specific target and location from JPL Horizons.
    We request a very small window around 'now' to get near real-time data.
    """
    url = "https://ssd.jpl.nasa.gov/api/horizons.api"
    
    # Standardize target
    if target.upper() == "ISS":
        target = "DES=25544"
    
    now = datetime.utcnow()
    # Request a short window to ensure we get a valid point for 'now'
    start_time = now.strftime("%Y-%m-%d %H:%M:%S")
    stop_time = (now + timedelta(minutes=2)).strftime("%Y-%m-%d %H:%M:%S")

    params = {
        "format": "text",
        "COMMAND": f"'{target}'",
        "OBJ_DATA": "NO",
        "MAKE_EPHEM": "YES",
        "EPHEM_TYPE": "OBSERVER",
        "CENTER": f"'{location}'",
        "START_TIME": f"'{start_time}'",
        "STOP_TIME": f"'{stop_time}'",
        "STEP_SIZE": "'1m'",
        "QUANTITIES": "'1,2,3,4,9'", # RA/Dec (App), Rates, Mag, El/Az
    }

    try:
        response = requests.get(url, params=params)
        response.raise_for_status()
        return response.text
    except Exception as e:
        print(f"Error fetching data from Horizons: {e}")
        return None

def parse_current_data(text):
    """
    Parses the ephemeris text for the first data row after $$SOE.
    Returns a dict with RA, Dec, dRA*cosD, dDec, and Alt.
    """
    if not text:
        return None
        
    lines = text.splitlines()
    data_line = None
    for i, line in enumerate(lines):
        if "$$SOE" in line:
            if i + 1 < len(lines):
                data_line = lines[i+1]
            break
            
    if not data_line:
        return None
        
    # Example format: 2026-Feb-27 00:39:00.000 *m  01 23 45.67 +12 34 56.7  0.123  0.456  ...
    # Note: Columns depend on QUANTITIES. 
    # Q1,2,3,4,9 -> RA (Astrom), RA (App), Dec (App), dRA*cosD, dDec, Az, Alt
    parts = data_line.split()
    
    # Finding columns is tricky with variable spacing. 
    # Usually: Date(0-1) Time(2) RA_App(3-5) Dec_App(6-8) dRA(9) dDec(10) Az(11) Alt(12)
    # But wait, date and time might be merged or split.
    
    try:
        # Simple extraction based on typical Horizons spacing
        # This is a bit fragile, a regex or fixed-width parser would be better but let's try this first.
        # We look for the parts after the timestamp.
        
        # Horizons rows start with: YYYY-Mon-DD HH:MM:SS.SSS
        # Then sometimes flags like '*' or 'm'
        # Then the data.
        
        # Let's find where the RA/Dec start. They usually follow the time.
        # Parts 0, 1, 2 are usually date/time.
        data_start_idx = 3
        while data_start_idx < len(parts) and any(c.isalpha() for c in parts[data_start_idx]):
            data_start_idx += 1
            
        ra_h = float(parts[data_start_idx])
        ra_m = float(parts[data_start_idx+1])
        dra_arcsec_hr = float(parts[data_start_idx+6])
        ddec_arcsec_hr = float(parts[data_start_idx+7])
        
        # Skip Az (idx+8)
        alt = float(parts[data_start_idx+9])
        
        return {
            "ra": ra_h + ra_m/60.0 + ra_s/3600.0,
            "dec": dec_d + (dec_m/60.0 + dec_s/3600.0) * (1 if dec_d >= 0 else -1),
            "dra_arcsec_sec": dra_arcsec_hr / 3600.0,
            "ddec_arcsec_sec": ddec_arcsec_hr / 3600.0,
            "alt": alt
        }
    except Exception as e:
        print(f"Error parsing data line: {e}")
        print(f"Line content: {data_line}")
        return None

def main():
    parser = argparse.ArgumentParser(description="Track objects using ASCOM telescope and JPL Horizons.")
    parser.add_argument("--target", default="499", help="Target body (default: 499)")
    parser.add_argument("--location", default="440", help="Observer location (default: 440)")
    parser.add_argument("--address", default="localhost:5959", help="Alpaca server address (default: localhost:5959)")
    parser.add_argument("--device_number", type=int, default=0, help="Telescope device number (default: 0)")
    parser.add_argument("--max_dra", type=float, default=2.5, help="Max RA track rate in deg/sec (default: 2.5)")
    parser.add_argument("--max_ddec", type=float, default=2.5, help="Max Dec track rate in deg/sec (default: 2.5)")
    parser.add_argument("--interval", type=int, default=10, help="Update interval in seconds (default: 10)")

    args = parser.parse_args()

    # Connect to telescope
    addr = DeviceAddress(args.address, args.device_number)
    telescope = Telescope(addr)

    try:
        print(f"Connecting to telescope at {args.address}...")
        telescope.Connected = True
        print(f"Connected to: {telescope.Name}")
        
        if not telescope.CanSetTracking:
            print("Error: Telescope does not support setting tracking rates.")
            sys.exit(1)

        while True:
            print(f"\n[{datetime.now().strftime('%H:%M:%S')}] Fetching current position for {args.target}...")
            raw_data = get_current_ephemeris(args.target, args.location)
            data = parse_current_data(raw_data)
            
            if not data:
                print("Failed to get or parse target data.")
                time.sleep(args.interval)
                continue
                
            print(f"Target Info: Alt={data['alt']:.2f}, RA={data['ra']:.4f}, Dec={data['dec']:.4f}")
            
            # Checks
            if data['alt'] < 0:
                print("Error: Object is below the horizon!")
                sys.exit(1)
            
            # Convert rates from arcsec/sec to deg/sec for limit check
            dra_deg_sec = abs(data['dra_arcsec_sec']) / 3600.0
            ddec_deg_sec = abs(data['ddec_arcsec_sec']) / 3600.0
            
            if dra_deg_sec > args.max_dra or ddec_deg_sec > args.max_ddec:
                print(f"Error: Required tracking rate exceeds limits!")
                print(f"Required: dRA={dra_deg_sec:.4f}, dDec={ddec_deg_sec:.4f} deg/s")
                print(f"Limits: RA={args.max_dra}, Dec={args.max_ddec} deg/s")
                sys.exit(1)
            
            # Slew and Track
            if not telescope.AtPark:
                print(f"Slewing to Target...")
                telescope.SlewToCoordinatesAsync(data['ra'], data['dec'])
                
                # Wait for slew to complete (simplified)
                while telescope.Slewing:
                    time.sleep(1)
                
                # Set tracking rates
                # ASCOM tracking rates are in arcseconds per SI second.
                print(f"Setting tracking rates: dRA={data['dra_arcsec_sec']:.4f}, dDec={data['ddec_arcsec_sec']:.4f} arcsec/s")
                telescope.RightAscensionRate = data['dra_arcsec_sec']
                telescope.DeclinationRate = data['ddec_arcsec_sec']
                telescope.Tracking = True
            else:
                print("Telescope is parked. Unpark to start tracking.")

            time.sleep(args.interval)

    except KeyboardInterrupt:
        print("\nTracking stopped by user.")
    except Exception as e:
        print(f"An error occurred: {e}")
    finally:
        if telescope.Connected:
            telescope.Tracking = False
            telescope.Connected = False
            print("Disconnected from telescope.")

if __name__ == "__main__":
    main()
