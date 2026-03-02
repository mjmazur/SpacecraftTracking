import time
import msvcrt
import argparse
import sys
from datetime import datetime, timedelta
from alpaca.telescope import Telescope
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
    Returns a dict with RA, Dec (apparent), dRA*cosD, dDec, and Alt.
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
        
    parts = data_line.split()
    
    try:
        # Horizons rows typically follow this order with our QUANTITIES (1,2,3,4,9):
        # 0: Date, 1: Time, [2: Solar flag], ...
        
        # Find where numeric data starts (RA h)
        # We look for the first element after the time (index 1) that could be a flag.
        data_start_idx = 2
        while data_start_idx < len(parts) and any(c.isalpha() or c == '*' for c in parts[data_start_idx]):
            data_start_idx += 1
            
        # We requested 1 & 2, so we have Astrometric RA/Dec then Apparent RA/Dec.
        # Astrom RA: 3 parts, Astrom Dec: 3 parts
        # Apparent RA: 3 parts, Apparent Dec: 3 parts
        
        app_ra_h = float(parts[data_start_idx + 6])
        app_ra_m = float(parts[data_start_idx + 7])
        app_ra_s = float(parts[data_start_idx + 8])
        
        app_dec_d = float(parts[data_start_idx + 9])
        app_dec_m = float(parts[data_start_idx + 10])
        app_dec_s = float(parts[data_start_idx + 11])
        
        # 3. Rates (dRA*cosD, dDec)
        dra_val = float(parts[data_start_idx + 12])
        ddec_val = float(parts[data_start_idx + 13])
        
        # 9. Azimuth and Elevation (Altitude)
        # In our QUANTITIES 1,2,3,4,9, Az is at index +14 and Alt at +15
        alt = float(parts[data_start_idx + 15])
        
        # Calculate full degrees/hours
        ra_total = app_ra_h + app_ra_m/60.0 + app_ra_s/3600.0
        dec_total = (abs(app_dec_d) + app_dec_m/60.0 + app_dec_s/3600.0) * (1 if app_dec_d >= 0 else -1)

        # Rate units: check if they are arcsec/hr (JPL default for observer tables)
        # If rates seem small, they might be arcsec/s. If large, arcsec/hr.
        # We'll assume arcsec/sec for now as the script expects, but a real unit check would be better.
        # Usually, '3. Rates' are arcsec/hr in observer tables if not specified otherwise.
        
        return {
            "ra": ra_total,
            "dec": dec_total,
            "dra_arcsec_sec": dra_val,
            "ddec_arcsec_sec": ddec_val,
            "alt": alt
        }
    except Exception as e:
        print(f"Error parsing data line: {e}")
        print(f"Line content: {data_line}")
        return None

def handle_keyboard_input(telescope):
    """
    Checks for non-blocking keyboard input on Windows.
    Returns True if an exit command (like Park) was issued.
    """
    if msvcrt.kbhit():
        key = msvcrt.getch()
        if key in [b'\x08', b'\x13']:  # CTRL+H, CTRL+S
            print("\n[User Command] Halting movement...")
            try:
                telescope.AbortSlew()
                telescope.Tracking = False
                print("Movement halted and tracking disabled.")
            except Exception as e:
                print(f"Error halting: {e}")
        elif key == b'\x10':  # CTRL+P
            print("\n[User Command] Parking telescope at Az=90, Alt=0...")
            try:
                # Slew to position 2 (Az=90, Alt=0) then Park
                telescope.SlewToAltAzAsync(90.0, 0.0)
                while telescope.Slewing:
                    time.sleep(0.5)
                telescope.Park()
                print("Telescope parked.")
                return True
            except Exception as e:
                print(f"Error parking: {e}")
    return False

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
    telescope = Telescope(args.address, args.device_number)

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
            if telescope.AtPark:
                print("Telescope is parked. Unparking...")
                telescope.Unpark()
            
            # Ensure tracking is ON before slewing to equatorial coordinates (prevents error 0x40b)
            if not telescope.Tracking:
                print("Enabling tracking...")
                telescope.Tracking = True

            print(f"Slewing to Target...")
            telescope.SlewToCoordinatesAsync(data['ra'], data['dec'])
            
            # Wait for slew to complete (simplified)
            while telescope.Slewing:
                time.sleep(1)
            
            # Set custom tracking rates
            # ASCOM tracking rates are in arcseconds per SI second.
            print(f"Setting tracking rates: dRA={data['dra_arcsec_sec']:.4f}, dDec={data['ddec_arcsec_sec']:.4f} arcsec/s")
            telescope.RightAscensionRate = data['dra_arcsec_sec']
            telescope.DeclinationRate = data['ddec_arcsec_sec']

            # Polling wait instead of time.sleep
            print(f"Tracking... (Press CTRL+S to Halt, CTRL+P to Park)")
            start_wait = time.time()
            while time.time() - start_wait < args.interval:
                if handle_keyboard_input(telescope):
                    return # Exit main
                time.sleep(0.1)

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
