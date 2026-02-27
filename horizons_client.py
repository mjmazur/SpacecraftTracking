import requests
import argparse
from datetime import datetime, timedelta
import sys

def get_horizons_ephemeris(target, location, start_time, stop_time, step_size):
    """
    Connects to JPL Horizons API and retrieves ephemeris data.
    """
    url = "https://ssd.jpl.nasa.gov/api/horizons.api"
    
    # Horizons API parameters
    # COMMAND: Target body (e.g., '25544' for ISS)
    # CENTER: Observer location (e.g., coord@399 for lat/lon on Earth)
    # MAKE_EPHEM: 'YES' for ephemeris
    # EPHEM_TYPE: 'OBSERVER'
    # START_TIME, STOP_TIME, STEP_SIZE: Time range and resolution
    # QUANTITIES: Which data points to return (defaulting to essential ones)
    
    # For ISS, we often need to prefix with 'DES=' or use the ID
    if target.upper() == "ISS":
        target = "DES=25544"
    
    # Elginfield Coordinates: 43.1925 N, 81.3158 W (represented as negative for West)
    # Format: 'coord@399' with site coords defined as SITE_COORD='lon,lat,alt'
    # However, Horizons also accepts some site names or we can use custom coordinates.
    # We'll use the specific coordinates for elginfield if requested.
    
    if location.lower() == "elginfield":
        # lon, lat, alt (km)
        location_param = "coord@399"
        site_coord = "-81.3158,43.1925,0.325"
    else:
        location_param = location

    params = {
        "format": "text",
        "COMMAND": f"'{target}'",
        "OBJ_DATA": "YES",
        "MAKE_EPHEM": "YES",
        "EPHEM_TYPE": "OBSERVER",
        "CENTER": f"'{location_param}'",
        "START_TIME": f"'{start_time}'",
        "STOP_TIME": f"'{stop_time}'",
        "STEP_SIZE": f"'{step_size}'",
        "QUANTITIES": "'1,2,3,9,20,23,24'",  # RA/Dec (Astr & App), Rates, Mag, Range, etc.
    }
    
    if location.lower() == "elginfield":
        params["SITE_COORD"] = f"'{site_coord}'"

    try:
        response = requests.get(url, params=params)
        response.raise_for_status()
        return response.text
    except requests.exceptions.RequestException as e:
        return f"Error connecting to JPL Horizons: {e}"

def parse_and_filter(text, filter_daylight=True):
    """
    Parses the Horizons output and filters out rows based on daylight status.
    """
    lines = text.splitlines()
    output = []
    in_data = False
    
    header = []
    footer = []
    data_rows = []

    for line in lines:
        if "$$SOE" in line:
            in_data = True
            output.append(line)
            continue
        if "$$EOE" in line:
            in_data = False
            output.append(line)
            continue
        
        if in_data:
            # Check for daylight indicator '*'
            # In standard observer tables, '*' indicates daylight.
            # It usually appears right after the timestamp.
            if filter_daylight and "*" in line[:25]: # Look in the first part of the line
                continue
            output.append(line)
        else:
            output.append(line)
            
    return "\n".join(output)

def main():
    # Set default times
    now = datetime.utcnow()
    default_start = now.strftime("%Y-%m-%d %H:%M")
    default_stop = (now + timedelta(days=1)).strftime("%Y-%m-%d %H:%M")

    parser = argparse.ArgumentParser(description="Fetch ephemeris data from JPL Horizons.")
    parser.add_argument("--target", default="499", help="Target body (default: 499)")
    parser.add_argument("--location", default="440", help="Observer location (default: 440)")
    parser.add_argument("--start", default=default_start, help=f"Start time (default: {default_start})")
    parser.add_argument("--stop", default=default_stop, help=f"Stop time (default: {default_stop})")
    parser.add_argument("--step", default="15m", help="Step size (default: 15m)")
    parser.add_argument("--no-filter", action="store_true", help="Do not filter out daylight entries")

    args = parser.parse_args()

    print(f"--- Fetching Ephemeris ---")
    print(f"Target:   {args.target}")
    print(f"Location: {args.location}")
    print(f"Start:    {args.start}")
    print(f"Stop:     {args.stop}")
    print(f"Step:     {args.step}")
    print(f"Filtering Daylight: {'No' if args.no_filter else 'Yes'}")
    print(f"--------------------------\n")

    result = get_horizons_ephemeris(args.target, args.location, args.start, args.stop, args.step)
    
    if "Error" not in result:
        filtered_result = parse_and_filter(result, filter_daylight=not args.no_filter)
        print(filtered_result)
    else:
        print(result)

if __name__ == "__main__":
    main()
