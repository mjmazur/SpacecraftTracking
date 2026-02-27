import requests
import argparse
import sys

def list_major_bodies():
    """
    Fetches and prints the list of major bodies from JPL Horizons.
    """
    url = "https://ssd.jpl.nasa.gov/api/horizons.api"
    params = {
        "format": "text",
        "COMMAND": "'MB'"
    }
    
    try:
        response = requests.get(url, params=params)
        response.raise_for_status()
        print("\n--- Major Bodies (ID: Name) ---")
        print(response.text)
    except Exception as e:
        print(f"Error listing major bodies: {e}")

def search_small_bodies(query):
    """
    Searches for small bodies (asteroids/comets) matching a query.
    """
    # The Horizons Lookup API is better for searching
    url = "https://ssd.jpl.nasa.gov/api/horizons_lookup.api"
    params = {
        "format": "json",
        "sstr": query
    }
    
    try:
        response = requests.get(url, params=params)
        response.raise_for_status()
        data = response.json()
        
        if not data.get("result"):
            print(f"\nNo objects found matching '{query}'.")
            return

        print(f"\n--- Search Results for '{query}' ---")
        for item in data["result"]:
            name = item.get("name", "Unknown")
            spkid = item.get("spkid", "N/A")
            des = item.get("des", "")
            print(f"ID: {spkid:<10} Name: {name} ({des})")
            
    except Exception as e:
        print(f"Error searching small bodies: {e}")

def main():
    parser = argparse.ArgumentParser(description="List or search for objects in JPL Horizons.")
    parser.add_argument("--search", help="Search for a specific object (small bodies)")
    parser.add_argument("--major", action="store_true", default=True, help="List major bodies (default: True)")

    args = parser.parse_args()

    if args.search:
        search_small_bodies(args.search)
    else:
        list_major_bodies()

if __name__ == "__main__":
    main()
