import urllib.request
import urllib.parse
import sys

url = "https://ssd.jpl.nasa.gov/api/horizons.api"

def test_cmd(cmd):
    params = {
        "format": "text",
        "COMMAND": cmd,
    }
    
    try:
        full_url = url + "?" + urllib.parse.urlencode(params)
        print("URL:", full_url)
        with urllib.request.urlopen(full_url) as response:
            text = response.read().decode('utf-8')
            print(f"OUTPUT for {cmd}:")
            lines = text.splitlines()
            for line in lines[:30]:
                print(line)
            print("...")
            for line in lines[-20:]:
                print(line)
    except Exception as e:
        print(f"ERROR: {cmd}: {e}")

test_cmd("SPACECRAFT")
