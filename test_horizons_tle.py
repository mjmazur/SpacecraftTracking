import requests

url = "https://ssd.jpl.nasa.gov/api/horizons.api"

params = {
    "format": "text",
    "COMMAND": "'-1024'",
    "OBJ_DATA": "YES",
    "MAKE_EPHEM": "YES",
    "EPHEM_TYPE": "VECTORS",
    "CENTER": "'500@399'",
    "REF_PLANE": "'FRAME'",
    "VEC_TABLE": "'2'",
    "START_TIME": "'2026-04-03'",
    "STOP_TIME": "'2026-04-04'",
    "STEP_SIZE": "'1d'",
    "OUT_UNITS": "'KM-S'",
}

response = requests.get(url, params=params)
print(response.text)
