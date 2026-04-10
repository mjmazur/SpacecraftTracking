# Spacecraft Tracking and Telescope Control

A suite of Python scripts for retrieving spacecraft ephemeris data from JPL Horizons and controlling ASCOM-compliant telescopes for automated tracking.

## Overview

This project provides tools to:
- Search and list celestial objects (planets, asteroids, comets, and spacecraft like the ISS) using the JPL Horizons API.
- Fetch detailed, observer-specific ephemeris data (RA/Dec, rates, magnitudes, etc.).
- Automatically slew and track objects using an ASCOM Alpaca telescope interface, including safety checks for horizon limits and maximum tracking rates.

## Features

- **JPL Horizons Integration**: Real-time access to high-precision orbital data.
- **Automated Tracking**: Slew to and track fast-moving objects (like the ISS) with an ASCOM telescope.
- **Interactive Controls**: Real-time keyboard shortcuts to halt movement or park the telescope safely.
- **Safety Measures**: Integrated checks to prevent slewing to objects below the horizon or exceeding hardware tracking limits.
- **Environment Management**: Optimized for use with Conda environments (e.g., `sattrack`).

## Installation

### Prerequisites

- Python 3.9 or later.
- Windows (required for the interactive keyboard shortcuts).
- An [ASCOM Alpaca](https://ascom-standards.org/Alpaca/Index.htm) compliant telescope driver or simulator (e.g., [ASCOM Remote](https://ascom-standards.org/Downloads/Index.htm)).

### Setup

1. **Clone the repository**:
   ```powershell
   git clone https://github.com/mjmaz/SpacecraftTracking.git
   cd SpacecraftTracking
   ```

2. **Create and activate the virtual environment**:
   ```powershell
   conda create -n sattrack python=3.10
   conda activate sattrack
   ```

3. **Install dependencies**:
   ```powershell
   pip install -r requirements.txt
   ```

## Usage: telescope_tracker.py

The main script for real-time tracking. It fetches ephemeris data every few seconds and updates the telescope's position and custom tracking rates.

### Basic Command
```powershell
python telescope_tracker.py --target ISS --address localhost:5959
```

### Interactive Keyboard Shortcuts (Windows)
While the tracker is running, the following shortcuts can be used to control the telescope instantly:
- **CTRL + S** or **CTRL + H**: **Halt Movement**. Immediately aborts any active slew and disables tracking.
- **CTRL + P**: **Park Telescope**. Slews the telescope to its parking location (**Azimuth 90°, Altitude 0°**) and then issues the Park command.

### Arguments
- `--target`: JPL Horizons target ID. Use `ISS` for the International Space Station, or numerically (e.g., `599` for Venus, `499` for Mars).
- `--location`: JPL Horizons observer code (default: `440` for Elginfield).
- `--address`: IP and Port of the Alpaca server (default: `localhost:5959`).
- `--device_number`: The integer index of the telescope in your Alpaca server (usually `0`).
- `--max_dra` / `--max_ddec`: Maximum allowed tracking rates in **degrees per second** (default: `2.5`). The script will refuse to track if the target exceeds these limits.
- `--interval`: How often (in seconds) to fetch a new ephemeris point and update the telescope position (default: `10`).
- `--slew_threshold`: The angular distance (in **arc-seconds**) at which the telescope should stop its slew and begin tracking (default: `30.0`).

## Usage: Utility Scripts

### 1. Listing and Searching Objects (`list_objects.py`)
- **List major bodies** (Planets, Moon, etc.):
  ```powershell
  python list_objects.py --major
  ```
- **Search for spacecraft or small bodies**:
  ```powershell
  python list_objects.py --search "ISS"
  ```

### 2. Fetching Ephemeris Data (`horizons_client.py`)
Fetches and displays a full ephemeris table for a given time range.
- **Example**:
  ```powershell
  python horizons_client.py --target 499 --start "2026-03-01 00:00" --stop "2026-03-01 12:00"
  ```

## Advanced Algorithms & Physics Modules

### 1. `horizons_tle.py` (Vector-Based Ephemeris Extraction)
This tool calculates TLE parameters algebraically directly from state variables instead of relying on default JPL element tables (which frequently omit drag/motion parameters for missions).
- **Cartesian Mechanics**: It queries the Horizons `VECTORS` payload (ICRF frame) to harvest $X$, $Y$, $Z$ and $V_X$, $V_Y$, $V_Z$.
- **SGP4 Inversion**: These variables are run backwards through `sgp4.ext.rv2coe` to formulate instantaneous Classical Orbital Elements (COE), effectively bypassing conventional TLE gaps while maintaining precise `astropy.time` epoch synchronization.

### 2. `query_horizons.py` (Advanced Ephemeris Profiler)
This script models absolute physical conditions and geometric transits for telescope predictions:
- **Apparent Magnitude Extinction**: Artemis II (and analog tracking objects) brightness is modeled utilizing a mathematical Diffuse Sphere projection standardizing a $20 m^2$ frame. It additionally integrates the **1989 Kasten Young Geometric Airmass Model** ($k=0.2$ visual baseline) to reliably attenuate visible brightness down toward the atmospheric horizon.
- **Lunar Limb Angular Separation**: Transits are predicted via a dual-API payload. The script grabs Topocentric Apparent RA/DEC coordinates for both the target spacecraft and the Moon concurrently. It merges them inside an `astropy.coordinates.SkyCoord` engine to calculate the absolute degree of spherical separation natively subtracting the dynamic lunar radius sequence!
- **Mechanical Rate of Motion**: The tracker intercepts raw Horizons Output Apparent Rate arrays (`dRA*cosD` and `dDEC/dt`), parses their strings out of standard `arcsec/hour` variables, and mathematically scales the true tracking rate limits sequentially into `arcsec/s` suitable for precision Alpaca Slew constraints.

## Technical Details

### Tracking Logic
The script uses the ASCOM `RightAscensionRate` and `DeclinationRate` properties to achieve smooth tracking of fast-moving objects. It assumes the mount is in **Equatorial mode**.

### Safety
- **Horizon Check**: If the target's altitude drops below 0°, the script will stop tracking and exit.
- **Connection Loss**: On any unexpected error or user interruption (CTRL+C), the script attempts to disable tracking and disconnect cleanly.

## License

This project is intended for educational and astronomical use. Please ensure compliance with local flight safety and laser/telescope regulations.
