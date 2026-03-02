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
- **Safety Measures**: Integrated checks to prevent slewing to objects below the horizon or exceeding hardware tracking limits.
- **Daylight Filtering**: Option to filter out daylight observations for optical tracking.
- **Environment Management**: Optimized for use with Conda environments (e.g., `sattrack`).

## Installation

### Prerequisites

- Python 3.9 or later.
- An [ASCOM Alpaca](https://ascom-standards.org/Alpaca/Index.htm) compliant telescope driver or simulator (e.g., [ASCOM Remote](https://ascom-standards.org/Downloads/Index.htm)).

### Setup

1. **Clone the repository**:
   ```powershell
   git clone https://github.com/mjmaz/SpacecraftTracking.git
   cd SpacecraftTracking
   ```

2. **Create and activate a virtual environment** (recommended):
   ```powershell
   conda create -n sattrack python=3.10
   conda activate sattrack
   ```

3. **Install dependencies**:
   ```powershell
   pip install -r requirements.txt
   ```
   *Note: Ensure `alpyca` (the Python Alpaca interface) is installed.*

## Usage

### 1. Listing and Searching Objects
Use `list_objects.py` to find the ID or name of the target you wish to track.

- **List major bodies**:
  ```powershell
  python list_objects.py --major
  ```
- **Search for small bodies/spacecraft**:
  ```powershell
  python list_objects.py --search "ISS"
  ```

### 2. Fetching Ephemeris Data
Use `horizons_client.py` to get a detailed observer table for a specific location.

- **Example for Mars from a specific location**:
  ```powershell
  python horizons_client.py --target 499 --location 440 --start "2026-03-01 00:00" --stop "2026-03-01 12:00"
  ```

### 3. Automated Telescope Tracking
Use `telescope_tracker.py` to connect to your telescope and begin tracking a target in real-time.

- **Track the ISS**:
  ```powershell
  python telescope_tracker.py --target ISS --address localhost:5959 --device_number 0
  ```
- **Arguments**:
  - `--target`: JPL Horizons target ID (e.g., `ISS`, `499` for Mars).
  - `--location`: JPL Horizons observer code (e.g., `440` for Elginfield).
  - `--address`: IP:Port of the Alpaca server.
  - `--max_dra / --max_ddec`: Maximum allowed tracking rates (deg/sec).
  - `--interval`: Update frequency in seconds.

## Project Structure

- `telescope_tracker.py`: Main script for telescope control and real-time tracking.
- `horizons_client.py`: Utility for fetching and parsing JPL Horizons ephemeris.
- `list_objects.py`: Utility for discovering target IDs and names.
- `requirements.txt`: Project dependencies.

## License

This project is intended for educational and enthusiast astronomical use. Please ensure compliance with local flight safety and laser/telescope regulations.
