#!/usr/bin/env python3
"""
Author : Emmanuel Gonzalez, Jeffrey Demieville
Date   : 2023-06-01
Purpose: Environmental Association updated for Plot-Level Values
"""
import os
import argparse
import pandas as pd
import glob
import json
import utm
from scipy.spatial.distance import cdist
import numpy as np
from multiprocessing import Pool
from datetime import datetime, timedelta
import re
import subprocess as sp
import shutil
import shlex
import math
import geopandas as gpd
from pathlib import Path
import time


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Environmental Association',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-o',
                        '--out_dir',
                        help='Output directory',
                        required=False,
                        default='environmental_association')

    parser.add_argument('-s',
                        '--season',
                        help='Season during which data were collected',
                        type=str,
                        choices=['10', '11', '12', '13', '14', '15', '16', '17', '18', '19'],
                        required=True)
    
    parser.add_argument('-c',
                        '--crop',
                        help='Crop name of data to download',
                        type=str,
                        choices=['sorghum', 'lettuce', 'cotton', 'soybean', 'sunflower', 'tepary', 'NA'],
                        required=True)

    parser.add_argument('-lev',
                        '--level',
                        help='Data level to download. Choices are 0, 1, 2, 3, or 4.',
                        metavar='str',
                        type=str,
                        default='1')

    parser.add_argument('-i',
                        '--instrument',
                        help='Instrument (sensor) used to collect phenotype data.',
                        type=str,
                        choices=['FLIR', 'PS2'],
                        required=True)

    parser.add_argument('-d',
                        '--data_path',
                        help='Path to directory containing CSV files.',
                        type=str,
                        default=None)
    
    parser.add_argument('-p',
                        '--plot_level',
                        help='Add flag if using level_2 FlirIrCamera plot-level data instead of level_1 individual detection data.',
                        action='store_true')
                
    parser.add_argument('-w',
                        '--weather',
                        help='Which weather station to use (AZMet, EnvironmentLogger, MeteorologicalSensor)',
                        type=str,
                        choices=['AZMet', 'EnvironmentLogger', 'MeteorologicalSensor'],
                        required=True)
    
    parser.add_argument('--nodownload',
                        help='Add flag to skip downloading sensor data in case it is already present.',
                        action='store_true')

    return parser.parse_args()


#-------------------------------------------------------------------------------
def get_date_position(data_path):
    '''
    Downloads gantry raw image data (level_0) and extracts position and timestamps for each data capture/acquisition. 

    Input:
        - data_path: Path containing the raw data (level_0)
    Output: 
        - Dataframe containing time, x position, y position, z positon, latitude, and longitude for each data capture/acquisition
    '''

    # Create empty dataframe and set up counter
    df = pd.DataFrame(columns=['time', 'x_position', 'y_position', 'z_position', 'latitude', 'longitude'])
    cnt = 0

    try:
        # Iterate through each metadata JSON file
        for jfile in glob.glob(data_path):
            # Update counter
            cnt += 1

            # Open JSON file
            with open(jfile) as f:
    
                data = json.load(f)['lemnatec_measurement_metadata']
                
                # Extract time, x gantry position, and y gantry position
                time = pd.to_datetime(data['gantry_system_variable_metadata']['time'])
                x_position = float(data['gantry_system_variable_metadata']['position x [m]'])
                y_position = float(data['gantry_system_variable_metadata']['position y [m]'])
                z_position = float(data['gantry_system_variable_metadata']['position z [m]'])

                # Apply offsets
                # Hardcoded values implemented 2020; Not valid for seasons before S11
                offset_x = -1.035
                offset_y = 1.684
                offset_z = 0.856

                x_position = x_position + offset_x
                y_position = y_position + offset_y              
                z_position = z_position + offset_z

                # Convert gantry coordinate to latitude, longitude
                lat, lon = scanalyzer_to_latlon(x_position, y_position)

                # Save data to a dataframe
                df.loc[cnt] = [time, x_position, y_position, z_position, lat, lon]

        # Sort dataframe by time
        df = df.sort_values('time')

        # Add capture sequence numbers
        df['capture_sequence'] = df['time'].argsort()

    except:
        pass
    
    return df


#-------------------------------------------------------------------------------
def utm_to_latlon(utm_x, utm_y):
    '''
    Convert coordinates from UTM 12N to lat/lon

    Input:
        - utm_x: UTM coordinate for gantry x position
        - utm_y: UTM coordinate for gantry y position
    Output: 
        - Latitude and longitude for the provided x, y gantry positions
    '''

    # Get UTM information from southeast corner of field
    SE_utm = utm.from_latlon(33.07451869, -111.97477775)
    utm_zone = SE_utm[2]
    utm_num  = SE_utm[3]

    return utm.to_latlon(utm_x, utm_y, utm_zone, utm_num)


#-------------------------------------------------------------------------------
def scanalyzer_to_latlon(gantry_x, gantry_y):
    '''
    Convert coordinates from gantry to lat/lon

    Input:
        - gantry_x: Raw gantry x position
        - gantry_y: Raw gantry y position
    Output: 
        - Latitude and longitude for the provided x, y gantry positions
    '''

    utm_x, utm_y = scanalyzer_to_utm(gantry_x, gantry_y)
    return utm_to_latlon(utm_x, utm_y)


#-------------------------------------------------------------------------------
def scanalyzer_to_utm(gantry_x, gantry_y):
    '''
    Convert coordinates from gantry to UTM 12N
    
    Input:
        - gantry_x: Raw gantry x position
        - gantry_y: Raw gantry y position
    Output: 
        - Easting and northing for the provided x, y gantry positions
    '''

    # TODO: Hard-coded
    # Linear transformation coefficients
    ay = 3659974.971; by = 1.0002; cy = 0.0078;
    ax = 409012.2032; bx = 0.009; cx = - 0.9986;

    utm_x = ax + (bx * gantry_x) + (cx * gantry_y)
    utm_y = ay + (by * gantry_x) + (cy * gantry_y)

    return utm_x, utm_y


#-------------------------------------------------------------------------------
def get_phenotype_df(df, data_path, data_type, season, crop):
    '''
    Get phenotype CSV, either thermal or PS2, from processed data (level_1) data on the CyVerse Data Store. 
    
    Input:
        - df: Dataframe containing the gantry time and positions
        - data_path: Path containing the processed data (level_1)
    Output: 
        - Merged dataframe containing the phenotype data in addition to the gantry time and positions in df
    '''

    pheno_df = pd.read_csv(data_path)
    # print(pheno_df)

    if data_type == 'PS2':
        geojson_path = download_geojson(season=season, crop=crop)
        # Open GeoJSON
        gdf = gpd.read_file(geojson_path).drop('plot', axis=1, errors='ignore').rename(columns={'ID': 'plot'})
        print(gdf)

        # Open phenotype information
        pheno_df = pheno_df.drop('Unnamed: 0', axis=1, errors='ignore').rename(columns={'Plot': 'plot'})
        print(pheno_df)

        # Ensure column types match
        gdf['plot'] = gdf['plot'].astype('str').str.zfill(4)
        pheno_df['plot'] = pheno_df['plot'].astype('str').str.zfill(4)

        # Split the 'center_point' into 'lat' and 'lon'
        gdf['lon'] = gdf['geometry'].centroid.x
        gdf['lat'] = gdf['geometry'].centroid.y

        # Extract specific columns
        gdf = gdf[['plot', 'lon', 'lat']]

        # Convert gdf to pandas dataframe
        gdf = pd.DataFrame(gdf)

        # Merge gdf and df
        pheno_df = gdf.merge(pheno_df, on='plot')

    # assuming df and pheno_df are pandas dataframes
    df_coords = df[['latitude', 'longitude']].to_numpy()
    pheno_coords = pheno_df[['lat', 'lon']].to_numpy()

    # calculate pairwise distances
    distances = cdist(pheno_coords, df_coords)

    # find index of minimum distance for each row
    min_indices = np.argmin(distances, axis=1)

    # get closest points from df
    closest_points = df.iloc[min_indices]

    # assuming df and pheno_df are pandas dataframes
    # and min_indices is the array of indices of closest points in df
    pheno_df['time'] = df['time'].iloc[min_indices].values
    pheno_df['x_position'] = df['x_position'].iloc[min_indices].values
    pheno_df['y_position'] = df['y_position'].iloc[min_indices].values
    pheno_df['z_position'] = df['z_position'].iloc[min_indices].values
    pheno_df['capture_sequence'] = df['capture_sequence'].iloc[min_indices].values

    return pheno_df.sort_values('time')


#-------------------------------------------------------------------------------
def get_phenotype_df_plot(df, data_path):
    '''
    Get phenotype CSV for plot-level thermal from processed data (level_2).

    This function:
      - attaches plot identity to each FLIR record
      - preserves capture-level timestamps
      - computes plot-level representative time (t_plot)
      - does NOT collapse rows or overwrite time

    Input:
        - df: metadata dataframe (not used for matching here, only retained for compatibility)
        - data_path: path to plot_thresholding_results.csv

    Output:
        - pheno_df with added columns:
            plot, plot_lat, plot_lon, t_plot
    '''

    # ------------------------------------------------------------------
    # 1. Load FLIR phenotype data (already plot-level FLIR points)
    # ------------------------------------------------------------------
    pheno_df = pd.read_csv(data_path)

    # Ensure time column is datetime
    if 'time' in pheno_df.columns:
        pheno_df['time'] = pd.to_datetime(pheno_df['time'])
    else:
        raise ValueError("Expected 'time' column in plot-level FLIR phenotype CSV")

    # ------------------------------------------------------------------
    # 2. Load plot geometry
    # ------------------------------------------------------------------
    geojson_path = download_geojson(season=season, crop=crop)

    plots_gdf = (
        gpd.read_file(geojson_path)
        .drop('plot', axis=1, errors='ignore')
        .rename(columns={'ID': 'plot'})
    )

    plots_gdf['plot'] = plots_gdf['plot'].astype(str).str.zfill(4)
    plots_gdf['plot_lon'] = plots_gdf.geometry.centroid.x
    plots_gdf['plot_lat'] = plots_gdf.geometry.centroid.y

    # ------------------------------------------------------------------
    # 3. Convert FLIR points to GeoDataFrame
    #    (expects center_lat / center_lon from thresholding results)
    # ------------------------------------------------------------------
    if not {'center_lat', 'center_lon'}.issubset(pheno_df.columns):
        raise ValueError("Expected 'center_lat' and 'center_lon' in plot-level FLIR CSV")

    pheno_gdf = gpd.GeoDataFrame(
        pheno_df,
        geometry=gpd.points_from_xy(pheno_df.center_lon, pheno_df.center_lat),
        crs=plots_gdf.crs
    )

    # ------------------------------------------------------------------
    # 4. Spatial join: assign plot to each FLIR measurement
    # ------------------------------------------------------------------
    pheno_gdf = gpd.sjoin(
        pheno_gdf,
        plots_gdf[['plot', 'plot_lat', 'plot_lon', 'geometry']],
        how='inner',
        predicate='within'
    )

    # drop geometry after spatial join to keep dataframe lightweight
    pheno_df = pd.DataFrame(pheno_gdf.drop(columns='geometry'))

    # ------------------------------------------------------------------
    # 5. Compute plot-level median timestamp (t_plot)
    #    and broadcast back to each row
    # ------------------------------------------------------------------
    t_plot_df = (
        pheno_df
        .groupby('plot')['time']
        .median()
        .reset_index(name='t_plot')
    )

    pheno_df = pheno_df.merge(t_plot_df, on='plot', how='left')

    # ------------------------------------------------------------------
    # 6. Final ordering (preserve expected behavior)
    # ------------------------------------------------------------------
    pheno_df = pheno_df.sort_values('time').reset_index(drop=True)

    return pheno_df


#-------------------------------------------------------------------------------
def associate_weather_spatiotemporal(pheno_df, weather_df, window_minutes=3):
    """
    pheno_df: output of get_phenotype_df_plot
    weather_df: output of get_environment_df

    Returns:
        pheno_df with added weather columns
    """
    weather_df = weather_df.copy()
    weather_df['time'] = pd.to_datetime(weather_df['time'])

    results = {}

    window = pd.Timedelta(minutes=window_minutes)

    for plot, group in pheno_df.groupby('plot'):
        t_plot = group['t_plot'].iloc[0]
        plot_lat = group['plot_lat'].iloc[0]
        plot_lon = group['plot_lon'].iloc[0]

        # --------------------------------------------------
        # 1. Temporal filter
        # --------------------------------------------------
        dt = (weather_df['time'] - t_plot).abs()
        candidates = weather_df[dt <= window]

        if not candidates.empty:
            # --------------------------------------------------
            # 2a. Spatial nearest neighbor
            # --------------------------------------------------
            plot_xy = np.array([[plot_lon, plot_lat]])
            weather_xy = candidates[['longitude', 'latitude']].to_numpy()

            dists = cdist(plot_xy, weather_xy)
            idx = dists.argmin()

            wx = candidates.iloc[idx]
            results[plot] = {
                'air_temp': wx['air_temp'],
                'rh': wx['relative_humidity'],
                'weather_time': wx['time'],
                'weather_match_type': 'spatiotemporal'
            }
        else:
            # --------------------------------------------------
            # 2b. Temporal interpolation fallback
            # --------------------------------------------------
            results[plot] = {
                'air_temp': (
                    weather_df
                    .set_index('time')['air_temp']
                    .sort_index()
                    .interpolate(method='time')
                    .get(t_plot, np.nan)
                ),
                'rh': (
                    weather_df
                    .set_index('time')['relative_humidity']
                    .sort_index()
                    .interpolate(method='time')
                    .get(t_plot, np.nan)
                ),
                'weather_time': t_plot,
                'weather_match_type': 'temporal_interpolation'
            }

    # ------------------------------------------------------
    # 3. Broadcast results back to all rows
    # ------------------------------------------------------
    weather_out = (
        pd.DataFrame.from_dict(results, orient='index')
        .reset_index()
        .rename(columns={'index': 'plot'})
    )

    pheno_df = pheno_df.merge(weather_out, on='plot', how='left')

    return pheno_df


#-------------------------------------------------------------------------------
def process_file(jfile):
    '''
    Processes a single Environmental Logger JSON file, correcting the timestamp and extracting the necessary environmental parameters. 
    
    Input:
        - jfile: Path to a single Environmental Logger JSON file
    Output: 
        - Dataframe containing corrected timestamp formats and environmental parameters from "environment_sensor_readings," including "timestamp," 
          "weather_station," and "sensor par." The "timestamp" value is in format yyyy.MM.dd-HH:mm:ss, which will be converted in subsequent functions. 
          The "weather_station values include sunDirection (degrees), airPressure (hPa), brightness (kilo Lux), relHumidity (relHumPerCent), temperature
          (DegCelsius), windDirection (degrees), precipitation (mm/h), windVelocity (m/s). The "sensor par" value includes photosynthetic active radiation (umol/(m^2*s)).  
    '''
        
    dfs = []

    try:
        with open(jfile) as f:
            print(f"Opening {jfile}...", flush=True)
            data = json.load(f)
            for item in data['environment_sensor_readings']:
                # Convert to appropriate datetime format
                date = pd.to_datetime(item['timestamp'], format="%Y.%m.%d-%H:%M:%S")
                
                # Create a dataframe from the data
                data = {key: float(value['value']) for key, value in item['weather_station'].items()}
                df = pd.DataFrame(data, index=[0])

                # Add datetime to the dataframe
                df['time'] = date

                # Add PAR to df 
                df['par'] = float(item['sensor par']['value'])
                
                print(f"appending df...", flush=True)
                dfs.append(df)

    except:
        print("An error occurred while reading the JSON file or processing the data.")
        dfs = []

    return dfs


#-------------------------------------------------------------------------------
def get_environment_df(data_path):
    """
    Serial version of the EnvironmentLogger extractor.
    Runs process_file() one file at a time for maximum transparency.
    """

    file_list = glob.glob(data_path)
    print(f"[EnvironmentLogger] Found {len(file_list)} files to process.", flush=True)

    all_dfs = []

    for idx, jfile in enumerate(file_list, start=1):
        print(f"[EnvironmentLogger] ({idx}/{len(file_list)}) Processing: {jfile}", flush=True)

        try:
            dfs = process_file(jfile)   # process_file returns a list of dfs
            num_rows = sum(len(df) for df in dfs)
            print(f"[EnvironmentLogger]     -> Extracted {num_rows} rows from file.", flush=True)

            all_dfs.extend(dfs)

        except Exception as e:
            print(f"[EnvironmentLogger] ERROR while processing {jfile}: {e}", flush=True)
            continue

    if not all_dfs:
        raise RuntimeError("No EnvironmentLogger dataframes extracted from any file.")

    print("[EnvironmentLogger] Combining all dataframes...", flush=True)
    env_df = pd.concat(all_dfs, ignore_index=True)
    print("[EnvironmentLogger] Combination complete.", flush=True)

    return env_df.sort_values("time")


#-------------------------------------------------------------------------------
def get_env_dates(date_string):
    '''
    Given a date string, the function finds two flanking dates (+/- 1 day). 
    
    Input:
        - date_string: String containing the date in yyyy-MM-dd format
    Output: 
        - List containing the date from date_string variable in addition to the two flanking dates (+/- 1 day)
    '''

    match = re.search(r'\d{4}-\d{2}-\d{2}', date_string)

    if match:
        date_str = match.group()
        date = datetime.strptime(date_str, '%Y-%m-%d').date()
        day_before = date - timedelta(days=1)
        day_after = date + timedelta(days=1)
    
        return [day_before, date, day_after]
    

#-------------------------------------------------------------------------------
def get_dict():
    '''
    Provides notation for CyVerse directories. 
    
    Input:
        - NA
    Output: 
        - A dictionary containing season, level, and sensor notations which will be used to query the CyVerse Data Store. 
    '''

    irods_dict = {
        'server_path': '/iplant/home/shared/phytooracle/',

        'season': {
            '10': 'season_10_lettuce_yr_2020',
            '11': 'season_11_sorghum_yr_2020',
            '12': 'season_12_sorghum_soybean_sunflower_tepary_yr_2021',
            '13': 'season_13_lettuce_yr_2022',
            '14': 'season_14_sorghum_yr_2022',
            '15': 'season_15_lettuce_yr_2022',
            '16': 'season_16_sorghum_yr_2023',
            '17': 'season_17_lettuce_yr_2023',
            '18': 'season_18_sorghum_yr_2024',
            '19': 'season_19_sorghum_cotton_yr_2025'
        },

        'level': {
            '0': 'level_0', 
            '1': 'level_1',
            '2': 'level_2',
            '3': 'level_3',
            '4': 'level_4'
        },

        'sensor': {
            'FLIR': 'flirIrCamera',
            'PS2': 'ps2Top',
            'RGB': 'stereoTop',
            '3D': 'scanner3DTop',
            'ENV': 'EnvironmentLogger',
            'MET': 'MeteorologicalSensor'
        }
    }

    return irods_dict


#-------------------------------------------------------------------------------
def get_file_list(data_path, sequence):
    '''
    Using the dictionary containing season, level, and sensor notations, this function finds all files matching the season, 
    level, and sensor paths, as well as an identifying sequence such as %.tar.gz. The % is similar to Linux's wild card "*"
    
    Input:
        - data_path: The CyVerse Data Store path created from dictionary
        - sequence: An identifying sequence, such as "%.tar.gz". The "%" character is similar to Linux's wild card "*" character.
    Output: 
        - List of files matching the season, level, sensor, and sequence
    '''
    result = sp.run(f'ilocate {os.path.join(data_path, "%", f"{sequence}")}', stdout=sp.PIPE, shell=True)
    files = result.stdout.decode('utf-8').split('\n')

    # Filter out ilocate errors and empty lines
    files = [f for f in files if f and not f.strip().startswith("ERROR:")]

    return files


#-------------------------------------------------------------------------------
def download_AZMet_rh(date_or_year: str):
    """
    Download AZMet raw hourly 'rh.txt' file for station 06 (Maricopa) based on the provided date/year string.
    - Accepts 'yyyy' or 'yyyy-MM-dd' and extracts the year.
    - Uses current path logic for years >= 2023; legacy logic otherwise.
    - Always saves to the current working directory as 'azmet_06<YY>_rh.txt'.
    - Uses wget via subprocess; raises RuntimeError on failure.

    Returns:
        Path: The resolved path to the downloaded file in the current directory.
    """
    # Extract year from 'yyyy' or 'yyyy-MM-dd'
    year = date_or_year[:4]
    try:
        yr_num = int(year)
    except ValueError:
        raise ValueError(f"Invalid year/date string: {date_or_year!r}. Expected 'yyyy' or 'yyyy-MM-dd'.")

    yy = year[-2:]  # last two digits
    filename = f"azmet_06{yy}_rh.txt"
    out_path = os.path.join(os.getcwd(), filename)

    # Path logic (per your spec: only the base differs)
    if yr_num >= 2023:
        url_base = "https://azmet.arizona.edu/azmet/data/06"
    else:
        if yr_num < 2003:
            raise ValueError("This application does not support AZMet data for years before 2003.")
        else:
            url_base = "https://cales.arizona.edu/azmet/data/06"  # legacy base

    url = f"{url_base}{yy}rh.txt"  # e.g., https://.../0625rh.txt

    # Build wget command (quiet + retries + timeout)
    cmd = [
        "wget",
        "-O", str(out_path),
        "-q",              # quiet; remove this flag to see progress
        "--tries=3",
        "--timeout=30",
        url,
    ]

    try:
        sp.run(cmd, check=True, capture_output=True, text=True)
        return out_path
    except FileNotFoundError:
        raise RuntimeError("wget not found on PATH. Please install wget and try again.")
    except sp.CalledProcessError as e:
        stderr = e.stderr or ""
        raise RuntimeError(
            f"wget failed (exit code {e.returncode}).\n"
            f"URL: {url}\n"
            f"STDERR: {stderr}"
        )


#-------------------------------------------------------------------------------
def download_MeteorologicalSensor_csv(season: str, crop: str, out_dir: str, skip_download: bool = False) -> str:
    irods_dict = get_dict()
    out_base = os.path.join(out_dir, irods_dict['season'][season], irods_dict['sensor']['MET'])
    os.makedirs(out_base, exist_ok=True)

    sequence = "all_sensors_long_merged.csv"
    
    # Check if file already exists
    existing = glob.glob(os.path.join(out_base, "**", sequence), recursive=True)
    if existing:
        print(f"Using existing MeteorologicalSensor CSV: {existing[0]}")
        wait_for_stable_file(existing[0], min_stable_secs=3, timeout=180)
        return existing[0]

    if skip_download:
        raise FileNotFoundError(f"--nodownload specified and file not found under {out_base}")

    # Perform download
    data_path = os.path.join(
        irods_dict['server_path'],
        irods_dict['season'][season],
        irods_dict['level']['2'],
        irods_dict['sensor']['MET'],
        crop
    )
    files = get_file_list(data_path, sequence)
    if not files:
        raise FileNotFoundError(f"No MeteorologicalSensor file found for season {season}, crop {crop}")
    
    twd = os.getcwd()
    download_files(item=files[0], out_path=out_base)
    os.chdir(twd)
    
    # Locate the file after download - Retry logic for locating the file after download
    candidates = []
    for attempt in range(20):  # ~10 seconds total
        candidates = glob.glob(os.path.join(out_base, "**", sequence), recursive=True)
        if candidates:
            break
        time.sleep(0.5)

    if not candidates:
        print(f"DEBUG: Contents of {out_base}: {os.listdir(out_base)}")
        raise FileNotFoundError(f"Downloaded MeteorologicalSensor file not found under {out_base}")
    
    wait_for_stable_file(candidates[0], min_stable_secs=3, timeout=180)
    return candidates[0]


#-------------------------------------------------------------------------------
def download_files(item, out_path):
    '''
    Uses iRODS to access the CyVerse Data Store. The function downloads data and extracts contents from ".tar" and "tar.gz" files if applicable.
    
    Input:
        - item: The list of CyVerse file paths to download locally.
        - out_path: Output directory where data will be saved. 
    Output: 
        - Downloaded and extracted data within the specified output directory 
    '''
        
    os.chdir(out_path)

    if not item or item.strip().startswith("ERROR:"):
        print(f"Skipping invalid item: {item}")
        return
    
    if 'dep' in item:
        print(f"Skipping deprecated item: {item}")
        return

    try:
        item = os.path.normpath(item)

        # Extract date string from item path (two formats supported)
        try:
            match_str = re.search(r'\d{4}\-\d{2}\-\d{2}\_\_\d{2}\-\d{2}\-\d{2}\-\d{3}', item)
            if match_str:
                date = match_str.group()
            else:
                # Try simpler date format
                match_str = re.search(r'\d{4}\-\d{2}\-\d{2}', item)
                if match_str:
                    date = datetime.strptime(match_str.group(), '%Y-%m-%d').date()
                else:
                    # No date found, raise a warning and ignore the date
                    print("No date found")
                    date = "unknown_date"
        except Exception as e:
            print(f"Could not extract date from {item}: {e}")
            date = "unknown_date"
            
        date = str(date)

        print(f"Found item {item}.")
        if not os.path.isdir(date):
            print(f"Making directory {date}.")
            os.makedirs(date)

        # Shell-safe quoting
        q_item = shlex.quote(os.path.basename(item))  # archive name after iget
        q_irods_item = shlex.quote(item)              # full iRODS source path
        q_date = shlex.quote(date)

        if '.tar.gz' in item:
            print(f"Downloading {item}.")
            ret = sp.call(f'iget -KPVT {q_irods_item}', shell=True)
            if ret != 0:
                print(f"iget failed for {item} (ret={ret}); skipping.")
                return

            print(f"Extracting {item}.")
            ret = sp.call(f'tar -xzvf {q_item} -C {q_date}', shell=True)
            if ret != 0:
                print(f"Reattempting to extract {item} without -z.")
                ret2 = sp.call(f'tar -xvf {q_item} -C {q_date}', shell=True)
                if ret2 != 0:
                    print(f"Extraction failed for {item}; skipping.")
            sp.call(f'rm -f {q_item}', shell=True)

        elif '.tar' in item:
            print(f"Downloading {item}.")
            ret = sp.call(f'iget -KPVT {q_irods_item}', shell=True)
            if ret != 0:
                print(f"iget failed for {item} (ret={ret}); skipping.")
                return

            print(f"Extracting {item}.")
            ret = sp.call(f'tar -xvf {q_item} -C {q_date}', shell=True)
            if ret != 0:
                print(f"Extraction failed for {item}; skipping.")
            sp.call(f'rm -f {q_item}', shell=True)

        else:
            # Non-archive: download directly into the date folder
            os.chdir(date)
            ret = sp.call(f'iget -KPVT {q_irods_item}', shell=True)
            if ret != 0:
                print(f"iget failed for {item} (ret={ret}); skipping.")
                return

    except Exception as e:
        # Don’t crash the whole run; just log and continue
        print(f"Failure during download_files for {item}: {e}")


#-------------------------------------------------------------------------------
def download_data(crop, season, level, sensor, sequence, cwd, outdir, download=True, download_first=False):
    '''
    Recursively runs `download_files` to download all data into a single output directory specified by the user.
    
    Input:
        - crop: The name of the crop data to download
        - season: The season numer to download
        - level: The level of data to download
        - sensor: The name of the sensor to download
        - sequence: The identifying sequence to download, such as ".tar" or ".tar.gz"
        - cwd: The current working directory
        - outdir: The output directory
        - download: Boolean value to specify whether to download data (True) or not (False)

    Output: 
        - Downloaded and extracted data within the specified output directory 
    '''

    try:
        irods_dict = get_dict()
        # Create iRODS path from components. 
        data_path = os.path.join(irods_dict['server_path'], irods_dict['season'][season], irods_dict['level'][level], irods_dict['sensor'][sensor])
        if crop != "NA":
            data_path = os.path.join(irods_dict['server_path'], irods_dict['season'][season], irods_dict['level'][level], irods_dict['sensor'][sensor], crop)
        # Get list of all files that match a character sequence.
        print(f'Searching for files matching "{os.path.join(data_path, sequence)}". Note: This process may take 1-5 minutes.')
        files = get_file_list(data_path, sequence)
        print('Matches obtained.')

        # Prepare to download data.
        out_path = os.path.join(outdir, irods_dict['season'][season], irods_dict['sensor'][sensor])
        if not os.path.isdir(out_path):
            os.makedirs(out_path)

        if download:
            os.chdir(out_path)

            if download_first:
                files = files[:1]

            # Download files.
            for item in files: 
                print(f'Downloading {item}.')
                download_files(item=item, out_path=os.path.join(cwd, out_path))

        return out_path
    
    except Exception as e:
        print(f"An error occurred while downloading data: {e}")

        
#-------------------------------------------------------------------------------
def get_vapor_pressure_deficit(air_temp: float, canopy_temp: float, relative_humidity: float) -> float:
    """
    Calculate leaf vapor pressure deficit (VPD) given atmospheric temperature, leaf temperature, and relative humidity.
    
    :param air_temp: The temperature of the atmosphere in degrees Celsius.
    :param canopy_temp: The temperature of the leaf in degrees Celsius.
    :param relative_humidity: The relative humidity as a percentage (0-100).
    :return: The calculated VPD in kPa.
    """
    # Convert temperatures to Kelvin
    canopy_temp += 273.15
    air_temp += 273.15
    
    # Calculate saturation vapor pressures using the Tetens equation
    if canopy_temp < 273.15:
        # Use different coefficients for temperatures below freezing
        svp_leaf = 0.61078 * math.exp((21.87 * (canopy_temp - 273.15)) / (canopy_temp - 7.66))
    else:
        svp_leaf = 0.61078 * math.exp((17.27 * (canopy_temp - 273.15)) / (canopy_temp - 35.85))
    
    if air_temp < 273.15:
        # Use different coefficients for temperatures below freezing
        svp_air = 0.61078 * math.exp((21.87 * (air_temp - 273.15)) / (air_temp - 7.66))
    else:
        svp_air = 0.61078 * math.exp((17.27 * (air_temp - 273.15)) / (air_temp - 35.85))
    
    # Calculate actual vapor pressure
    avp_air = svp_air * (relative_humidity / 100)
    
    # Calculate VPD
    vpd = svp_leaf - avp_air
    
    return vpd


#-------------------------------------------------------------------------------
def get_geojson_path(season, crop):
    geojson_dict = {
        '10': {
            'lettuce': '/iplant/home/shared/phytooracle/season_10_lettuce_yr_2020/level_0/season10_multi_latlon_geno.geojson'
        },
        '11': {
            'sorghum': '/iplant/home/shared/phytooracle/season_11_sorghum_yr_2020/level_0/season11_multi_latlon_geno.geojson'
        },
        '12': {
            'sorghum': '/iplant/home/shared/phytooracle/season_12_sorghum_soybean_sunflower_tepary_yr_2021/level_0/season12_multi_latlon_geno_updated.geojson',
            'soybean': '/iplant/home/shared/phytooracle/season_12_sorghum_soybean_sunflower_tepary_yr_2021/level_0/season12_multi_latlon_geno_updated.geojson',
            'sunflower': '/iplant/home/shared/phytooracle/season_12_sorghum_soybean_sunflower_tepary_yr_2021/level_0/season12_multi_latlon_geno_updated.geojson',
            'tepary': '/iplant/home/shared/phytooracle/season_12_sorghum_soybean_sunflower_tepary_yr_2021/level_0/season12_multi_latlon_geno_updated.geojson'
        },
        '13': {
            'lettuce': '/iplant/home/shared/phytooracle/season_13_lettuce_yr_2022/level_0/season13_multi_latlon_geno.geojson'
        },
        '14': {
            'sorghum': '/iplant/home/shared/phytooracle/season_14_sorghum_yr_2022/level_0/season14_multi_latlon_geno_correction_labeled.geojson'
        },
        '15': {
            'lettuce': '/iplant/home/shared/phytooracle/season_15_lettuce_yr_2022/level_0/season15_multi_latlon_geno.geojson'
        },
        '16': {
            'sorghum': '/iplant/home/shared/phytooracle/season_16_sorghum_yr_2023/level_0/season16_multi_latlon_geno_correction_relabeled.geojson'
        },
        '17': {
            'lettuce': '/iplant/home/shared/phytooracle/season_17_lettuce_yr_2023/level_0/season17_multi_latlon_geno.geojson'
        },
        '18': {
            'sorghum': '/iplant/home/shared/phytooracle/season_18_sorghum_yr_2024/level_0/season18_multi_latlon_geno_sorghum_v2.geojson',
            'cotton': '/iplant/home/shared/phytooracle/season_18_sorghum_yr_2024/level_0/season18_multi_latlon_geno_cotton_v2.geojson'
        },
        '19': {
            'sorghum': '/iplant/home/shared/phytooracle/season_19_sorghum_cotton_yr_2025/level_0/season19_multi_latlon_geno_sorghum.geojson',
            'cotton': '/iplant/home/shared/phytooracle/season_19_sorghum_cotton_yr_2025/level_0/season19_multi_latlon_geno_cotton.geojson'
        }
    }

    try:
        return geojson_dict[season][crop]
    except KeyError:
        raise ValueError(f"No geojson path found for season '{season}' and crop '{crop}'")


#-------------------------------------------------------------------------------
def download_geojson(season, crop):
    irods_path = get_geojson_path(season=season, crop=crop)
    sp.call(f'iget -fKPVT {irods_path}', shell=True)

    return os.path.basename(irods_path)


#-------------------------------------------------------------------------------
def wait_for_stable_file(path, min_stable_secs=3, timeout=120):
    """Wait until `path` exists, is non-zero, and its size is stable for `min_stable_secs`."""
    start = time.time()
    last_size = -1
    last_change = time.time()

    while time.time() - start < timeout:
        if os.path.isfile(path):
            size = os.path.getsize(path)
            if size > 0:
                if size != last_size:
                    last_size = size
                    last_change = time.time()
                else:
                    if time.time() - last_change >= min_stable_secs:
                        return True
        time.sleep(0.5)

    raise TimeoutError(
        f"File '{path}' did not become stable/readable within {timeout} seconds. "
        f"Last observed size={last_size} bytes.")


#-------------------------------------------------------------------------------
def main():
    args = get_args()

    # Create output directory
    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)

    # Get the working directory
    wd = os.getcwd()

    # Download sensor data
    if args.nodownload:
        irods_dict = get_dict()
        data_path = os.path.join(args.out_dir, irods_dict['season'][args.season], irods_dict['sensor'][args.instrument])
    else: 
        data_path = download_data(
                                crop = "NA",
                                season = args.season,
                                level = '0',
                                sensor = args.instrument,
                                sequence = '%/%.tar' if args.season=='10' else '%/%.tar.gz',
                                cwd = wd,
                                outdir = args.out_dir)
    os.chdir(wd)

    # Find dates for this season
    path_list = [path for path in glob.glob(os.path.join(data_path, '*')) if '2222' not in path]

    # Handle MeteorologicalSensor station
    MeteorologicalSensor_csv_path = None
    if args.weather == 'MeteorologicalSensor':
        print("Downloading MeteorologicalSensor data")
        MeteorologicalSensor_csv_path = download_MeteorologicalSensor_csv(
            season=args.season,
            crop=args.crop,
            out_dir=args.out_dir
        )
        print("MeteorologicalSensor data downloaded")
        
    # Iterate through all dates within this season
    for path in path_list:
        try:
            if args.season == '10':
                date_string  = re.search(r"\b\d{4}-\d{2}-\d{2}\b", path).group() if re.search(r"\b\d{4}-\d{2}-\d{2}\b", path) else None
            else:
                date_string = re.search(r"\b\d{4}-\d{2}-\d{2}__\d{2}-\d{2}-\d{2}-\d{3}\b", path).group() if re.search(r"\b\d{4}-\d{2}-\d{2}__\d{2}-\d{2}-\d{2}-\d{3}\b", path) else None

            print(f'Date string: {date_string}')

            if not args.crop == 'NA':
                date_species = '_'.join([date_string, args.crop])
            else:
                date_species = date_string

            date_list = get_env_dates(date_string = date_species)

            # Download weather data
            if args.weather == 'EnvironmentLogger':
                print("Downloading EnvironmentLogger data")
                for date in date_list:
                    env_path = download_data(
                                    crop = "NA",
                                    season = args.season,
                                    level = '0',
                                    sensor = 'ENV',
                                    sequence = f'{date}.tar.gz',
                                    cwd = wd,
                                    outdir = args.out_dir)
                print("EnvironmentLogger data downloaded")
            elif args.weather == 'AZMet':
                print('Downloading AZMet data')
                # Date string expected to be in format yyyy-MM-dd; parse this to get the correct AZMet hourly data
                AZMet_data = download_AZMet_rh(date_species)
                print("AZMet data downloaded")
            elif args.weather == 'MeteorologicalSensor':
                env_path = MeteorologicalSensor_csv_path
            else:
                raise ValueError(f"Unsupported weather station: {args.weather}.")
            
            os.chdir(wd)

            # Get gantry metadata
            print("Setting meta_df")
            meta_df = get_date_position(data_path = os.path.join(data_path, date_string, '*', '*', '*', '*.json') if args.season == '10' else os.path.join(data_path, date_string, '*', '*', '*.json'))
    
            # Determining the sequence to use based on specified instrument (sensor) name
            if args.instrument == 'PS2':
                sensor_seq = f'{date_species}/%_aggregation_out.tar'
            elif args.instrument == 'FLIR':
                if args.plot_level:
                    sensor_seq = f'{date_species}/%_plot_temps_out.tar'
                else:
                    sensor_seq = f'{date_species}/%_detect_out.tar'
            else:
                raise ValueError(f"Unsupported instrument: {args.instrument}.")
            print("sensor_seq: ", sensor_seq)

            # Download phenotype data
            if args.data_path:
                csv_path = args.data_path
            else:
                print("Downloading phenotype data")
                csv_path = download_data(
                                        crop = args.crop,
                                        season = args.season,
                                        level = args.level,
                                        sensor = args.instrument,
                                        sequence = sensor_seq,
                                        cwd = wd,
                                        outdir = args.out_dir)

            os.chdir(wd)

            # Open phenotype data
            if args.plot_level:
                candidates = glob.glob(os.path.join(csv_path, date_string, '*', '*plot_thresholding_results.csv'))
                if not candidates:
                    print(f"No plot_thresholding_results.csv found for {date_string}; skipping this date.")
                    continue
                print("Setting pheno_df")
                pheno_df = get_phenotype_df_plot(
                    df=meta_df, 
                    data_path=candidates[0]
                    )
            else:
                candidates = glob.glob(os.path.join(csv_path, date_string, '*', '*.csv'))
                if not candidates:
                    print(f"No csv found for {date_string}; skipping this date.")
                    continue
                print("Setting pheno_df")
                pheno_df = get_phenotype_df(
                    df = meta_df, 
                    data_path = candidates[0], 
                    data_type=args.instrument,
                    season=args.season,
                    crop=args.crop
                    )

            if args.plot_level and args.instrument == 'FLIR':
                pheno_df = associate_weather_spatiotemporal(
                    pheno_df=pheno_df,
                    weather_df=env_df,
                    window_minutes=3
                )

            print("Setting env_df")
            if args.weather == "EnvironmentLogger":
                env_df = get_environment_df(data_path = os.path.join(env_path, '*', '*', '*', '*.json') if args.season == '10' else os.path.join(env_path, '*', '*', '*.json'))
            elif args.weather == "AZMet":
                env_df = pd.read_csv(
                    AZMet_data,
                    sep=",",           # or sep=";", sep="\t"
                    comment="#",
                    na_values=["NA", "NaN", "-999", "-99"],
                    header=None,                # no header in file
                    names=["year", "day_of_year", "hour_of_day", "temperature", "relHumidity", "vapor_pressure_deficit_azmet", 
                           "solar_radiation", "precipitation", "4in_soil_temperature", "20in_soil_temperature", "wind_speeg_avg", 
                           "wind_vector_magnitude", "wind_vector_direction", "wind_direction_stdev", "wind_speed_max", 
                           "reference_evapotranspiration_Eto", "actual_vapor_pressure", "dewpoint_hourly_avg"],
                )
                # Add in a time column to more closely match with gantry data
                env_df["time"] = (
                    pd.to_datetime(env_df["year"], format="%Y")                # Jan 1 of year
                    + pd.to_timedelta(env_df["day_of_year"] - 1, unit="D")     # add DOY offset
                    + pd.to_timedelta(env_df["hour_of_day"] % 24, unit="h")    # hour (24 → 0)
                    + pd.to_timedelta((env_df["hour_of_day"] == 24).astype(int), unit="D")  # rollover day if hour=24
                )
            elif args.weather == 'MeteorologicalSensor': 
                env_df = pd.read_csv(
                    env_path,
                    header=0,
                    low_memory=False
                )
                env_df.columns = env_df.columns.str.strip()
                
                required_cols = {"timestamp"}
                missing_cols = required_cols - set(env_df.columns)
                if missing_cols:
                    raise RuntimeError(
                        f"MeteorologicalSensor CSV missing required columns: {missing_cols}. "
                        f"Columns present: {list(env_df.columns)}"
                    )
                
                if "sensor_type" in env_df.columns:
                    env_df = env_df[env_df["sensor_type"] == "weather_station"].copy()

                
                # Add in a time column to more closely match with gantry data
                env_df["time"] = (
                    pd.to_datetime(env_df['timestamp'], format="%Y-%m-%d %H:%M:%S.%f", errors="coerce")
                    .dt.tz_localize("UTC")
                    .dt.tz_convert("America/Phoenix")
                    .dt.tz_localize(None) 
                )
                
                drop_columns = ['NS_tilt_angle', 'EW_tilt_angle', 'plot_id', 'sensor_type', 'tarball_name', 'number_lightnings', 'lightning_distance',
                            'meta_date', 'meta_time', 'meta_crop', 'name', 'date', 'date_int', 'tension', 'value', 'timestamp']
                env_df = env_df.drop(drop_columns, axis=1)

                env_df = env_df.dropna(subset=["time"])

                
                # Rename some columns to match the original gantry data for consistency
                env_df = env_df.rename(columns={"air_temperature": "temperature", "relative_humidity": "relHumidity"})
            else:
                raise ValueError(f"Unsupported weather station: {args.weather}.")

            print("sorting env_df")
            env_df = env_df.sort_values("time")
            if env_df.empty:
                raise RuntimeError(
                    "No valid weather rows after cleaning: timestamp parse failures or all -9999 placeholders.\n"
                    f"Example bad timestamps (first 5): {env_df['timestamp'].head(5).tolist() if 'timestamp' in env_df else 'N/A'}"
                )
            if pheno_df.empty:
                raise RuntimeError(
                    "No valid phenotype rows after cleaning—check CSV selection and time derivation in get_date_position()."
                )

            
            
            print("pheno_df columns:", list(pheno_df.columns))
            print("env_df columns:", list(env_df.columns))
            print("pheno_df time dtype:", pheno_df["time"].dtype if "time" in pheno_df.columns else "missing")
            print("env_df time dtype:", env_df["time"].dtype if "time" in env_df.columns else "missing")
            print("pheno_df time NaT rate:", pheno_df["time"].isna().mean() if "time" in pheno_df.columns else "n/a")
            print("env_df time NaT rate:", env_df["time"].isna().mean() if "time" in env_df.columns else "n/a")
            
            # Merge the phenotype and weather dataframes on the "time" column, finding the closest match in env_df for each row in pheno_df
            print("merging dfs")
            result = pd.merge_asof(pheno_df, env_df, on='time', direction='nearest')

            # Calculate additional columns based on instrument (sensor) type
            if args.instrument == 'FLIR':
                print("calculate additional columns...")
                if args.plot_level:
                    result['canopy_temperature_depression'] = result['temperature'] - result['plot_plant_temp']
                    result['vapor_pressure_deficit'] = result.apply(lambda x: get_vapor_pressure_deficit(x['temperature'], x['plot_plant_temp'], x['relHumidity']), axis=1)
                else:
                    result['canopy_temperature_depression'] = result['temperature'] - result['median']
                    result['vapor_pressure_deficit'] = result.apply(lambda x: get_vapor_pressure_deficit(x['temperature'], x['median'], x['relHumidity']), axis=1)

            # Drop potentially erroneous column from EnvironmentLogger weather station
            if args.weather == 'EnvironmentLogger' and 'brightness' in result.columns:
                result = result.drop('brightness', axis=1)

            # Save CSV to defined output directory
            base_name = '_'.join([date_species, args.instrument, args.weather])
            suffix = 'plot_level_environmental_association.csv' if args.plot_level else 'individual_level_environmental_association.csv'
            out_path = os.path.join(args.out_dir, f'{base_name}_{suffix}')
            result.to_csv(out_path, index=False)

            # Clean up input data
            if args.weather == 'EnvironmentLogger':
                shutil.rmtree(env_path)

        except Exception as e:
            print(f"An error occurred while processing path: {path}: {e}. Continuing with next path.")
            try:
                if 'env_path' in locals() and env_path and os.path.isdir(env_path):
                    shutil.rmtree(env_path)
            except Exception as ce:
                print(f"Cleanup skipped or failed for {env_path}: {ce}")
    
    # Clean up input data
    if args.weather == 'MeteorologicalSensor':
        shutil.rmtree(Path(env_path).resolve().parent)
    if args.weather == 'AZMet':
        p = Path(AZMet_data)
        if p.is_file():
            p.unlink()
        else:
            print(f"Warning: file not found: {p}")


# --------------------------------------------------
if __name__ == '__main__':
    main()
