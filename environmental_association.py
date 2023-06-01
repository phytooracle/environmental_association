#!/usr/bin/env python3
"""
Author : Emmanuel Gonzalez, Jeffrey Demieville
Date   : 2023-06-01
Purpose: Environmental Association
"""

import os
import sys
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
                        choices=['10', '11', '12', '13', '14', '15', '16'],
                        required=True)
    
    parser.add_argument('-c',
                        '--crop',
                        help='Crop name of data to download',
                        type=str,
                        choices=['sorghum', 'lettuce'],
                        required=True)
    
    parser.add_argument('-i',
                        '--instrument',
                        help='Instrument (sensor) used to collect phenotype data.',
                        type=str,
                        choices=['FLIR', 'PS2'],
                        required=True)
    
    return parser.parse_args()


#-------------------------------------------------------------------------------
def get_date_position(data_path):
    '''
    Downloads gantry raw data (level_0) and extracts position and timestamps for each data capture/acquisition. 

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
                offset_x = -1.035
                x_position = x_position + offset_x

                offset_y = 1.684
                y_position = y_position + offset_y
                
                offset_z = 0.856
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
def get_phenotype_df(df, data_path):
    '''
    Get phenotype CSV, either thermal or PS2, from processed data (level_1) data on the CyVerse Data Store. 
    
    Input:
        - df: Dataframe containing the gantry time and positions
        - data_path: Path containing the processed data (level_1)
    Output: 
        - Merged dataframe containing the phenotype data in addition to the gantry time and positions in df
    '''

    pheno_df = pd.read_csv(data_path)

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
    with open(jfile) as f:
        data = json.load(f)
        for item in data['environment_sensor_readings']:
            # Convert to appropriate datetime format
            date = pd.to_datetime(item['timestamp'], format="%Y.%m.%d-%H:%M:%S")
            
            # Create a dataframe from the data
            data = {key: float(value['value']) for key, value in item['weather_station'].items()} #f"{key}_{value['unit']}"
            df = pd.DataFrame(data, index=[0])

            # Add datetime to the dataframe
            df['time'] = date

            # Add PAR to df 
            df['par'] = float(item['sensor par']['value'])
            
            dfs.append(df)
    return dfs


#-------------------------------------------------------------------------------
def get_environment_df(data_path):
    '''
    Uses multiprocessing to run the function `process_file` to extract multiple Environmental Logger JSON files and combine them into a single dataframe. 
    
    Input:
        - data_path: Path containing the raw data (level_0)
    Output: 
        - Merged dataframe containing the timestamp and environmental parameters from multiple Environmental Logger JSON files
    '''
        
    with Pool() as pool:
        results = pool.map(process_file, glob.glob(data_path))
    dfs = [df for result in results for df in result]
    # Combine all dataframes in the list into one
    env_df = pd.concat(dfs, ignore_index=True)
    return env_df.sort_values('time')


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
            '15': 'season_15_lettuce_yr_2022'
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
            'ENV': 'EnvironmentLogger'
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

    return files


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

    if not 'deprecated' in item:

        try:
            item = os.path.normpath(item)

            try:

                match_str = re.search(r'\d{4}-\d{2}-\d{2}__\d{2}-\d{2}-\d{2}-\d{3}', item)
                date = match_str.group()
                # date = datetime.strptime(match_str.group(), '%Y-%m-%d').date()
            except:
                match_str = re.search(r'\d{4}-\d{2}-\d{2}', item)
                date = datetime.strptime(match_str.group(), '%Y-%m-%d').date()
                date = str(date)

            print(f"Found item {item}.")

            if not os.path.isdir(date):
                print(f"Making directory {date}.")
                os.makedirs(date)


            if '.tar.gz' in item: 
                print(f"Downloading {item}.")
                sp.call(f'iget -KPVT {item}', shell=True)

                print(f"Extracting {item}.")
                ret = sp.call(f'tar -xzvf {os.path.basename(item)} -C {date}', shell=True)
                # ret = sp.call(f'tar -c --use-compress-program=pigz -f {os.path.basename(item)}', shell=True) #-C {date} 

                if ret != 0:
                    print(f"Reattempting to extract {item}.")
                    sp.call(f'tar -xvf {os.path.basename(item)} -C {date}', shell=True)

                sp.call(f'rm {os.path.basename(item)}', shell=True)

            elif '.tar' in item:
                print(f"Downloading {item}.")
                sp.call(f'iget -KPVT {item}', shell=True)
                
                print(f"Extracting {item}.")
                sp.call(f'tar -xvf {os.path.basename(item)} -C {date}', shell=True)
                sp.call(f'rm {os.path.basename(item)}', shell=True)

            else:
                os.chdir(date)
                sp.call(f'iget -KPVT {item}', shell=True)
            
        except:
            pass

        
#-------------------------------------------------------------------------------
def download_data(crop, season, level, sensor, sequence, cwd, outdir, download=True):
    '''
    Recursively runs `download_files` to download all data into a single output directory specified by the user.
    
    Input:
        - crop: The name of the crop data to download, either "lettuce" or "sorghum"
        - season: The season numer to download, either 14, 15, or 16
        - level: The level of data to download, either 0, 1, 2
        - sensor: The name of the sensor to download, either RGB, FLIR, PS2, or 3D
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

            # Download files.
            for item in files: 
                print(f'Downloading {item}.')
                download_files(item=item, out_path=os.path.join(cwd, out_path))
                
            # os.chdir(cwd)

        return out_path
    
    except Exception as e:
        # code to handle the exception
        print(f"An error occurred while downloading data: {e}")


# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = get_args()

    # Create output directory
    if not os.path.isdir(args.out_dir):

        os.makedirs(args.out_dir)

    # Get the working directory
    wd = os.getcwd()

    # Download sensor data
    data_path = download_data(
                            crop = "NA",
                            season = args.season,
                            level = '0',
                            sensor = args.instrument,
                            sequence = '%/%.tar.gz',
                            cwd = wd,
                            outdir = args.out_dir)
    os.chdir(wd)

    # Find dates for this season
    path_list = glob.glob(os.path.join(data_path, '*'))

    # Iterate through all dates within this season
    for path in path_list:

        split_path = path.split(os.sep)
        date_string = split_path[-1]
        date_species = '_'.join([date_string, args.crop])
        date_list = get_env_dates(date_string = date_species)

        # Download Environment Logger data
        for date in date_list:

            env_path = download_data(
                            crop = "NA",
                            season = args.season,
                            level = '0',
                            sensor = 'ENV',
                            sequence = f'{date}.tar.gz',
                            cwd = wd,
                            outdir = args.out_dir)
            
            os.chdir(wd)

        # Get gantry metadata
        meta_df = get_date_position(data_path = os.path.join(data_path, date_string, '*', '*', '*.json'))

        # Determining the sequence to use based on specified instrument (sensor) name
        if args.instrument == 'FLIR':

            sensor_seq = f'{date_species}/%_detect_out.tar'
        
        elif args.instrument == 'PS2':

            sensor_seq = f'{date_species}/%_aggregation_out.tar'

        else:

            raise ValueError(f"Unsupported instrument: {args.instrument}.")

        # Download phenotype data
        csv_path = download_data(
                                crop = args.crop,
                                season = args.season,
                                level = '1',
                                sensor = args.instrument,
                                sequence = sensor_seq,
                                cwd = wd,
                                outdir = args.out_dir)
        os.chdir(wd)

        # Open phenotype data
        pheno_df = get_phenotype_df(df = meta_df, data_path = glob.glob(os.path.join(csv_path, date_string, '*', '*.csv'))[0])

        # Open environmental logger data
        env_df = get_environment_df(data_path = os.path.join(env_path, '*', '*', '*.json'))

        # Merge the phenotype and environmental logger dataframes on the "time" column, finding the closest match in env_df for each row in pheno_df
        result = pd.merge_asof(pheno_df, env_df, on='time', direction='nearest')

        # Calculate additional columns based on instrument (sensor) type
        if args.instrument == 'FLIR':

            # Calculate normalized canopy temperature
            result['normalized_temp'] = result['median'] - result['temperature']

        # Drop potentially erroneous column
        result = result.drop('brightness', axis=1)

        # Save CSV to defined output directory
        result.to_csv(os.path.join(args.out_dir, '_'.join([date_species, 'environmental_association.csv'])), index=False)

        # Clean up input data
        shutil.rmtree(env_path)
        shutil.rmtree(data_path)


# --------------------------------------------------
if __name__ == '__main__':

    main()