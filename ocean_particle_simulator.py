import numpy as np
import xarray as xr
import pandas as pd
import sys
import os
import matplotlib.pyplot as plt
import json
import warnings
warnings.filterwarnings('ignore')
from cartopy import config
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.path as mpath
import csv
import math
import random
import logging
import argparse
from subprocess import run

# Local NASA ECCO python package.
# By default ECCOv4-py is in home directory
# If not then specify in PYTHONPATH environment variable
# https://ecco-v4-python-tutorial.readthedocs.io/Installing_Python_and_Python_Packages.html
sys.path.append(f'{os.path.expanduser("~")}/ECCOv4-py')
try:
    import ecco_v4_py as ecco
except Exception as e:
    print('Please add ECCOv4-py to System Environment PYTHONPATH!', file=sys.stderr)
    raise(e)

# convert velocity from meter/second to degree/month
MPS_TO_DEG_PER_MONTH = (365.0/12) * 24 * 3600 / (40075017.0/360)

# Average sinking speed, assuming 4 months per layer
KVEL = 0.25

# index for next particle, TODO: concurrency
particle_id = 0
def next_particle_id():
    global particle_id
    index = particle_id
    particle_id += 1
    return index

# Specify where ECCO data set is. By default it's in ~/eccodata
# https://ecco-v4-python-tutorial.readthedocs.io/Downloading_the_ECCO_v4_state_estimate.html
# Currently use Version 4 Release 4
def configure_base_dir(base_dir=None):
    eccodata_dir = ''
    if base_dir and os.path.isdir(base_dir):
        eccodata_dir = base_dir
    elif 'ECCODATA_DIR' in os.environ and os.path.isdir(os.path.expanduser(os.environ['ECCODATA_DIR'])):
        eccodata_dir = os.path.expanduser(os.environ['ECCODATA_DIR'])
    elif os.path.isdir(os.path.expanduser('~/eccodata')):
        eccodata_dir = os.path.expanduser('~/eccodata')
    elif os.path.isdir('/eccodata'):
        eccodata_dir = '/eccodata'
    else:
        raise Exception('Cannot find eccodata directory')
    logging.info(f'Setting eccodata to {eccodata_dir}')
    return eccodata_dir


# Load ECCO dataset to memory by the year
# This is the most memory and CPU intensive operations, call it as less as possible
def load_ecco_ds(year, base_dir, vars=['UVEL', 'VVEL']):
    # ECCO_dir = base_dir + '/Version4/Release3_alt'
    ECCO_dir = base_dir + '/Version4/Release4'
    grid_dir = f'{ECCO_dir}/nctiles_grid'
    logging.info(f'Loading {year} {vars}')
    # ecco_grid = ecco.load_ecco_grid_nc(grid_dir, 'ECCOv4r3_grid.nc')
    ecco_grid = ecco.load_ecco_grid_nc(grid_dir, 'ECCO-GRID.nc')
    day_mean_dir = f'{ECCO_dir}/nctiles_monthly'
    ecco_vars = ecco.recursive_load_ecco_var_from_years_nc(
        day_mean_dir, vars_to_load=vars, years_to_load=year , dask_chunk=False
        )
    ecco_ds = xr.merge((ecco_grid , ecco_vars))
    return ecco_ds


def usage():
    parser = argparse.ArgumentParser(description='Compute particle movements within one tile')
    parser.add_argument('inputfile')
    # TODO k_base1, month_base1
    parser.add_argument('--kvel', type=float, default=0.25, help='average sinking speed (layer/month)')
    parser.add_argument('--from-year', type=int, default=1992, help='Starting year')
    parser.add_argument('--to-year', type=int, default=2017, help='End year')
    parser.add_argument('--fudge-pct', type=int, default=30, help='Percentage of factor to disturb')
    parser.add_argument('--debug', action='store_true', help='Debug Mode')
    parser.add_argument('--test', action='store_true', help='Test Mode')
    action_group = parser.add_mutually_exclusive_group()
    action_group.add_argument('--only-plot', action='store_true', help='Only Plot')
    action_group.add_argument('--only-compute', action='store_true', help='Only output results csv')
    parser.add_argument('--png-ym', help='year:month')
    parser.add_argument('--keep-png', action='store_true', help='year:month')
    plot_group = parser.add_mutually_exclusive_group()
    plot_group.add_argument('--plot-1tile', type=int, default=10)
    plot_group.add_argument('--plot-all-tiles', action='store_true', help='Plot all tiles by tiles')
    plot_group.add_argument('--plot-all-lonlat', action='store_true', help='Plot all tiles by lon-lat')
    refresh_group = parser.add_mutually_exclusive_group()
    refresh_group.add_argument('--no-refresh', action='store_true', help='No new particle is added')
    refresh_group.add_argument('--refresh-random', action='store_true', help='Refresh random within the tile')
    refresh_group.add_argument('--refresh-original', action='store_true', help='Refresh from original particle')
    args = parser.parse_args()
    return args

####################################################
# Compute Section
####################################################


# Check if (x,y) is within the 90x90 tile.
def outOfTile(x,y):
    if y < 0:
        return True
    elif y >= 90:
        return True
    elif x < 0:
        return True
    elif x >= 90:
        return True
    else:
        return False

# A particle may move from one tile to another,
# this function converts coordinates between tiles
# Return: (newtile, newx, newy)
def adjustTile(tile, ix, jy):
    newtile, newi, newj = tile, ix, jy
    if (ix < 0) and (jy >= 90): # top left
        if tile == 0: #=> 11
            newtile, newi, newj = 11, 180-jy, ix+90
        elif tile == 1: #=>10
            newtile, newi, newj = 10, 180-jy, ix+90
        elif tile == 2: #undefined
            pass
        elif tile == 3: #=>1
            newtile, newi, newj = 1, ix+90, jy-90
        elif tile == 4: #=>2
            newtile, newi, newj = 2, ix+90, jy-90
        elif tile == 5: # undefined
            pass
        elif tile == 6: # undefined
            pass
        elif tile == 7: #undefined
            pass
        elif tile == 8: #=>10
            newtile, newi, newj = 10, ix+90, jy-90
        elif tile == 9: #=>11
            newtile, newi, newj = 11, ix+90, jy-90
        elif tile == 10: #undefined
            pass
        elif tile == 11: #=>2
            newtile, newi, newj = 2, jy-90, -ix
        elif tile == 12: #=>1
            newtile, newi, newj = 1, jy-90, -ix
    elif (0 <= ix < 90) and (jy >= 90): # top
        if tile == 0: #=> 1
            newtile, newi, newj = (1, ix, jy-90)
        elif tile == 1: #=> 2
            newtile, newi, newj = (2, ix, jy-90)
        elif tile == 2: #=> 6
            newtile, newi, newj = (6, jy-90, 90-ix)
        elif tile == 3: #=> 4
            newtile, newi, newj = (4, ix, jy-90)
        elif tile == 4: #=> 5
            newtile, newi, newj = (5, ix, jy-90)
        elif tile == 5: #=> 6
            newtile, newi, newj = (6, ix, jy-90)
        elif tile == 6: #=> 10
            newtile, newi, newj = (10, jy-90, 90-ix)
        elif tile == 7: #=> 10
            newtile, newi, newj = (10, ix, jy-90)
        elif tile == 8: #=> 11
            newtile, newi, newj = (11, ix, jy-90)
        elif tile == 9: #=> 12
            newtile, newi, newj = (12, ix, jy-90)
        elif tile == 10: #=> 2
            newtile, newi, newj = 2, jy-90, 90-ix
        elif tile == 11: #=> 1
            newtile, newi, newj = 1, jy-90, 90-ix
        elif tile == 12: #=> 0
            newtile, newi, newj = 0, jy-90, 90-ix
    elif (ix >= 90) and (jy >= 90): # top right
        if tile == 0: #=>4
            newtile, newi, newj = 4, ix-90, jy-90
        elif tile == 1: #=>5
            newtile, newi, newj = 5, ix-90, jy-90
        elif tile == 2: # undefineable
            pass
        elif tile == 3: #=>8
            newtile, newi, newj = 8, 180-jy, ix-90
        elif tile == 4: #=>7
            newtile, newi, newj = 7, 180-jy, ix-90
        elif tile == 5: # undefineable
            pass
        elif tile == 6: # undefineable
            pass
        elif tile == 7: #=>11
            newtile, newi, newj = 11, ix-90, jy-90
        elif tile == 8: #=>12
            newtile, newi, newj = 12, ix-90, jy-90
        elif tile == 9: # undefineable
            pass
        elif tile == 10: #=>1
            newtile, newi, newj = 1, jy-90, 180-ix
        elif tile == 11: #=>0
            newtile, newi, newj = 0, jy-90, 180-ix
        elif tile == 12: # undefineable
            pass
    elif (ix < 0) and (0 <= jy < 90): # left
        if tile == 0: #=>12
            newtile, newi, newj = 12, 90-jy, 90+ix
        elif tile == 1: #=>11
            newtile, newi, newj = 11, 90-jy, 90+ix
        elif tile == 2: #=> 10
            newtile, newi, newj = 10, 90-jy, 90+ix
        elif tile == 3: #=> 0
            newtile, newi, newj = 0, 90+ix, jy
        elif tile == 4: #=> 1
            newtile, newi, newj = 1, 90+ix, jy
        elif tile == 5: #=> 2
            newtile, newi, newj = 2, 90+ix, jy
        elif tile == 6: #=>2
            newtile, newi, newj = 2, 90-jy, 90+ix
        elif tile == 7: #=> 6
            newtile, newi, newj = 6, 90+ix, jy
        elif tile == 8: #=> 7
            newtile, newi, newj = 7, 90+ix, jy
        elif tile == 9: #=> 8
            newtile, newi, newj = 8, 90+ix, jy
        elif tile == 10: #=>6
            newtile, newi, newj = 6, 90-jy, 90+ix
        elif tile == 11: #=> 10
            newtile, newi, newj = 10, 90+ix, jy
        elif tile == 12: #=> 11
            newtile, newi, newj = 11, 90+ix, jy
    elif (ix >= 90) and (0 <= jy < 90): # right
        if tile == 0: #=>3
            newtile, newi, newj = 3, ix-90, jy
        elif tile == 1: #=>4
            newtile, newi, newj = 4, ix-90, jy
        elif tile == 2: #=> 5
            newtile, newi, newj = 5, ix-90, jy
        elif tile == 3: #=>9
            newtile, newi, newj = 9, 90-jy, ix-90
        elif tile == 4: #=>8
            newtile, newi, newj = 8, 90-jy, ix-90
        elif tile == 5: #=>7
            newtile, newi, newj = 7, 90-jy, ix-90
        elif tile == 6: #=> 7
            newtile, newi, newj = 7, ix-90, jy
        elif tile == 7: #=> 8
            newtile, newi, newj = 8, ix-90, jy
        elif tile == 8: #=> 9
            newtile, newi, newj = 9, ix-90, jy
        elif tile == 9: # undefined
            pass
        elif tile == 10: #=> 11
            newtile, newi, newj = 11, ix-90, jy
        elif tile == 11: #=> 12
            newtile, newi, newj = 12, ix-90, jy
        elif tile == 12: # undefined
            pass
        elif tile == 10: #=>11
            newtile, newi, newj = 11, ix-90, jy
    elif (ix < 0) and (jy < 0): # bottom left
        if tile == 0: #undefined
            pass
        elif tile == 1: #=> 12
            newtile, newi, newj = 12, -jy, 90+ix
        elif tile == 2: #=> 11
            newtile, newi, newj = 11, -jy, 90+ix
        elif tile == 3: #undefined
            pass
        elif tile == 4: #=>0
            newtile, newi, newj = 0, 90+ix, 90+jy
        elif tile == 5: #=>1
            newtile, newi, newj = 1, 90+ix, 90+jy
        elif tile == 6: #undefined
            pass
        elif tile == 7: #=>5 amb
            newtile, newi, newj = 5, 90+ix, 90+jy
        elif tile == 8: #=>5
            newtile, newi, newj = 5, 90+jy, -ix
        elif tile == 9: #=>4
            newtile, newi, newj = 4, 90+jy, -ix
        elif tile == 10: # ambiguous
            newtile, newi, newj = 6, 90+ix, 90+jy
        elif tile == 11: #=>7
            newtile, newi, newj = 7, 90+ix, 90+jy
        elif tile == 12: #=>8
            newtile, newi, newj = 8, 90+ix, 90+jy
    elif (0 <= ix < 90) and (jy < 0): # bottom
        if tile == 0: # undefined
            pass
        elif tile == 1: #=> 0
            newtile, newi, newj = 0, ix, 90+jy
        elif tile == 2: #=> 1
            newtile, newi, newj = 1, ix, 90+jy
        elif tile == 3: # undefined
            pass
        elif tile == 4: #=> 3
            newtile, newi, newj = 3, ix, 90+jy
        elif tile == 5: #=> 4
            newtile, newi, newj = 4, ix, 90+jy
        elif tile == 6: #=> 5
            newtile, newi, newj = 5, ix, 90+jy
        elif tile == 7: #=> 5
            newtile, newi, newj = 5, jy+90, 90-ix
        elif tile == 8: #=> 4
            newtile, newi, newj = 4, jy+90, 90-ix
        elif tile == 9: #=> 3
            newtile, newi, newj = 3, jy+90, 90-ix
        elif tile == 10: #=> 7
            newtile, newi, newj = 7, ix, 90+jy
        elif tile == 11: #=> 8
            newtile, newi, newj = 8, ix, 90+jy
        elif tile == 12: #=> 9
            newtile, newi, newj = 9, ix, 90+jy
    elif (ix >= 90) and (jy < 0): # bottom right
        if tile == 0: #undefined
            pass
        elif tile == 1: #=>3
            newtile, newi, newj = 3, ix-90, 90+jy
        elif tile == 2: #=>4
            newtile, newi, newj = 4, ix-90, 90+jy
        elif tile == 3: #undefined
            pass
        elif tile == 4: #=>9
            newtile, newi, newj = 9, -jy, ix-90
        elif tile == 5: #=>8
            newtile, newi, newj = 8, -jy, ix-90
        elif tile == 6: #undefined
            pass
        elif tile == 7: #=>4
            newtile, newi, newj = 4, jy+90, 180-ix
        elif tile == 8: #=>3
            newtile, newi, newj = 3, jy+90, 180-ix
        elif tile == 9: #undefined
            pass
        elif tile == 10: #=>8
            newtile, newi, newj = 8, ix-90, 90+jy
        elif tile == 11: #=>9
            newtile, newi, newj = 9, ix-90, 90+jy
        elif tile == 12: #undefined
            pass
    else: # within the tile
        pass
    return newtile, newi, newj


def notMoving(uvel, vvel):
    return uvel == 0 and vvel == 0


def get_hfacw(ecco_ds, tile, xi, yj, k):
    # hfacw = ecco_ds.hFacW.values[k,tile,int(yoge),int(xoge)]
    hfacw = float(ecco_ds.hFacW.isel(k=int(k), tile=tile, i_g=int(xi), j=int(yj)).values)
    return hfacw

# Use hFacW to check whether one position is beached.
# Call beached(ecco_ds, particle) or beached(ecco_ds, tile, xoge, yoge, k)
def beached(ecco_ds, tile_or_particle, xoge=None, yoge=None, k=None):
    if type(tile_or_particle) == int:
        tile = tile_or_particle
    else: # particle
        tile = tile_or_particle['tile']
        xoge = tile_or_particle['xoge']
        yoge = tile_or_particle['yoge']
        k = tile_or_particle['k']

    if outOfTile(xoge, yoge):
        return False

    hfacw = get_hfacw(ecco_ds, tile, xoge, yoge, k)
    return int(hfacw) == 0


# Disturb with exponential algorithm, i.e. small flucuation may result big disturbance
def disturb_exp(uvel, vvel, ix, jy, fudge):
    xfactor = random.uniform(-1.0, 1.0) * float(fudge) / 100.0
    yfactor = random.uniform(-1.0, 1.0) * float(fudge) / 100.0

    vel = math.hypot(uvel, vvel)
    xfactor *= 10 ** vel
    yfactor *= 10 ** vel

    if 65 <= ix <= 70:
        # xfactor = -abs(xfactor)
        yfactor = -0.99

    dx = xfactor * uvel * MPS_TO_DEG_PER_MONTH
    dy = yfactor * vvel * MPS_TO_DEG_PER_MONTH
    return dx, dy

# Simple disturbance propotional to uvel/vvel
def disturb_simple(uvel, vvel, kvel, fudge):
    xfactor = random.uniform(-1.0, 1.0) * float(fudge) / 100.0
    yfactor = random.uniform(-1.0, 1.0) * float(fudge) / 100.0
    kfactor = random.uniform(-1.0, 1.0) 
    if fudge < 100.0:
        kfactor *= float(fudge) / 100.0 

    dx = xfactor * uvel * MPS_TO_DEG_PER_MONTH
    dy = yfactor * vvel * MPS_TO_DEG_PER_MONTH
    dk = kfactor * kvel
    return dx, dy, dk

def disturb(uvel, vvel, kvel, ix, jy, k, fudge):
    # return disturb_exp(uvel, vvel, ix, jy, fudge)
    return disturb_simple(uvel, vvel, kvel, fudge)

def nudge(uvel, vvel):
    dx = random.uniform(-0.5, 0.5) if uvel else 0
    dy = random.uniform(-0.5, 0.5) if vvel else 0
    return dx, dy

def mps_to_degreePerMonth():
    return (365.0/12) * 24 * 3600 / (40075017.0/360)

# Move the particle 1 month by its position and vel
# If fudge is set, then add a disturb within (-0.5, 0.5)
# If retry is set, then retry if the particle moves out of tile
def move_1month(ecco_ds, particle, fudge=0, retry=0):
    tile = particle['tile']
    x0,y0,k0 = particle['xoge'], particle['yoge'], particle['k']
    uvel,vvel,kvel = particle['uvel'], particle['vvel'], particle['kvel']

    # new position by ecco data
    x = float(x0) + uvel * MPS_TO_DEG_PER_MONTH
    y = float(y0) + vvel * MPS_TO_DEG_PER_MONTH
    k = k0
    if (k0+kvel) < ecco_ds.hFacW.sizes['k']:
        k = k0 + kvel

    # randomize the move
    dx, dy, dk = disturb(uvel, vvel, kvel, x, y, k, fudge)
    newtile, newx, newy, newk = tile, x+dx, y+dy, k+dk
    if outOfTile(newx, newy):
        newtile, newx, newy = adjustTile(tile, newx, newy)

    # Could retry the move in condition
    for run in range(retry):
        if (not outOfTile(newx, newy)) and beached(ecco_ds, newtile, newx, newy, k):
            logging.debug(f"    {newx}, {newy}, retry {run+1}")
            dx, dy = nudge(uvel, vvel)
            newx, newy = x+dx, y+dy
        if outOfTile(newx, newy):
            newtile, newx, newy = adjustTile(tile, newx, newy)
            logging.debug(f"    {newx}, {newy}, retry {run+1}")
            dx, dy = nudge(uvel, vvel)
            newx, newy = x+dx, y+dy
        else:
            break

    particle['tile'] = newtile
    particle['xoge'], particle['yoge'], particle['k'] = newx, newy, newk


# Read input file in format of tile,x,y,k
# default tile is 10, default k is 0
def read_input(input_file, tile=10, k=0):
    locations = pd.read_csv(input_file)
    particles = []
    for i, row in locations.iterrows():
        particle = {'id': next_particle_id(),
                    'tile': int(row['tile']) if 'tile' in row else tile,
                    'xoge': float(row['xoge']),
                    'yoge': float(row['yoge']),
                    'k': float(row['k']) if 'k' in row else k,
                    'state': 'New'
                   }
        particles.append(particle)
    return particles


# Write results of every particle every month to the disk
def write_results(input_file, results, extra=''):
    # output path is adds results_ to input file name
    input_dir = os.path.dirname(input_file)
    input_fpattern, fext = os.path.splitext(os.path.basename(input_file))
    output_pattern = f'results_{input_fpattern}'
    if input_dir:
        output_pattern = f'{input_dir}/{output_pattern}'
    if extra:
        output_pattern = f'{output_pattern}_{extra}'

    rotate_files(output_pattern)

    output_path = output_pattern + '.csv'
    with open(output_path, mode='w+', newline='') as out_file:
        out_writer = csv.writer(out_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for result in results:
            out_writer.writerow(result)
    logging.info(f' Results are written to {output_path}')
    return output_path


# Update kvel by 1 month
# velocity on k is influences by various factors
def update_kvel(ecco_ds, particle):
    k = particle['k']

    # start from KVEL
    kv = KVEL

    # k represents different depth
    # https://grace.jpl.nasa.gov/data/get-data/ocean-bottom-pressure/
    # There are 46 vertical levels with layer thickness of 10 m in the upper 150 m and 21 layers above 300 m.
    if k < 15:
        pass
    elif 15 <= k < 21:
        kv /= (1 + (k-15) * 0.6)
    elif 21 <= k < 46:
        # assuming avg ocean floor of 3700m, so each k increase 8m
        kv /= (4 + (k-21) * 0.8)
    elif 46 <= k < 50:
        kv /= 150.

    # biological accumulation areas gain density, see MBARI and modis map
    # k layer: 2-11, primarily whenever it passes k=6, 80% gain 20% mass
    lat = float(ecco_ds.YC.isel(
        tile=particle['tile'], j=int(particle['yoge']), i=int(particle['xoge'])).values)
    if 2 <= k <= 11:
        if particle['month'] in [5,6,7,8,9,10]: # summer
            if 45 <= lat <= 75 or -45 <= lat <= -30:
                if random.random() < 0.8:
                    kv *=10
        else: # winter
            if 32 <= lat <= 45 or -60 <= lat <= -40:
                if random.random() < 0.8:
                    kv *=10

    # deep water wave could have more turbulance
    # assuming 1% chance to get large wave
    if 1 <= k <= 10 and random.random() < 0.01:
        kv *= (1.0 + random.uniform(-1.5, 1.5))

    # in south/north pole, uvel/vvel declines very quickly, so reduce the kvel
    # if particle['tile'] in [0,3,6,9,12,4,8,11]:
    #     kv /= 3.0

    particle['kvel'] = kv
    return kv


# uvel/vvel are updated based on ECCO dataset
# kvel is updated based on various factors
def update_velocities(ecco_ds, particle):
    if beached(ecco_ds, particle):
        particle['uvel'] = particle['vvel'] = particle['kvel'] = 0.
        return
    # uvel = ecco_ds.UVEL.values[month,int(k), tile,int(yoge),int(xoge)]
    # vvel = ecco_ds.VVEL.values[month,int(k), tile,int(yoge),int(xoge)]
    particle['uvel'] = float(ecco_ds.UVEL.isel(time=particle['month'],
                                               k=int(particle['k']),
                                               j=int(particle['yoge']),
                                               i_g=int(particle['xoge']),
                                               tile=particle['tile']
                                              ).values)
    particle['vvel'] = float(ecco_ds.VVEL.isel(time=particle['month'],
                                               k=int(particle['k']),
                                               j_g=int(particle['yoge']),
                                               i=int(particle['xoge']),
                                               tile=particle['tile']
                                              ).values)

    update_kvel(ecco_ds, particle)

def find_initial_position(results, id):
    columns = results[0]
    idx = columns.index('id')
    # assuming results is an ordered list, so the first matching index is original position
    for result in results:
        if result[idx] == id:
            return result
    else:
        return None

# Find a random initial particle within the tile
def find_random_in_tile(results, tile):
    tile_initial_list = []
    columns = results[0]
    yearidx = columns.index('year')
    monthidx = columns.index('month')
    year = results[1][yearidx]
    month = results[1][monthidx]

    for result in results:
        if (year == result[yearidx]) and (month == result[monthidx]):
            tile_initial_list.append(result)

    if tile_initial_list:
        pick = int(len(tile_initial_list) * random.random())
        return tile_initial_list[pick]
    else:
        return None


# After the particle is beached, it can choose between
# - no refreshing new particle
# - refresh new particle from initial location
# - refresh a random particle from the tile
def refresh_particle(particle, results):
    # Don't add new particles in later years
    if particle['year'] >= 2005:
        return None

    columns = results[0]
    tileidx = columns.index('tile')
    xogeidx = columns.index('xoge')
    yogeidx = columns.index('yoge')

    if args.no_refresh:
        return None
    elif args.refresh_original:
        # Add new particle from the beached particle's original location
        result = find_initial_position(results, particle['id'])
    elif args.refresh_random:
        # Add new particle randomly from the tile
        result = find_random_in_tile(results, particle['tile'])
    else:
        # By default, use random
        result = find_random_in_tile(results, particle['tile'])

    if result:
        new_particle = {
            'id': next_particle_id(),
            'tile': result[tileidx],
            'xoge': result[xogeidx],
            'yoge': result[yogeidx],
            'k': 0,
            'state': 'New'
        }
        logging.info(f" {particle['year']}/{particle['month']}"
                        f" Refresh {particle['id']} => {new_particle['id']:}"
                        f" {result[tileidx]},{int(result[xogeidx])},{int(result[yogeidx])}")
        return new_particle
    else: # beached from beginning, or something wrong
        return None

def add_to_results(particle, results):
    result = []
    for column in results[0]:
        result.append(particle[column])
    results.append(result)


# Compute particle position of the month, and prepare for the next month
def particle_position(ecco_ds, particle, results, fudge=0):
    # When calling this function, expect the particle to be up-to-date
    # on positions, i.e. index, tile, xoge, yoge, k, and year, month,
    # but vels (uvel, vvel, kvel) need refreshing

    new_particle = None

    # Ignore out-of-tile ones -- only limited tiles are processed
    if outOfTile(particle['xoge'], particle['yoge']):
        logging.info(f"    particle {particle['id']} is out of tile")
        particle['state'] = 'OutOfTile'
        return new_particle

    # If beached then refresh particle.
    elif beached(ecco_ds, particle):
        if particle['state'] != 'Beached':
            new_particle = refresh_particle(particle, results)
        particle['state'] = 'Beached'
    else:
        particle['state'] = 'OK'

    update_velocities(ecco_ds, particle)

    logging.debug(f' {particle}')

    # everything is up-to-date, save it
    add_to_results(particle, results)

    # move to next month's position(tile, x,y,k) based on pos and vel this month
    move_1month(ecco_ds, particle, fudge=fudge, retry=0)

    return new_particle


def hypot(uvel_ds, vvel_ds):
    # uvel_ds and vvel_ds have different coordinates
    # >>> list(uvel_ds.coords)
    # ['i_g', 'k', 'j', 'tile', 'Z', 'dxC', 'rAw', 'dyG', 'PHrefC', 'drF', 'hFacW', 'maskW', 'timestep', 'time']
    # >>> list(vvel_ds.coords)
    # ['j_g', 'k', 'i', 'tile', 'Z', 'rAs', 'dxG', 'dyC', 'PHrefC', 'drF', 'hFacS', 'maskS', 'timestep', 'time']
    # vel_ds = np.hypot(uvel_ds, vvel_ds)
    vel_ds = uvel_ds
    return vel_ds

# Main function to compute the simulation, and write to results_xxx.csv file
def compute(args):
    input_file = args.inputfile
    particles = read_input(input_file)

    global KVEL
    KVEL = float(args.kvel)

    columns = ['id', 'year', 'month', 'tile', 'xoge', 'yoge', 'k', 'uvel', 'vvel', 'kvel', 'state']
    results = [columns]
    base_dir = configure_base_dir()
    for year in range(args.from_year, args.to_year+1):
        ecco_ds = load_ecco_ds(int(year), base_dir)
        for month in range(12):
            for particle in particles:
                particle['year'] = year
                particle['month'] = month
                new_particle = particle_position(ecco_ds, particle, results, fudge=args.fudge_pct)
                if new_particle:
                    particles.append(new_particle)
    extra_info = f'f{args.fudge_pct}kv{args.kvel}'
    result_file = write_results(input_file, results, extra=extra_info)
    return result_file

####################################################
# Plot Section
####################################################

# Decide particle color by its k layer:
# Surface particle (k=0) is blue, k=49 is red
# https://matplotlib.org/stable/tutorials/colors/colors.html
def color_by_k(k):
    # return str(k*5/255.)
    # colors = "bgrcmykw"
    # color = colors[int(k/7)]
    kgrey = int(k * 256 / 50.0)
    # blue --> red during sinkg
    rr = format(kgrey, '02X')
    gg = '00'
    bb = format(255-kgrey, '02X')
    color = f'#{rr}{gg}{bb}'
    return color

# Plot a single tile on year and month, based on results
def plot_tile(ecco_ds, tile, year, month, results, kplot=0):
    uvel_ds = ecco_ds.UVEL.isel(tile=tile, time=month, k=kplot)
    # vvel_ds = ecco_ds.VVEL.isel(tile=tile, time=month, k=kplot)
    # tile_to_plot = hypot(uvel_ds, vvel_ds)
    tile_to_plot = uvel_ds.where(ecco_ds.hFacW.isel(tile=tile,k=kplot) !=0, np.nan)
    plt.imshow(tile_to_plot, cmap='jet', origin='lower', vmin=-0.25, vmax=0.25);
    plt.colorbar()
    plt.xlim([0,90])
    plt.ylim([0,90])
    plt.xlabel('x-dimension of u grid')
    plt.ylabel('y-dimension of v grid')
    plt.title(f'Anna Du Ocean Particle Simulator\nTile {tile} {year}-{(month+1):02}')

    results_match = results[(results.year == year) & (results.month == month) & (results.tile == tile)]
    xoges, yoges, colors = [], [], []
    for index, result in results_match.iterrows():
        logging.debug(f'  {result.tile}{int(result.xoge)},{int(result.yoge)}')
        xoges.append(result.xoge)
        yoges.append(result.yoge)
        colors.append(color_by_k(result.k))
    plt.scatter(xoges, yoges, c=colors)
    plt.tight_layout(pad=0)

# Plot all tiles using the NASA layout
# https://ecco-group.org/images/ecco_tiles.png
def plot_all_tiles(ecco_ds, year, month, results, outfile, k=0):
    logging.info(f'k={k}, tiles=all, {year}-{month}, {outfile}')
    tiles = range(13)
    fig = plt.figure(figsize=(30,30))
    plt.title(f'Anna Du Ocean Particle Simulator\n{year}-{month+1:02}')
    grid5x5 = {0:21, 1:16, 2:11, 3:22, 4:17, 5:12, 6:7, 7:8, 8:9, 9:10, 10:3, 11:4, 12:5}
    for tile in tiles:
        fig = plt.subplot(5,5,grid5x5[tile])
        plot_tile(ecco_ds, tile, year, month, results, kplot=k)
    plt.savefig(outfile)
    # plt.show()
    return outfile

# Plot a single tile of year-month and save to disk
def plot_1tile(ecco_ds, year, month, results, outfile, tile=10, k=0):
    logging.info(f'k={k}, tile={tile}, {year}-{month}, {outfile}')
    fig = plt.figure(figsize = (9,9))
    plot_tile(ecco_ds, tile, year, month, results, kplot=k)
    plt.savefig(outfile)
    return outfile


# convert (tile, x, y) to (longitude, latitude)
def lon_lat(ecco_ds, tile, i, j):
    lon = float(ecco_ds.XC.isel(tile=tile, j=int(j), i=int(i)).values)
    lat = float(ecco_ds.YC.isel(tile=tile, j=int(j), i=int(i)).values)
    return lon, lat


# Plot the global map with longitude-latitude
def plot_all_lonlat(ecco_ds, year, month, results, outfile):
    logging.info(f'{year}-{month}, {outfile}')
    fig = plt.figure(figsize=(30,15))
    uvel_ds = ecco_ds.UVEL.isel(time=month, k=0)
    tile_to_plot = uvel_ds.where(ecco_ds.hFacW.isel(k=0) !=0)
    ecco.plot_proj_to_latlon_grid(ecco_ds.XC, ecco_ds.YC, tile_to_plot,
                plot_type = 'pcolormesh', projection_type = 'PlateCarree',
                cmap='jet', dx=1, dy=1, show_colorbar=True,
                show_grid_labels=True,
                cmin=-0.25, cmax=0.25)
    plt.xlim([-180,180])
    plt.ylim([-90,90])
    plt.xlabel('Longitue')
    plt.ylabel('Latitude')
    plt.title(f"Anna Du Ocean Particle Simulator\n{year}-{month+1:02}")

    results_match = results[(results.year == year) & (results.month == month)]
    lons, lats, colors = [], [], []
    for index, result in results_match.iterrows():
        lon, lat = lon_lat(ecco_ds, result.tile, result.xoge, result.yoge)
        logging.debug(f' {result.tile},{int(result.xoge)},{int(result.yoge)} => ({lon},{lat})')
        lons.append(lon)
        lats.append(lat)
        colors.append(color_by_k(result.k))
    plt.scatter(lons, lats, c=colors)

    plt.savefig(outfile)
    return outfile


# Rotate the previous results csv and mp4 files
def rotate_files(file_pattern):
    file_mp4 = file_pattern+'.mp4'
    file_csv = file_pattern+'.csv'
    if (not os.path.exists(file_mp4)) and (not os.path.exists(file_csv)):
        return

    for i in range(1000):
        mp4_backup = f'{file_pattern}_{i}.mp4'
        csv_backup = f'{file_pattern}_{i}.csv'
        if (not os.path.exists(mp4_backup)) and (not os.path.exists(csv_backup)):
            if os.path.exists(file_mp4):
                logging.info(f'{file_mp4} => {mp4_backup}')
                os.rename(file_mp4, mp4_backup)
            if os.path.exists(file_csv):
                logging.info(f'{file_csv} => {csv_backup}')
                os.rename(file_csv, csv_backup)
            break


# Generate mp4 file with ffmpeg
def gen_mp4(file_pattern, keep_png=False):
    # rotate_files(file_pattern)
    cmd = f'ffmpeg -r 24 -f image2 -s 1920x1080 -i {file_pattern}_%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p -y {file_pattern}.mp4'
    run(cmd, shell=True)
    if not keep_png:
        logging.info(f'rm {file_pattern}_*.png')
        run(f'rm {file_pattern}_*.png', shell=True)
    logging.info(f' Generated {file_pattern}.mp4')


# Main function to plot the result
# Takes results csv from compute() to plot the video
def visualize(args, result_csv, years=[], months=[]):
    base_dir = configure_base_dir()
    results = pd.read_csv(result_csv)
    count = 0
    file_pattern, fext = os.path.splitext(result_csv)
    plot_years = years if years else results.year.unique()
    for year in np.sort(plot_years):
        ecco_ds = load_ecco_ds(int(year), base_dir)
        plot_months = months if months else range(12)
        for month in np.sort(plot_months):
            outfile = f'{file_pattern}_{count:03}.png'
            if args.plot_all_lonlat:
                plot_all_lonlat(ecco_ds, year, month, results, outfile)
            elif args.plot_all_tiles:
                # TODO: use ecco.plot_tiles instead
                plot_all_tiles(ecco_ds, year, month, results, outfile, k=0)
            else:
                plot_1tile(ecco_ds, year, month, results, outfile, tile=args.plot_1tile)
            count += 1
    
    if not (years and months):
        gen_mp4(file_pattern, keep_png=args.keep_png)


# Test function
def test(args):
    base_dir = configure_base_dir()
    result_csv = args.inputfile
    results = pd.read_csv(result_csv)
    count = 0
    k = 0
    tiles = [10, 2]
    fname, fext = os.path.splitext(result_csv)
    file_pattern = f'{fname}-k{k}'
    year, month = 2005, 10
    ecco_ds = load_ecco_ds(int(year), base_dir)
    outfile = f'{file_pattern}_{count:03}.png'
    plot_all_tiles(ecco_ds, year, month, results, outfile, k=0)


# main entry point
def main(args):
    if args.only_plot:
        result_file = args.inputfile
    else:
        result_file = compute(args)

    years, months = [], []
    if args.png_ym:
        y,m = args.png_ym.split(':')
        years = [int(y)]
        months = [int(m)-1]
    if not args.only_compute:
        visualize(args, result_file, years=years, months=months)


if __name__ == '__main__':
    args = usage()
    loglevel = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(level=loglevel)

    if args.test:
        test(args)
    else:
        main(args)
