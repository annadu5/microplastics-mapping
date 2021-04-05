import numpy as np
import xarray as xr
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
import random
import logging
import argparse

sys.path.append(f'{os.path.expanduser("~")}/ECCOv4-py')
try:
    import ecco_v4_py as ecco
except Exception as e:
    print('Please add ECCOv4-py to System Environment PYTHONPATH!', file=sys.stderr)
    raise(e)


logging.basicConfig(level=logging.DEBUG)

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


def load_ecco_ds(year, base_dir, vars=['UVEL', 'VVEL']):
    # ECCO_dir = base_dir + '/Version4/Release3_alt'
    ECCO_dir = base_dir + '/Version4/Release4'
    grid_dir = f'{ECCO_dir}/nctiles_grid'
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
    parser.add_argument('-k', '--k', type=int, default=0, help='k layer')
    parser.add_argument('--tile', type=int, default=10, help='tile number [0,12]')
    # TODO add --from-year --to-year
    args = parser.parse_args()
    return args


def outOfTile(x,y):
    if y < 0:
        # anything else?
        return True
    elif y > 90:
        return True
    elif x < 0:
        return True
    elif x > 90:
        return True
    else:
        return False


def notMoving(uvel, vvel):
    return uvel == 0 and vvel == 0


def beached(ecco_ds, xoge, yoge, tile, k):
    hfacw = ecco_ds.hFacW.values[k,tile,int(yoge),int(xoge)]
    logging.debug(hfacw)
    return int(hfacw) == 0


def disturb(x, y):
    return (random.uniform(-0.5, 0.5), random.uniform(-0.5, 0.5))


# Move the particle 1 month by its position and vel
# If fudge is set, then add a disturb within (-0.5, 0.5)
# If retry is set, then retry if the particle moves out of tile
def move_1month(ecco_ds, x0, y0, uvel, vvel, tile, k, fudge=False, retry=0):
    month_vel_to_pixel = (365.0/12) * 24 * 3600 / (40075017.0/360)
    x = float(x0) + uvel * month_vel_to_pixel
    y = float(y0) + vvel * month_vel_to_pixel

    if fudge:
        # Fudge around half a pixel
        dx, dy = disturb(x, y)
    else:
        dx = dy = 0

    for run in range(retry):
        if outOfTile(x+dy, y+dy) or beached(ecco_ds, x+dx, y+dy, tile, k):
            logging.debug(f"    r{x+dx}, {y+dy}, etry {run+1}")
            dx, dy = disturb(x, y)
        else:
            break

    return (x+dx, y+dy)


counter = 0
def particle_positions(particle, xoge0, yoge0, year_range, tile=10, k=0):
    monthly = []
    xoge = xoge0
    yoge = yoge0
    global counter
    base_dir = configure_base_dir()
    for year in year_range:  # good stretch is 1999 to 2009  BELOW UNDER K and TILES NEED TO MAKE VARS
        ecco_ds = load_ecco_ds(int(year), base_dir)

        for month in range (12):   # 12 vs 11
            counter += 1
            
            logging.debug(f" tile {tile} k {k} PARTICLE {particle} {year}/{month} @ ({xoge},{yoge})")

            uvel = ecco_ds.UVEL.values[month,tile,k,int(yoge),int(xoge)] # Here the first threeo of these are correct
            vvel = ecco_ds.VVEL.values[month,tile,k,int(yoge),int(xoge)] # m/s needs to be converted into a distance -- this is a VELOCITY

            logging.debug(f"    (uvel,vvel)=({uvel},{vvel})    ecco_ds.hFacW.dims: {ecco_ds.hFacW.dims}")

            if notMoving(uvel, vvel):
                logging.debug("    NOT MOVING")

            if beached(ecco_ds, xoge, yoge, tile, k):
                logging.debug("    beached!")


            monthly.append([tile, k, particle, year, month, xoge, yoge, uvel, vvel])

            xoge, yoge = move_1month(ecco_ds, xoge, yoge, uvel, vvel, tile, k, fudge=True, retry=4)
            if outOfTile(xoge, yoge):
                logging.debug("    particle will be out of tile")
                return monthly
    return monthly
                    
# Above counter increments the particle count                


args = usage()
input_file = args.inputfile
k = args.k
tile = args.tile

results = [['tile', 'k', 'particle', 'year', 'month', 'xoge', 'yoge', 'uvel', 'vvel']]
with open(input_file) as csv_file:
    csv_reader = csv.DictReader(csv_file, delimiter=',')
    particle = 0
    for row in csv_reader:
        xoge = int(row['xoge'])
        yoge = int(row['yoge'])

        # With recursive_load_ecco_var_from_years_nc one can specify one or more variables to load, one or more years to load, while also requesting only a subset of tiles and vertical levels.
      # ecco.recursive_load_ecco_var_from_years_nc(d
        particle_result = particle_positions(particle, xoge, yoge, range(1992,2015), tile=args.tile, k=args.k)
        results.extend(particle_result)

        particle += 1

# output path is adds results_ to input file
output_path = 'results_' + os.path.basename(input_file)
if os.path.dirname(input_file):
    output_path = os.path.dirname(input_file) + output_path
with open(output_path, mode='w+', newline='') as out_file:
    out_writer = csv.writer(out_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for result in results:
        out_writer.writerow(result)
