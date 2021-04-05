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
    parser.add_argument('-k', '--k', type=int, default=0, help='k layer')
    parser.add_argument('--tile', type=int, default=10, help='tile number [0,12]')
    parser.add_argument('--from-year', type=int, default=1992, help='Starting year')
    parser.add_argument('--to-year', type=int, default=2015, help='End year')
    parser.add_argument('--disturbing', action='store_true', help='Add turbulence to particle movement')
    parser.add_argument('--debug', action='store_true', help='Debug Mode')
    args = parser.parse_args()
    return args


def outOfTile(x,y):
    if y < 0:
        # anything else?
        return True
    elif y >= 90:
        return True
    elif x < 0:
        return True
    elif x >= 90:
        return True
    else:
        return False


def notMoving(uvel, vvel):
    return uvel == 0 and vvel == 0


def get_hfacw(ecco_ds, k, tile, xi, yj):
    # hfacw = ecco_ds.hFacW.values[k,tile,int(yoge),int(xoge)]
    hfacw = float(ecco_ds.hFacW.isel(k=k, tile=tile, i_g=int(xi), j=int(yj)).values)
    return hfacw

def beached(ecco_ds, xoge, yoge, tile, k):
    hfacw = get_hfacw(ecco_ds, k, tile, xoge, yoge)
    return int(hfacw) == 0


def disturb(uvel, vvel):
    dx = random.uniform(-0.5, 0.5) if uvel else 0
    dy = random.uniform(-0.5, 0.5) if vvel else 0
    return dx, dy


# Move the particle 1 month by its position and vel
# If disturbing is set, then add a disturb within (-0.5, 0.5)
# If retry is set, then retry if the particle moves out of tile
def move_1month(ecco_ds, x0, y0, uvel, vvel, tile, k, disturbing=False, retry=0):
    month_vel_to_pixel = (365.0/12) * 24 * 3600 / (40075017.0/360)
    x = float(x0) + uvel * month_vel_to_pixel
    y = float(y0) + vvel * month_vel_to_pixel

    dx = dy = 0
    if disturbing:
        # Fudge around half a pixel
        dx, dy = disturb(uvel, vvel)

    for run in range(retry):
        if outOfTile(x+dx, y+dy):
            logging.debug(f"    {x+dx}, {y+dy}, retry {run+1}")
            dx, dy = disturb(uvel, vvel)
        elif beached(ecco_ds, x+dx, y+dy, tile, k):
            logging.debug(f"    {x+dx}, {y+dy}, retry {run+1}")
            dx, dy = disturb(uvel, vvel)
        else:
            break

    return (x+dx, y+dy)


def read_input(input_file):
    with open(input_file) as csv_file:
        csv_reader = csv.DictReader(csv_file, delimiter=',')
        particles = []
        index = 0
        for row in csv_reader:
            particle = {'index': index, 'xoge': int(row['xoge']), 'yoge': int(row['yoge']), 'state': 'ok'}
            particles.append(particle)
            index += 1
    return particles


def write_results(input_file, results):
    # output path is adds results_ to input file
    output_path = 'results_' + os.path.basename(input_file)
    if os.path.dirname(input_file):
        output_path = os.path.dirname(input_file) + output_path
    with open(output_path, mode='w+', newline='') as out_file:
        out_writer = csv.writer(out_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for result in results:
            out_writer.writerow(result)
    logging.info(f' Results are written to {output_path}')


def get_vel(ecco_ds, k, month, tile, xi, yj):
    # uvel = ecco_ds.UVEL.values[month,k, tile,int(yoge),int(xoge)]
    # vvel = ecco_ds.VVEL.values[month,k, tile,int(yoge),int(xoge)]
    uvel = float(ecco_ds.UVEL.isel(time=month, k=k, j=int(yj), i_g=int(xi), tile=tile).values)
    vvel = float(ecco_ds.VVEL.isel(time=month, k=k, j_g=int(yj), i=int(xi), tile=tile).values)
    return uvel, vvel


def particle_position(ecco_ds, particle, tile, k, year, month, results, disturbing=False):
    if particle['state'] == 'OutOfTile':
        return False
    xoge = float(particle['xoge'])
    yoge = float(particle['yoge'])

    uvel, vvel = get_vel(ecco_ds, k, month, tile, xoge, yoge)

    logging.info(f" tile {tile} k {k} PARTICLE {particle['index']} {year}/{month} @ ({xoge},{yoge}) vel: ({uvel}, {vvel})")

    if notMoving(uvel, vvel):
        logging.debug("    Not Moving")

    if beached(ecco_ds, xoge, yoge, tile, k):
        particle['state'] = 'Beached'
        logging.debug("    beached!")

    results.append([tile, k, particle['index'], year, month, xoge, yoge, uvel, vvel])

    xoge, yoge = particle['xoge'], particle['yoge'] = move_1month(ecco_ds, xoge, yoge, uvel, vvel, tile, k, disturbing=disturbing, retry=4)
    if outOfTile(xoge, yoge):
        particle['state'] = 'OutOfTile'
        logging.debug("    particle will be out of tile")
    return True

def main(args):
    tile = args.tile
    k = args.k
    input_file = args.inputfile
    particles = read_input(input_file)

    results = [['tile', 'k', 'particle', 'year', 'month', 'xoge', 'yoge', 'uvel', 'vvel']]
    base_dir = configure_base_dir()
    for year in range(args.from_year, args.to_year):
        ecco_ds = load_ecco_ds(int(year), base_dir)
        for month in range(12):
            for particle in particles:
                particle_position(ecco_ds, particle, tile, k, year, month, results, disturbing=args.disturbing)
    write_results(input_file, results)


if __name__ == '__main__':
    args = usage()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    main(args)
