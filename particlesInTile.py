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
import random
import logging
import argparse
from subprocess import run

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
    parser.add_argument('--only-plot', action='store_true', help='Only Plot')
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

def tileTo(tile, ix, jy):
    newtile, newi, newj = tile, ix, jy
    if (ix < 0) and (jy >= 90): # top left
        pass
    elif (0 <= ix < 90) and (jy >= 90): # top
        if tile == 2:
            newtile, newi, newj = (6, jy-90, 90-ix)
        elif tile == 10:
            newtile, newi, newj = 2, jy-90, 90-ix
    elif (ix >= 90) and (jy >= 90): # top right
        if tile == 2:
            pass # undefineable
        elif tile == 10:
            newtile, newi, newj = 1, jy-90, 180-ix
    elif (ix < 0) and (0 <= jy < 90): # left
        if tile == 2:
            newtile, newi, newj = 10, 90-jy, 90+ix
        elif tile == 10:
            newtile, newi, newj = 6, 90-jy, 90+ix
    elif (ix >= 90) and (0 <= jy < 90): # right
        if tile == 2:
            newtile, newi, newj = 5, ix-90, jy
        elif tile == 10:
            newtile, newi, newj = 11, ix-90, jy
    elif (ix < 0) and (jy < 0): # bottom left
        if tile == 2:
            newtile, newi, newj = 11, -jy, 90+ix
        elif tile == 10:
            newtile, newi, newj = 6, 90+ix, 90+jy # also ambuguous
    elif (0 <= ix < 90) and (jy < 0): # bottom
        if tile == 2:
            newtile, newi, newj = 1, ix, 90+jy
        elif tile == 10:
            newtile, newi, newj = 7, ix, 90+jy
    elif (ix >= 90) and (jy < 0): # bottom right
        if tile == 2:
            newtile, newi, newj = 4, ix-90, 90+jy
        elif tile == 10:
            newtile, newi, newj = 8, ix-90, 90+jy
    return newtile, newi, newj


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


def read_input(input_file, tile, k):
    with open(input_file) as csv_file:
        csv_reader = csv.DictReader(csv_file, delimiter=',')
        particles = []
        index = 0
        for row in csv_reader:
            particle = {'index': index, 'xoge': int(row['xoge']), 'yoge': int(row['yoge']), 'state': 'ok'}
            
            particle['tile'] = int(row['tile']) if 'tile' in row else tile
            particle['k'] = int(row['k']) if 'k' in row else k

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
    return output_path


def get_vel(ecco_ds, k, month, tile, xi, yj):
    # uvel = ecco_ds.UVEL.values[month,k, tile,int(yoge),int(xoge)]
    # vvel = ecco_ds.VVEL.values[month,k, tile,int(yoge),int(xoge)]
    uvel = float(ecco_ds.UVEL.isel(time=month, k=k, j=int(yj), i_g=int(xi), tile=tile).values)
    vvel = float(ecco_ds.VVEL.isel(time=month, k=k, j_g=int(yj), i=int(xi), tile=tile).values)
    return uvel, vvel


def particle_position(ecco_ds, particle, year, month, results, disturbing=False):
    if particle['state'] == 'OutOfTile':
        return False
    tile = particle['tile']
    k = particle['k']
    xoge = float(particle['xoge'])
    yoge = float(particle['yoge'])

    uvel, vvel = get_vel(ecco_ds, k, month, tile, xoge, yoge)

    logging.info(f" tile {tile} k {k} PARTICLE {particle['index']} {year}/{month} @ ({xoge},{yoge}) vel: ({uvel}, {vvel})")

    if notMoving(uvel, vvel):
        logging.debug("    Not Moving")

    if beached(ecco_ds, xoge, yoge, tile, k):
        particle['state'] = 'Beached'
        logging.debug("    beached!")

    results.append([particle['index'], year, month, tile, xoge, yoge, k, uvel, vvel])

    xoge, yoge = particle['xoge'], particle['yoge'] = move_1month(ecco_ds, xoge, yoge, uvel, vvel, tile, k, disturbing=disturbing, retry=4)
    if outOfTile(xoge, yoge):
        particle['tile'], particle['xoge'], particle['yoge'] = tileTo(tile, xoge, yoge)
    # only limited tile-to-tile movement is defined
    if outOfTile(particle['xoge'], particle['yoge']):
        particle['state'] = 'OutOfTile'
        logging.debug("    particle will be out of tile")
    return True


def hypot(uvel_ds, vvel_ds):
    # uvel_ds and vvel_ds have different coordinates
    # >>> list(uvel_ds.coords)
    # ['i_g', 'k', 'j', 'tile', 'Z', 'dxC', 'rAw', 'dyG', 'PHrefC', 'drF', 'hFacW', 'maskW', 'timestep', 'time']
    # >>> list(vvel_ds.coords)
    # ['j_g', 'k', 'i', 'tile', 'Z', 'rAs', 'dxG', 'dyC', 'PHrefC', 'drF', 'hFacS', 'maskS', 'timestep', 'time']
    # vel_ds = np.hypot(uvel_ds, vvel_ds)
    vel_ds = uvel_ds
    return vel_ds


def plot_vel(ecco_ds, tile, k, year, month, results, outfile):
    fig = plt.figure(figsize=(9,9))
    uvel_ds = ecco_ds.UVEL.isel(tile=tile, time=month, k=k)
    vvel_ds = ecco_ds.VVEL.isel(tile=tile, time=month, k=k)
    logging.info(f'k={k}, tile={tile}, {year}-{month}')
    tile_to_plot = hypot(uvel_ds, vvel_ds)
    tile_to_plot = tile_to_plot.where(ecco_ds.hFacW.isel(tile=tile,k=k) !=0, np.nan)
    plt.imshow(tile_to_plot, origin='lower', vmin=-0.25, vmax=0.25);
    plt.colorbar()
    plt.title(f'{year}-{month+1}')
    results_month = results[(results.year == year) & (results.month == month)]
    for index, result in results_month.iterrows():
        logging.debug(f'    {int(result.xoge)},{int(result.yoge)}')
        plt.scatter(result.xoge, result.yoge, color='black')
    plt.savefig(outfile)
    # plt.show()


def rotate_file(file_pattern):
    file_to_save = file_pattern+'.mp4'
    if os.path.isfile(file_to_save):
        for i in range(100):
            file_backup = f'{file_pattern}_{i}.mp4'
            if not os.path.exists(file_backup):
                os.rename(file_to_save, file_backup)
                break


def gen_mp4(file_pattern, keep_png=False):
    rotate_file(file_pattern)
    cmd = f'ffmpeg -r 24 -f image2 -s 1920x1080 -i {file_pattern}_%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p -y {file_pattern}.mp4'
    run(cmd, shell=True)
    if not keep_png:
        run(f'rm {file_pattern}_*.png', shell=True)
    logging.info(f' Generated {file_pattern}.mp4')


def visualize(result_csv, k):
    base_dir = configure_base_dir()
    results = pd.read_csv(result_csv)
    count = 0
    tile = 10
    fname, fext = os.path.splitext(result_csv)
    file_pattern = f'{fname}-k{k}'
    for year in np.sort(results.year.unique()):
        ecco_ds = load_ecco_ds(int(year), base_dir)
        for month in range(12):
            outfile = f'{file_pattern}_{count:03}.png'
            plot_vel(ecco_ds, tile, k, year, month, results, outfile)
            count += 1
    
    gen_mp4(file_pattern, keep_png=args.debug)


def compute(args):
    input_file = args.inputfile
    particles = read_input(input_file, args.tile, args.k)

    results = [['particle', 'year', 'month', 'tile', 'xoge', 'yoge', 'k', 'uvel', 'vvel']]
    base_dir = configure_base_dir()
    for year in range(args.from_year, args.to_year):
        ecco_ds = load_ecco_ds(int(year), base_dir)
        for month in range(12):
            for particle in particles:
                particle_position(ecco_ds, particle, year, month, results, disturbing=args.disturbing)
    result_file = write_results(input_file, results)
    return result_file

def main(args):
    if args.only_plot:
        result_file = args.inputfile
    else:
        result_file = compute(args)
    visualize(result_file, args.k)


if __name__ == '__main__':
    args = usage()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    main(args)
