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

sys.path.append(f'{os.path.expanduser("~")}/ECCOv4-py')
try:
    import ecco_v4_py as ecco
except Exception as e:
    print('Please add ECCOv4-py to System Environment PYTHONPATH!', file=sys.stderr)
    raise(e)

particle_index = 0
def next_particle_index():
    global particle_index
    index = particle_index
    particle_index += 1
    return index

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
    parser.add_argument('--fudge-pct', type=int, default=50, help='Percentage of factor to disturb')
    parser.add_argument('--debug', action='store_true', help='Debug Mode')
    parser.add_argument('--test', action='store_true', help='Test Mode')
    parser.add_argument('--only-plot', action='store_true', help='Only Plot')
    parser.add_argument('--png-ym', help='year:month')
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

def adjustTile(tile, ix, jy):
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
    else: # within the tile
        pass
    return newtile, newi, newj


def notMoving(uvel, vvel):
    return uvel == 0 and vvel == 0


def get_hfacw(ecco_ds, k, tile, xi, yj):
    # hfacw = ecco_ds.hFacW.values[k,tile,int(yoge),int(xoge)]
    hfacw = float(ecco_ds.hFacW.isel(k=int(k), tile=tile, i_g=int(xi), j=int(yj)).values)
    return hfacw

def beached(ecco_ds, xoge, yoge, tile, k):
    hfacw = get_hfacw(ecco_ds, k, tile, xoge, yoge)
    return int(hfacw) == 0


def disturb_exp(uvel, vvel, ix, jy, fudge):
    month_vel_to_pixel = (365.0/12) * 24 * 3600 / (40075017.0/360)
    xfactor = random.uniform(-1.0, 1.0) * float(fudge) / 100.0
    yfactor = random.uniform(-1.0, 1.0) * float(fudge) / 100.0

    vel = math.hypot(uvel, vvel)
    xfactor *= 10 ** vel
    yfactor *= 10 ** vel

    if 65 <= ix <= 70:
        # xfactor = -abs(xfactor)
        yfactor = -0.99

    dx = xfactor * uvel * month_vel_to_pixel
    dy = yfactor * vvel * month_vel_to_pixel
    return dx, dy

def disturb_simple(uvel, vvel, ix, jy, fudge):
    month_vel_to_pixel = (365.0/12) * 24 * 3600 / (40075017.0/360)
    xfactor = random.uniform(-1.0, 1.0) * float(fudge) / 100.0
    yfactor = random.uniform(-1.0, 1.0) * float(fudge) / 100.0

    dx = xfactor * uvel * month_vel_to_pixel
    dy = yfactor * vvel * month_vel_to_pixel
    return dx, dy

def disturb(uvel, vvel, ix, jy, fudge):
    # return disturb_exp(uvel, vvel, ix, jy, fudge)
    return disturb_simple(uvel, vvel, ix, jy, fudge)

def nudge(uvel, vvel):
    dx = random.uniform(-0.5, 0.5) if uvel else 0
    dy = random.uniform(-0.5, 0.5) if vvel else 0
    return dx, dy


# Move the particle 1 month by its position and vel
# If fudge is set, then add a disturb within (-0.5, 0.5)
# If retry is set, then retry if the particle moves out of tile
def move_1month(ecco_ds, tile, x0, y0, k0, uvel, vvel, kvel=1, fudge=0, retry=0):
    month_vel_to_pixel = (365.0/12) * 24 * 3600 / (40075017.0/360)
    x = float(x0) + uvel * month_vel_to_pixel
    y = float(y0) + vvel * month_vel_to_pixel
    k = k0
    if (k0+kvel) < ecco_ds.hFacW.sizes['k']:
        k = k0 + kvel

    dx, dy = disturb(uvel, vvel, x, y, fudge)
    newtile, newx, newy = tile, x+dx, y+dy
    if outOfTile(newx, newy):
        newtile, newx, newy = adjustTile(tile, newx, newy)

    for run in range(retry):
        if (not outOfTile(newx, newy)) and beached(ecco_ds, newx, newy, newtile, k):
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

    return (newtile, newx, newy, k)


def read_input(input_file, tile=10, k=0):
    locations = pd.read_csv(input_file)
    particles = []
    for i, row in locations.iterrows():
        particle = {'index': next_particle_index(),
                    'tile': int(row['tile']) if 'tile' in row else tile,
                    'xoge': float(row['xoge']),
                    'yoge': float(row['yoge']),
                    'k': float(row['k']) if 'k' in row else k
                   }
        particles.append(particle)
    return particles


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


def get_vel(ecco_ds, month, tile, xi, yj, k):
    # uvel = ecco_ds.UVEL.values[month,int(k), tile,int(yoge),int(xoge)]
    # vvel = ecco_ds.VVEL.values[month,int(k), tile,int(yoge),int(xoge)]
    uvel = float(ecco_ds.UVEL.isel(time=month, k=int(k), j=int(yj), i_g=int(xi), tile=tile).values)
    vvel = float(ecco_ds.VVEL.isel(time=month, k=int(k), j_g=int(yj), i=int(xi), tile=tile).values)
    return uvel, vvel

def refresh_particle(particle, results):
    # assuming an ordered list
    for result in results:
        if result[0] == particle['index'] and result[1] == 1992 and result[2] == 0:
            particle['tile'] = result[3]
            particle['xoge'] = result[4]
            particle['yoge'] = result[5]
            particle['k'] = result[6]
            return

# def add_to_results(particle, results):
    # results.append([particle['index'], year, month, tile, xoge, yoge, k, uvel, vvel, kvel])

def particle_position(ecco_ds, particle, year, month, results, fudge=0):
    if outOfTile(particle['xoge'], particle['yoge']):
        return False
    tile = particle['tile']
    k = float(particle['k'])
    xoge = float(particle['xoge'])
    yoge = float(particle['yoge'])

    uvel, vvel = get_vel(ecco_ds, month, tile, xoge, yoge, k)

    kvel = 0.3 # go down 1 k every 10 month
    if beached(ecco_ds, xoge, yoge, tile, k):
        # kvel = 0
        refresh_particle(particle, results)
        tile = particle['tile']
        k = float(particle['k'])
        xoge = float(particle['xoge'])
        yoge = float(particle['yoge'])

    logging.info(f" tile {tile} k {k} PARTICLE {particle['index']} {year}/{month} @ ({xoge},{yoge}) vel: ({uvel}, {vvel}, {kvel})")

    if notMoving(uvel, vvel):
        logging.debug("    Not Moving")

    results.append([particle['index'], year, month, tile, xoge, yoge, k, uvel, vvel, kvel])

    tile, xoge, yoge, k = move_1month(ecco_ds, tile, xoge, yoge, k, uvel, vvel, kvel, fudge=fudge, retry=0)

    particle['tile'] = tile
    particle['xoge'] = xoge
    particle['yoge'] = yoge
    particle['k'] = k
    # only limited tile-to-tile movement is defined
    if outOfTile(particle['xoge'], particle['yoge']):
        logging.debug("    particle will be out of tile")
    elif beached(ecco_ds, xoge, yoge, tile, k):
        logging.debug("    beached!")

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

# https://matplotlib.org/stable/tutorials/colors/colors.html
def color_by_k(k, kvel):
    # return str(k*5/255.)
    colors = "bgrcmykw"
    color = 'w' if kvel == 0 else colors[int(k/7)]
    return color

def plot_tile(ecco_ds, tile, kplot, year, month, results):
    uvel_ds = ecco_ds.UVEL.isel(tile=tile, time=month, k=kplot)
    vvel_ds = ecco_ds.VVEL.isel(tile=tile, time=month, k=kplot)
    tile_to_plot = hypot(uvel_ds, vvel_ds)
    tile_to_plot = tile_to_plot.where(ecco_ds.hFacW.isel(tile=tile,k=kplot) !=0, np.nan)
    plt.imshow(tile_to_plot, cmap='jet', origin='lower', vmin=-0.25, vmax=0.25);
    plt.colorbar()
    plt.xlim([0,90])
    plt.ylim([0,90])
    plt.xlabel('x-dimension of u grid')
    plt.ylabel('y-dimension of v grid')

    plt.title(f'Tile {tile} {year}-{str(month+1).zfill(2)}')
    # results_match = results[(results.year == year) & (results.month == month) & (results.tile == tile) & (results.k == k)]
    results_match = results[(results.year == year) & (results.month == month) & (results.tile == tile)]
    for index, result in results_match.iterrows():
        logging.debug(f'    {int(result.xoge)},{int(result.yoge)}')
        plt.scatter(result.xoge, result.yoge, color=color_by_k(result.k, result.kvel))
    plt.tight_layout(pad=0)

def plot_tiles(ecco_ds, tiles, k, year, month, results, outfile):
    logging.info(f'k={k}, tiles={tiles}, {year}-{month}, {outfile}')
    fig = plt.figure(figsize=(12,12))
    plt.title(f'{year}-{month+1}')
    for tile in tiles:
        if tile == 10:
            fig = plt.subplot(222)
            plot_tile(ecco_ds, tile, k, year, month, results)
        elif tile == 2:
            fig = plt.subplot(223)
            plot_tile(ecco_ds, tile, k, year, month, results)
    plt.savefig(outfile)
    # plt.show()
    return outfile

def plot_tile_10(ecco_ds, k, year, month, results, outfile):
    tile = 10
    logging.info(f'k={k}, tile={tile}, {year}-{month}, {outfile}')
    fig = plt.figure(figsize = (9,9))
    plot_tile(ecco_ds, tile, k, year, month, results)
    plt.savefig(outfile)
    return outfile

def plot_tiles_2_10(ecco_ds, k, year, month, results, outfile):
    fig = plt.figure(figsize=(12,6))

    # tile 2
    tile = 2
    fig = plt.subplot(122)
    uvel_ds = ecco_ds.UVEL.isel(tile=tile, time=month, k=k)
    vvel_ds = ecco_ds.VVEL.isel(tile=tile, time=month, k=k)
    tile_to_plot = hypot(uvel_ds, vvel_ds)
    tile_to_plot = tile_to_plot.where(ecco_ds.hFacW.isel(tile=tile,k=k) !=0, np.nan)
    plt.imshow(tile_to_plot, origin='lower', vmin=-0.25, vmax=0.25);
    plt.colorbar()
    plt.title(f'Tile {tile} {year}-{str(month+1).zfill(2)}')
    results_match = results[(results.year == year) & (results.month == month) & (results.tile == tile) & (results.k == k)]
    for index, result in results_match.iterrows():
        logging.debug(f'    {int(result.xoge)},{int(result.yoge)}')
        plt.scatter(result.xoge, result.yoge, color='magenta')
    plt.tight_layout(pad=0)

    # tile 10
    tile = 10
    fig = plt.subplot(121)
    uvel_ds = ecco_ds.UVEL.isel(tile=tile, time=month, k=k)
    vvel_ds = ecco_ds.VVEL.isel(tile=tile, time=month, k=k)
    tile_to_plot = hypot(uvel_ds, vvel_ds)
    tile_to_plot = tile_to_plot.where(ecco_ds.hFacW.isel(tile=tile,k=k) !=0, np.nan)
    arr_to_plot = tile_to_plot
    # arr_to_plot = np.swapaxes(tile_to_plot.values, 0, 1)
    # arr_to_plot = tile_to_plot.copy()
    # # This has a lot of computation
    # for i in range(90):
    #     for j in range(90):
    #         arr_to_plot[i][j] = tile_to_plot[j][89-i]
    plt.imshow(arr_to_plot, origin='lower', vmin=-0.25, vmax=0.25);
    plt.colorbar()
    plt.title(f'Tile {tile} {year}-{str(month+1).zfill(2)}')
    results_match = results[(results.year == year) & (results.month == month) & (results.tile == tile) & (results.k == k)]
    for index, result in results_match.iterrows():
        logging.debug(f'    {int(result.xoge)},{int(result.yoge)}')
        plt.scatter(result.yoge, 89-result.xoge, color='magenta')
    plt.tight_layout(pad=0)

    plt.savefig(outfile)
    # plt.show()
    return outfile


def rotate_files(file_pattern):
    file_mp4 = file_pattern+'.mp4'
    file_csv = file_pattern+'.csv'
    if (not os.path.exists(file_mp4)) and (not os.path.exists(file_csv)):
        return

    for i in range(100):
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


def gen_mp4(file_pattern, keep_png=False):
    # rotate_files(file_pattern)
    cmd = f'ffmpeg -r 24 -f image2 -s 1920x1080 -i {file_pattern}_%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p -y {file_pattern}.mp4'
    run(cmd, shell=True)
    if not keep_png:
        logging.info(f'rm {file_pattern}_*.png')
        run(f'rm {file_pattern}_*.png', shell=True)
    logging.info(f' Generated {file_pattern}.mp4')


def visualize(result_csv, k, years=[], months=[]):
    base_dir = configure_base_dir()
    results = pd.read_csv(result_csv)
    count = 0
    tiles = [10, 2]
    file_pattern, fext = os.path.splitext(result_csv)
    plot_years = years if years else results.year.unique()
    for year in np.sort(plot_years):
        ecco_ds = load_ecco_ds(int(year), base_dir)
        plot_months = months if months else range(12)
        for month in np.sort(plot_months):
            outfile = f'{file_pattern}_{count:03}.png'
            plot_tile_10(ecco_ds, k, year, month, results, outfile)
            # plot_tiles(ecco_ds, tiles, k, year, month, results, outfile)
            # plot_tiles_2_10(ecco_ds, k, year, month, results, outfile)
            count += 1
    
    if not (years and months):
        gen_mp4(file_pattern, keep_png=args.debug)


def compute(args):
    input_file = args.inputfile
    particles = read_input(input_file, args.tile, args.k)

    results = [['particle', 'year', 'month', 'tile', 'xoge', 'yoge', 'k', 'uvel', 'vvel', 'kvel']]
    base_dir = configure_base_dir()
    for year in range(args.from_year, args.to_year):
        ecco_ds = load_ecco_ds(int(year), base_dir)
        for month in range(12):
            for particle in particles:
                particle['year'] = year
                particle['month'] = month
                particle_position(ecco_ds, particle, year, month, results, fudge=args.fudge_pct)
    extra_info = f'f{args.fudge_pct}'
    result_file = write_results(input_file, results, extra=extra_info)
    return result_file

def test(args):
    base_dir = configure_base_dir()
    result_csv = args.inputfile
    results = pd.read_csv(result_csv)
    count = 0
    k = args.k
    tiles = [10, 2]
    fname, fext = os.path.splitext(result_csv)
    file_pattern = f'{fname}-k{k}'
    year, month = 2005, 10
    ecco_ds = load_ecco_ds(int(year), base_dir)
    outfile = f'{file_pattern}_{count:03}.png'
    plot_tiles(ecco_ds, tiles, k, year, month, results, outfile)


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
    visualize(result_file, args.k, years=years, months=months)


if __name__ == '__main__':
    args = usage()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    if args.test:
        test(args)
    else:
        main(args)
