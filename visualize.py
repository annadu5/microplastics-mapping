import numpy as np
import xarray as xr
import pandas as pd
import logging
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
import argparse
from subprocess import run

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

# base_dir of `eccodata`
def load_ecco_ds(year, base_dir, vars=['UVEL', 'VVEL']):
    # ECCO_dir = base_dir + '/Version4/Release3_alt'
    ECCO_dir = base_dir + '/Version4/Release4'
    grid_dir = f'{ECCO_dir}/nctiles_grid'
    # ecco_grid = ecco.load_ecco_grid_nc(grid_dir, 'ECCOv4r3_grid.nc')
    ecco_grid = ecco.load_ecco_grid_nc(grid_dir, 'ECCO-GRID.nc')
    day_mean_dir = f'{ECCO_dir}/nctiles_monthly'
    ecco_vars = ecco.recursive_load_ecco_var_from_years_nc(
        day_mean_dir, vars_to_load=vars, years_to_load=year #, dask_chunk=False
        )
    ecco_ds = xr.merge((ecco_grid , ecco_vars))
    return ecco_ds


def plot_vel(ecco_ds, tile, k, year, month, results, outfile):
    fig = plt.figure(figsize=(9,9))
    lons = np.copy(ecco_ds.XC.sel(tile=tile))
    # we must convert the longitude coordinates from
    # [-180 to 180] to [0 to 360]
    # because of the crossing of the international date line.
    lons[lons < 0] = lons[lons < 0]+360
    lats = ecco_ds.YC.sel(tile=tile)
    tile_to_plot = ecco_ds.UVEL.isel(tile=tile, time=month, k=0)
    tile_to_plot= tile_to_plot.where(ecco_ds.hFacW.isel(tile=tile,k=0) !=0, np.nan)
    plt.pcolor(lons, lats, tile_to_plot, vmin=-.25, vmax=.25, cmap='jet')
    # plt.pcolor(lons, lats, tile_to_plot)
    # plt.imshow(tile_to_plot, origin='lower');
    plt.colorbar()
    results_month = results[(results.year == year) & (results.month == month)]
    print(type(results_month))
    for index, result in results_month.iterrows():
        plt.scatter(result.xoge, result.yoge, color='black')
    plt.savefig(outfile)
    # plt.show()


def gen_mp4(file_pattern):
    cmd = f'ffmpeg -r 30 -f image2 -s 1920x1080 -i {file_pattern}_%03d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p {file_pattern}.mp4'
    run(cmd, shell=True)


def main():
    base_dir = configure_base_dir()
    # ecco_ds = load_ecco_ds(2005, base_dir)
    # plot_vel(ecco_ds, 10, 0, 0)
    particles_results_file = 'results_klawinput.csv'
    results = pd.read_csv(particles_results_file)
    count = 0
    tile = 10
    k = 0
    file_pattern = f'results_{tile}_{k}'
    for year in np.sort(results.year.unique()):
        ecco_ds = load_ecco_ds(int(year), base_dir)
        for month in range(12):
            outfile = f'{file_pattern}_{count:03}.png'
            plot_vel(ecco_ds, tile, k, year, month, results, outfile)
            count += 1
    
    gen_mp4(file_pattern)

if __name__ == '__main__':
    main()