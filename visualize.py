import numpy as np
import xarray as xr
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
def load_ecco_ds(year, base_dir):
    grid_dir = f'{base_dir}/Version4/Release3_alt/nctiles_grid'
    ecco_grid = ecco.load_ecco_grid_nc(grid_dir, 'ECCOv4r3_grid.nc')
    day_mean_dir = f'{base_dir}/Version4/Release3_alt/nctiles_monthly'
    ecco_vars = ecco.recursive_load_ecco_var_from_years_nc(
        day_mean_dir, vars_to_load=['UVEL', 'VVEL'], years_to_load=year, dask_chunk=False)
    ecco_ds = xr.merge((ecco_grid , ecco_vars))
    return ecco_ds


def plot_vel(ecco_ds, tile, k, month):
    fig = plt.figure(figsize=(9,9))
    lons = np.copy(ecco_ds.XC.sel(tile=tile))
    # we must convert the longitude coordinates from
    # [-180 to 180] to [0 to 360]
    # because of the crossing of the international date line.
    lons[lons < 0] = lons[lons < 0]+360
    lats = ecco_ds.YC.sel(tile=tile)
    tile_to_plot = ecco_ds.UVEL.isel(tile=tile, time=month)
    tile_to_plot= tile_to_plot.where(ecco_ds.hFacW.isel(tile=tile,k=month) !=0, np.nan)
    plt.pcolor(lons, lats, tile_to_plot, vmin=-.25, vmax=.25, cmap='jet')
    plt.colorbar()
    plt.savefig(f'{tile}.png')

def main():
    base_dir = configure_base_dir()
    ecco_ds = load_ecco_ds(1992, base_dir)
    plot_vel(ecco_ds, 10, 0, 0)

if __name__ == '__main__':
    main()