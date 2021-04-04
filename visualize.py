import numpy as np
import xarray as xr
import sys
import matplotlib.pyplot as plt
import json
import warnings
warnings.filterwarnings('ignore')
sys.path.append('/home/ec2-user/ECCOv4-py')
from cartopy import config
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import ecco_v4_py as ecco
import matplotlib.path as mpath
import csv
import random
import logging
import argparse

logging.basicConfig(level=logging.DEBUG)


# base_dir = '/eccodata'
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
    ecco_ds = load_ecco_ds(1992, '/eccodata')
    plot_vel(ecco_ds, 10, 0, 0)

if __name__ == '__main__':
    main()