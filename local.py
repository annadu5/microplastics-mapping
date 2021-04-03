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

logging.basicConfig(level=logging.DEBUG)

base_dir = '/eccodata'
grid_dir = f'{base_dir}/Version4/Release3_alt/nctiles_grid'
ecco_grid = ecco.load_ecco_grid_nc(grid_dir, 'ECCOv4r3_grid.nc')
day_mean_dir = f'{base_dir}/Version4/Release3_alt/nctiles_monthly'


def move_1month(x0, y0, uvel, vvel, fudge=False):
    month_vel_to_pixel = (365.0/12) * 24 * 3600 / (40075017.0/360)
    dotx = float(x0) + uvel * month_vel_to_pixel
    doty = float(y0) + vvel * month_vel_to_pixel
    if fudge:
        # Fudge around half a pixel
        dotx += random.uniform(-0.5, 0.5)
        doty += random.uniform(-0.5, 0.5)
    return (dotx, doty)


def outOfTile(x,y):
    if y < 0:
        logging.debug("stuck on bottom")
        return True
    elif y > 90:
        logging.debug("stuck on top")
        return True
    elif x < 0:
        logging.debug("stuck on left")
        return True
    elif x > 90:
        logging.debug("stuck on right")
        return True
    else:
        return False


counter = 0
def particle_positions(particle, xoge0, yoge0, year_range, tile=10, k=0):
    monthly = []
    xoge = xoge0
    yoge = yoge0
    global counter
    for year in year_range:  # good stretch is 1999 to 2009  BELOW UNDER K and TILES NEED TO MAKE VARS
        ecco_vars = ecco.recursive_load_ecco_var_from_years_nc(
                day_mean_dir, vars_to_load=['UVEL', 'VVEL'], years_to_load=year, dask_chunk=False)
                #k_subset = [0], \
                #tiles_to_load = [10], \

        for month in range (12):   # 12 vs 11
            counter += 1
            
            logging.debug(" ")
            logging.debug(" ")
            logging.debug(" --------")
            logging.debug(f" PARTICLE {particle} STARTING COORD ({xoge},{yoge}) {year}/{month+1}")

            
            ecco_ds = xr.merge((ecco_grid , ecco_vars))
            ecco_ds.attrs = []

            uvel = ecco_ds.UVEL.values[month,tile,k,int(yoge),int(xoge)] # Here the first threeo of these are correct
            vvel = ecco_ds.VVEL.values[month,tile,k,int(yoge),int(xoge)] # m/s needs to be converted into a distance -- this is a VELOCITY

            logging.debug(f' (x,y)=({xoge},{yoge}), (uvel,vvel)=({uvel},{vvel})')

            logging.debug(f" ecco_ds.hFacW.dims: {ecco_ds.hFacW.dims}")
            

            if abs(uvel) == 0 and abs(vvel) == 0:
                    logging.debug("NOT MOVING")
                                        
            else:
                # THIS ONE IS NOW RIGHT# was backwards  [k,tile,j,i_g] (time, tile, k j i 
                beachhfacw = ecco_ds.hFacW.values[k,tile,int(yoge),int(xoge)]
                logging.debug(beachhfacw)
                if (int(beachhfacw) != 0):
                    plt.scatter(xoge,yoge,color='magenta')
                    logging.debug("notbeached")
                else:
                    logging.debug("beached!")


            monthly.append([tile, k, particle, year, month, xoge, yoge, uvel, vvel])
            # below need to have it open a new file automatically
            #newcsvfile =  'particles_out\\' + str("{:02d}".format(particle)) + 'eccodatasetoutput.csv'
            #f= open(newcsvfile,"rw")
            #f.close() 
            #with open('particles_out/' + str("{:02d}".format(particle)) + 'eccodatasetoutput.csv', mode='a+', newline='') as out_file:
                #out_writer = csv.writer(out_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                #out_writer.writerow([str(counter),str(particle),str(year),str(month),str(xoge), str(yoge)])
            plt.savefig(f'particles_out/particle_{particle}_year_{year}_month_{month}.png')

            #fig.savefig('C:\ECCOv4\output108\\' + str(counter) + '__' + str(year) +'_' + str("{:02d}".format(month+1)) + '_particle_' + str("{:02d}".format(particlecounter)) + '.png')
            #logging.debug(str(counter),str(particle),str(year),str("{:02d}".format(month+1)),str(xoge), str(yoge))
            xoge, yoge = move_1month(xoge, yoge, uvel, vvel, fudge=False)
            if outOfTile(xoge, yoge):
                return monthly
    return monthly
                    
# Above counter increments the particle count                




results = []
#input_file = 'klawinput.csv'
input_file = 'testing.csv'
with open(input_file) as csv_file:
    csv_reader = csv.DictReader(csv_file, delimiter=',')
    particle = 0
    for row in csv_reader:
        #logging.debug(particle)
        xoge = int(row['xoge'])
        yoge = int(row['yoge'])

        # With recursive_load_ecco_var_from_years_nc one can specify one or more variables to load, one or more years to load, while also requesting only a subset of tiles and vertical levels.
      # ecco.recursive_load_ecco_var_from_years_nc(d
        particle_result = particle_positions(particle, xoge, yoge, range(1992,2015))
        results.extend(particle_result)

        particle += 1

with open('particles_out/' + 'eccodataset_output.csv', mode='w+', newline='') as out_file:
    out_writer = csv.writer(out_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for result in results:
        out_writer.writerow(result)
