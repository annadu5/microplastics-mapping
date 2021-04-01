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
import logging

logging.basicConfig(level=logging.DEBUG)

base_dir = '/eccodata'
grid_dir = f'{base_dir}/Version4/Release3_alt/nctiles_grid'
ecco_grid = ecco.load_ecco_grid_nc(grid_dir, 'ECCOv4r3_grid.nc')
day_mean_dir = f'{base_dir}/Version4/Release3_alt/nctiles_monthly'

counter = 0
def particle_positions(particle, xoge0, yoge0, year_range, tile=10):
    monthly = []
    xoge = xoge0
    yoge = yoge0
    global counter
    for year in year_range:  # good stretch is 1999 to 2009  BELOW UNDER K and TILES NEED TO MAKE VARS
        ecco_vars = ecco.recursive_load_ecco_var_from_years_nc(
                day_mean_dir, vars_to_load=['UVEL', 'VVEL'], years_to_load=year, dask_chunk=False)
                #k_subset = [0], \
                #tiles_to_load = [10], \

        for month in range (0, 12):   # 12 vs 11
            counter += 1
            
            logging.debug(" ")
            logging.debug(" ")
            logging.debug("--------")
            logging.debug("PARTICLE : " + str(particle) + " STARTING COORD : " + str(particle) + "  YEAR : " + str(year) + "  MONTH : " + str(month+1))

            #ecco_vars = \
            #    ecco.recursive_load_ecco_var_from_years_nc(day_mean_dir,\
            #                                                vars_to_load=['UVEL', 'VVEL'], \
            #                                                years_to_load=year,\
            #                                                dask_chunk=False)
            

            
            ecco_ds = xr.merge((ecco_grid , ecco_vars))

            ecco_ds.attrs = []

            uvel_dataset = xr.merge((ecco_grid , ecco_vars))
            uvel_dataset.attrs = []

            uvel_A = ecco_ds.UVEL
            uvel_B = ecco_ds['UVEL']
            uvel_arr = uvel_A.values

            vvel_A = ecco_ds.VVEL
            vvel_B = ecco_ds['VVEL']
            vvel_arr = vvel_A.values    

            #-------------------
            # fig=plt.figure(figsize=(9, 9))
            # tile_num=10
            # # pull out lats and lons
            # lons = np.copy(ecco_ds.XC.sel(tile=tile_num))

            # lons[lons < 0] = lons[lons < 0]+360
            # lats = ecco_ds.YC.sel(tile=tile_num)
            # tile_to_plot = ecco_ds.UVEL.isel(tile=tile_num, time=month, k=0)


            # tile_to_plot= tile_to_plot.where(ecco_ds.hFacW.isel(tile=tile_num,k=0) !=0, np.nan)

            # plt.pcolor(lons, lats, tile_to_plot, vmin=-.25, vmax=.25, cmap='jet')
            # plt.colorbar()
            #-----------------------

            # the lines commented back here work fine! (below)

            # fig=plt.figure(figsize=(10, 8.5),dpi=90)
            # plt. layout='latlon'
            # ud_masked = uvel_dataset.UVEL.where(uvel_dataset.hFacW > 0, np.nan)

            # ud_masked.isel(k=0,tile=10, time=month).plot(cmap='jet', vmin=-.26,vmax=.26)   #, vmin=-.1,vmax=.1)

            # ecco_plt_tiles(ud_masked(layout='latlon', rotate_to_latlon=True)
            # plt.title('ECCO v4r3 Velocity Chart by Anna Du');
            # plt.title(' ');




            alist = uvel_A.values[month,tile,0,int(yoge),int(xoge)] # Here the first threeo of these are correct
            blist = vvel_A.values[month,tile,0,int(yoge),int(xoge)] # m/s needs to be converted into a distance -- this is a VELOCITY

            # HERE NEED TO CONVERT ABOVE INTO A DISTANCE TO ADD TO DOTCXOGE DOTYOGE
            # 2592000 seconds in one month
            # Now need to convert this into a real distance that can be plotted on a 90 x 90 grid
            # one latitude line is 69 miles x 90 pts
            # the earth is 40,075 km around at the equator aka circumference == aka linear distance of the 360 degree latitude in ECCO
            # so divide this by 4 (4x90=360)
            # so that is 40075 / 4 = 10,018.75 km (x 1000 m ) = 10,018,750 is the amt of meters in one 90 degree ECCO tile map
            # so you need to divide this number by 90 and then you get how many meters are in each one of the 90 array pixels = 111,319.4444
            # this number is the amount of meters in one pixel
            # as an example lets multoply .086 times the amt of seconds in a month (2592000) = 222,912 = meters that this moves in ONE MONTH
            # so it is reasonable to assume that we can move several pixels in one month so in a year we can move 20 pixels or more, maybe dozens
            # LOGIC = you take the velocity of UVEL or VVEL (alist or blist) multiply it by 
            #
                    
 
            # alist = uvel_A.values[month,10,0,int(xoge),int(yoge)]   # old way of doing this
            # blist = vvel_A.values[month,10,0,int(xoge),int(yoge)]

            logging.debug(f'alist={alist}, blist={blist}')

            dotxoge = float(xoge + alist*23.2843418683) #*fudgefactor)  THIS IS ADDING A VELOCITY TO A DISTNACE - NEED TO REDO THE THINGY ABOVE
            dotyoge = float(yoge + blist*23.2843418683)

            logging.debug(f'dotxoge={dotxoge}, dotyoge={dotyoge}')
            
            xoge = int(dotxoge) 
            yoge = int(dotyoge)
                
            # if blist == 0.0:
            #     blist = vvel_A.values[month,10,0,int(xoge),int(yoge)]

            #if (yoge > 10):
            #    plt.scatter(xoge,yoge,color='black')
            # KEEEP THIS EXAMPLE BELOW
            logging.debug(f"ecco_ds.hFacW.dims: {ecco_ds.hFacW.dims}")
            

            if abs(alist) == 0 and abs(blist) == 0:
                    logging.debug("NOT MOVING")
                    plt.scatter(dotxoge,dotyoge,color='black',marker='D',s=40)         
                                        
            elif abs(dotyoge) <= 0:
                    logging.debug("stuck on bottom")
                    plt.scatter(dotxoge,dotyoge,color='yellow',marker='D',s=40)
                    return monthly
            elif abs(dotyoge) >= 90:
                    logging.debug("stuck on top")
                    plt.scatter(dotxoge,dotyoge,color='yellow',marker='D',s=40)
                    return monthly

            elif abs(dotxoge) <= 0:
                    logging.debug("stuck on leftside")
                    plt.scatter(dotxoge,dotyoge,color='yellow',marker='D',s=40)
                    return monthly
            elif abs(dotxoge) >= 90:
                    logging.debug("stuck on rightside")
                    plt.scatter(dotxoge,dotyoge,color='yellow',marker='D',s=40)
                    return monthly

            else:
                # THIS ONE IS NOW RIGHT# was backwards  [k,tile,j,i_g] (time, tile, k j i 
                beachhfacw = ecco_ds.hFacW.values[0,tile,int(yoge),int(xoge)]
                logging.debug(beachhfacw)
                if (int(beachhfacw) != 0):
                    plt.scatter(dotxoge,dotyoge,color='magenta')
                    logging.debug("notbeached")
                else:
                    logging.debug("beached!")
                    plt.scatter(dotxoge,dotyoge,color='black',marker='D',s=40)


            monthly.append([counter, particle, year, month, xoge, yoge])
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
        particle_result = particle_positions(particle, xoge, yoge, range(1999,2009))
        results.extend(particle_result)

        particle += 1

with open('particles_out/' + 'eccodataset_output.csv', mode='w+', newline='') as out_file:
    out_writer = csv.writer(out_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for result in results:
        out_writer.writerow(result)
