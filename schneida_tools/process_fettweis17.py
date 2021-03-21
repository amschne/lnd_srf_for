#!/usr/bin/env python

""" Process surface mass balance (SMB) results from Fettweis et al., 2017

    Writes Portable Network Graphics files that map the Greenland Ice
    Sheet SMB
"""

import os
import argparse
import fnmatch

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import colorcet as cc
import cartopy.feature as cfeature
import cartopy.crs as ccrs

from schneida_tools import gris_dem

SECTOR = 0

def contourf_smb(dates, lat, lon, smb_mm_day):
    #smb_m_yr = (365.25 * smb_mm_day[:]) / 1000.
    dem = gris_dem.GrISDEM(os.path.join('data_raw','gimpdem_90m_v01.1.tif'))
    for t, date in enumerate(dates):
        
        #fig, ax = plt.subplots()
        mp = dem.setup_map()
        mp.add_feature(cfeature.LAND, color='white')#'#555759')
        mp.add_feature(cfeature.OCEAN, color='white')#'#555759')
        if int(date) != 2012121088 and int(date) != 2012110976 and int(date) != 2012020736 and int(date) != 2012022016 and int(date) != 2012022784: # problem spots
            coll = mp.contourf(lon[:], lat[:],
                           np.clip(smb_mm_day[t,SECTOR], -7.99, 7.99),
                           levels=np.arange(-8, 8.1, 1),
                           cmap='cet_CET_D10',
                           vmin=-8, vmax=8,
                           transform=ccrs.PlateCarree())
            dem.draw_contours(mp, os.path.join('data_clean', "gimpdem_90m_v01.1.nc"))

            '''
            mp.readshapefile('greenland_coast', 'coast',
                     color='#C6BEB5')
            mp.drawmapscale(-36, 62, -45, 70, 20, barstyle='simple', units='km',
                        fontsize=9, yoffset=None, labelstyle='simple',
                        fontcolor='white', fillcolor1='k', fillcolor2='w',
                        linewidth=0.5)
            '''             
            fig = plt.gcf()
            cbar = fig.colorbar(coll, ax=mp,
                            orientation = 'horizontal',
                            ticks=np.arange(-8, 9, 2),
                            boundaries=np.arange(-8,9,1),
                            values=np.arange(-7.5, 8, 1))
            cbar.set_label('SMB (mm w.eq. / day)')
            year = str(date)[:4]
            mon = str(date)[4:6]
            plt.title('%s-%s' % (year, mon))
         
            plt.savefig(os.path.join('results','graphics','smb_%s_sector%d.png')
                    % (date, SECTOR), dpi=300)
      
        plt.close()

def get_data(input_dir, start_year, stop_year):
    data_stream = list()
    anl_years = np.arange(start_year, stop_year)
    for i, file_name in enumerate(sorted(os.listdir(input_dir))):
        if fnmatch.fnmatch(file_name, 'ICE.*.01-12.nc'):
            # Data stream file
            year = np.array([int(file_name.split('.')[1])])
            if np.isin(year, anl_years)[0]:
                data_stream.append(os.path.join(input_dir, file_name))
                
    return data_stream

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", help='Path to data from '
                                          'Fettweis et al., 2017')
    parser.add_argument("--start_year", help='YYYY', type=int,
                        default=2014)
    parser.add_argument("--stop_year", help='YYYY', type=int,
                        default=2017)
    args = parser.parse_args()
    return args

def run(verbose=True):
    plt.style.use('agu_half_vertical')
    plt.style.use('uci_darkblue')
    args = get_args()
    data_stream = get_data(args.input_dir, args.start_year, args.stop_year)
    
    for i, fi in enumerate(data_stream):
        f = Dataset(fi)
        if verbose:
            print(f)
            print(f.variables.keys())
        date = f.variables['DATE']
        sector = f.variables['SECTOR']
        smb = f.variables['SMB']
        lat = f.variables['LAT']
        lon = f.variables['LON']
        mask = f.variables['MSK']
        print(date)
        print(smb)
        print(sector)
        print(mask)
        if verbose:
            for d in f.dimensions.items():
                print(d)
        gis_smb = 0.01 * mask[:] * smb[:]
        contourf_smb(date, lat, lon, gis_smb)
        
        f.close()

def main():
    run()

if __name__=='__main__':
    main()
