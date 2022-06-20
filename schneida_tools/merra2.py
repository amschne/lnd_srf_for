#!/usr/bin/env python

"""
"""

import os
import fnmatch

import numpy as np
from netCDF4 import Dataset

from schneida_args import get_args

import ipdb

SURF_AIR_TEMP= 'TLML'
TOT_PREC = 'PRECTOT'

class MERRA2(object):
    def __init__(self, first_year=1980, last_year=1989):
        args = get_args()
        self.raw_data_path = args.merra2_raw_data_path
        self.sorted_files = dict()
        for i, year in enumerate(np.arange(first_year, last_year + 1)):
            self.sorted_files[year] = sorted(os.listdir(os.path.join(
                                                             self.raw_data_path,
                                                             '%d' % year)))

    def get_data(self, pattern):
        """ Use fnmatch to find relevant data files and append to list
        
            Return data_list, a list of MERRA-2 (netCDF4) data instances
        """
        rootgrps = list()
        for year, sorted_files in self.sorted_files.items():
            print('Opening %s' % os.path.join(self.raw_data_path, '%d' % year,
                                              pattern))
            for j, file_name in enumerate(sorted_files):
                if fnmatch.fnmatch(file_name, pattern):
                    rootgrps.append(Dataset(os.path.join(self.raw_data_path,
                                                         '%d' % year,
                                                         file_name)))
        return rootgrps
    
    def get_tair(self, pattern='MERRA2_100.tavgM_2d_flx_Nx.*.nc4'):
        
        self.t_air = self.get_data(pattern)
        
    def get_precip(self, pattern='MERRA2_100.tavgM_2d_flx_Nx.*.nc4'):
        """ Get total precipitation
        """
        self.precip = self.get_data(pattern)

def get_temporal_mean(rootgrp_list, var, compute=True, results_dir='results'):
    """
    """
    file_path = os.path.join(results_dir,
                             '%s_MERRA2_temporal_mean.nc' % var)
    if compute: # Create new netCDF file
        rootgrp_out = Dataset(file_path, 'w', format='NETCDF4')
        
        merra2_init_mon = rootgrp_list[0]
        # Create dimensions
        lat = rootgrp_out.createDimension('lat',
                                          len(merra2_init_mon.dimensions['lat']))
        lon = rootgrp_out.createDimension('lon',
                                          len(merra2_init_mon.dimensions['lon']))
        # Create variables
        latitudes = rootgrp_out.createVariable('lat', 'f8', ('lat'))
        longitudes = rootgrp_out.createVariable('lon', 'f8', ('lon'))
        mean_var = rootgrp_out.createVariable(var, 'f8', ('lat', 'lon'))
        
        # Add attributes
        latitudes.units = merra2_init_mon.variables['lat'].units
        latitudes.long_name = merra2_init_mon.variables['lat'].long_name
        
        longitudes.units = merra2_init_mon.variables['lon'].units
        longitudes.long_name = merra2_init_mon.variables['lon'].long_name
        
        mean_var.long_name = merra2_init_mon.variables[var].long_name
        mean_var.units = merra2_init_mon.variables[var].units
        
        # Write data
        latitudes[:] = merra2_init_mon.variables['lat'][:]
        longitudes[:] = merra2_init_mon.variables['lon'][:]
        
        print('Computing MERRA-2 temporal mean "%s"...' % var)
        # Initialize temporal mean numpy masked array
        merra2_time_sum = np.ma.zeros(merra2_init_mon.variables[var][0].shape)
        for i, merra2_month in enumerate(rootgrp_list):
            merra2_time_sum += np.ma.mean(merra2_month.variables[var][:], axis=0)
        mean_var[:,:] = merra2_time_sum / float(i + 1)
        
        rootgrp_out.close()
    
    # Open and return temporal mean dataset
    return(Dataset(file_path, 'r'))

def close_rootgrps(rootgrps):
    """ Close files in MERRA-2 rootgrp list
    """
    for i, rootgrp in enumerate(rootgrps):
        rootgrp.close()

def run():
    #clean_data()
    #test()
    data = MERRA2()
    data.get_tair()
    mean_temperature = get_temporal_mean(data.t_air, SURF_AIR_TEMP, compute=True)
    close_rootgrps(data.t_air)
    
    data.get_precip()
    mean_precip = get_temporal_mean(data.precip, TOT_PREC, compute=True)
    close_rootgrps(data.precip)
    
def main():
    run()

if __name__=='__main__':
    main()