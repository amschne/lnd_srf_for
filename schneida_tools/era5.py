#!/usr/bin/env python

''' Map near surface meteorological variables from the fifth
    generation of the European Centre for Medium-Range Weather Forecasts (ECMWF)
    atmospheric reanalyses (ERA5) (Delhasse et al., 2020).

    References:
    Delhasse, A., Kittel, C., Amory, C., Hofer, S., van As, D., S. Fausto, R.,
        and Fettweis, X.: Brief communication: Evaluation of the near-surface
        climate in ERA5 over the Greenland Ice Sheet, The Cryosphere, 14,
        957–965, https://doi.org/10.5194/tc-14-957-2020, 2020.
'''

import os
import argparse
import fnmatch
from netCDF4 import Dataset
import numpy as np

from schneida_tools import ncks_mk_time_rec_dmn
from schneida_tools.schneida_args import get_args

import ipdb

class ERA5(object):
    def __init__(self):
        args = get_args()
        self.clean_data_path = args.era5_clean_data_path
        self.sorted_files = sorted(os.listdir(self.clean_data_path))
    
    def get_data(self, pattern):
        """ Use fnmatch to find relevant data files and append to list
        
            Return data_list, a list of ERA5 (netCDF4) data instances
        """
        rootgrps = list()
        print('Opening %s' % os.path.join(self.clean_data_path, pattern))
        for i, file_name in enumerate(self.sorted_files):
            if fnmatch.fnmatch(file_name, pattern):
                rootgrps.append(Dataset(os.path.join(self.clean_data_path,
                                                     file_name)))
        return rootgrps
    
    '''
    def get_lwdown(self, pattern='LWdown_ERA5_CRU_*_v1.1.nc'):
        """ Get downwelling longwave radiation
        """
        self.lw_down = self.get_data(pattern)
        
    def get_psurf(self, pattern='PSurf_ERA5_CRU_*_v1.1.nc'):
        """ Get surface pressure
        """
        self.p_surf = self.get_data(pattern)
    
    def get_qair(self, pattern='Qair_ERA5_CRU_*_v1.1.nc'):
        """ Get humidity
        """
        self.q_air = self.get_data(pattern)
        
    def get_rainf(self, pattern='Rainf_ERA5_CRU+GPCC_*_v1.1.nc'):
        """ Get rainfall
        """
        self.rainf = self.get_data(pattern)
        
    def get_snowf(self, pattern='Snowf_ERA5_CRU+GPCC_*_v1.1.nc'):
        """ Get snowfall
        """
        self.snowf = self.get_data(pattern)
    
    def get_swdown(self, pattern='SWdown_ERA5_CRU_*_v1.1.nc'):
        """ Get downwelling shortwave radiation
        """
        self.sw_down = self.get_data(pattern)
    '''    
    def get_tair(self, pattern='adaptor.mars.internal*.nc'):
        """ Get air temperature
        """
        self.t_air = self.get_data(pattern)
    def get_precip(self, pattern='adaptor.mars.internal*.nc'):
        """ Get total precipitation
        """
        self.precip = self.get_data(pattern)
    '''
    def get_wind(self, pattern='Wind_ERA5_CRU_*_v1.1.nc'):
        """ Get wind
        """
        self.wind = self.get_data(pattern)
    '''
        
def get_temporal_mean(rootgrp_list, var, compute=True, results_dir='results'):
    """
    """
    file_path = os.path.join(results_dir,
                             '%s_ERA5_temporal_mean.nc' % var)
    if compute: # Create new netCDF file
        rootgrp_out = Dataset(file_path, 'w', format='NETCDF4')
        
        era5_init_mon = rootgrp_list[0]
        # Create dimensions
        lat = rootgrp_out.createDimension('lat',
                                          len(era5_init_mon.dimensions['latitude']))
        lon = rootgrp_out.createDimension('lon',
                                          len(era5_init_mon.dimensions['longitude']))
        # Create variables
        latitudes = rootgrp_out.createVariable('lat', 'f8', ('lat'))
        longitudes = rootgrp_out.createVariable('lon', 'f8', ('lon'))
        mean_var = rootgrp_out.createVariable(var, 'f8', ('lat', 'lon'))
        
        # Add attributes
        #latitudes.standard_name = era5_init_mon.variables['latitude'].standard_name
        latitudes.units = era5_init_mon.variables['latitude'].units
        #latitudes.axis = era5_init_mon.variables['latitude'].axis
        latitudes.long_name = era5_init_mon.variables['latitude'].long_name
        
        #longitudes.standard_name = era5_init_mon.variables['longitude'].standard_name
        longitudes.units = era5_init_mon.variables['longitude'].units
        #longitudes.axis = era5_init_mon.variables['longitude'].axis
        longitudes.long_name = era5_init_mon.variables['longitude'].long_name
        
        mean_var.long_name = era5_init_mon.variables[var].long_name
        #mean_var.standard_name = era5_init_mon.variables[var].standard_name
        mean_var.units = era5_init_mon.variables[var].units
        
        # Write data
        latitudes[:] = era5_init_mon.variables['latitude'][:]
        longitudes[:] = era5_init_mon.variables['longitude'][:]
        
        print('Computing ERA5 temporal mean "%s"...' % var)
        # Initialize temporal mean numpy masked array
        era5_time_sum = np.ma.zeros(era5_init_mon.variables[var][0].shape)
        for i, era5_month in enumerate(rootgrp_list):
            era5_time_sum += np.ma.mean(era5_month.variables[var][:], axis=0)
        mean_var[:,:] = era5_time_sum / float(i + 1)
        
        rootgrp_out.close()
    
    # Open and return temporal mean dataset
    return(Dataset(file_path, 'r'))

def close_rootgrps(rootgrps):
    """ Close files in ERA5 rootgrp list
    """
    for i, rootgrp in enumerate(rootgrps):
        rootgrp.close()

def clean_data():
    """ 1. Gather raw data
        2. Call ncks to change the time dimension to the record dimension.
    """
    args = get_args()
    
    era5_files = os.listdir(args.era5_raw_data_path)
    ncks_mk_time_rec_dmn.call_ncks(args.era5_raw_data_path, era5_files,
                                   args.era5_clean_data_path)

def test():
    data = ERA5()
    '''
    data.get_lwdown()
    close_rootgrps(data.lw_down)
    
    data.get_psurf()
    close_rootgrps(data.p_surf)
    
    data.get_qair()
    close_rootgrps(data.q_air)
    
    data.get_rainf()
    close_rootgrps(data.rainf)
    
    data.get_snowf()
    close_rootgrps(data.snowf)
    
    data.get_swdown()
    close_rootgrps(data.sw_down)
    '''
    data.get_tair()
    close_rootgrps(data.t_air)
    data.get_precip()
    close_rootgprs(data.precip)
    '''
    data.get_wind()
    close_rootgrps(data.wind)
    '''
def run():
    #clean_data()
    #test()
    data = ERA5()
    data.get_tair()
    mean_temperature = get_temporal_mean(data.t_air, 't2m', compute=True)
    close_rootgrps(data.t_air)
    data.get_precip()
    mean_precip = get_temporal_mean(data.precip, 'tp', compute=True)
    close_rootgrps(data.precip)
    
    ipdb.set_trace()
    
def main():
    run()

if __name__=='__main__':
    main()
