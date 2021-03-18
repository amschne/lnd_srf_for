#!/usr/bin/env python

''' Map near surface meteorological variables from the "water, energy, and
    climate change" (WATCH) forcing data methodology applied to the fifth
    generation of the European Centre for Medium-Range Weather Forecasts (ECMWF)
    atmospheric reanalyses (ERA5), aka the WFDE5 (Cucchi et al., 2020).

    References:
    Cucchi, M., Weedon, G. P., Amici, A., Bellouin, N., Lange, S., Müller
        Schmied, H., Hersbach, H., and Buontempo, C.: WFDE5: bias-adjusted ERA5
        reanalysis data for impact studies, Earth Syst. Sci. Data, 12,
        2097–2120, https://doi.org/10.5194/essd-12-2097-2020, 2020.
'''

import os
import argparse
import fnmatch
from netCDF4 import Dataset
import numpy as np

from schneida_tools import ncks_mk_time_rec_dmn
from schneida_tools.schneida_args import get_args

class WFDE5(object):
    def __init__(self):
        args = get_args()
        self.clean_data_path = args.wfde5_clean_data_path
        self.sorted_files = sorted(os.listdir(self.clean_data_path))
    
    def get_data(self, pattern):
        """ Use fnmatch to find relevant data files and append to list
        
            Return data_list, a list of WFDE5 (netCDF4) data instances
        """
        rootgrps = list()
        print('Opening %s' % os.path.join(self.clean_data_path, pattern))
        for i, file_name in enumerate(self.sorted_files):
            if fnmatch.fnmatch(file_name, pattern):
                rootgrps.append(Dataset(os.path.join(self.clean_data_path,
                                                     file_name)))
        return rootgrps
    
    def get_lwdown(self, pattern='LWdown_WFDE5_CRU_*_v1.1.nc'):
        """ Get downwelling longwave radiation
        """
        self.lw_down = self.get_data(pattern)
        
    def get_psurf(self, pattern='PSurf_WFDE5_CRU_*_v1.1.nc'):
        """ Get surface pressure
        """
        self.p_surf = self.get_data(pattern)
    
    def get_qair(self, pattern='Qair_WFDE5_CRU_*_v1.1.nc'):
        """ Get humidity
        """
        self.q_air = self.get_data(pattern)
        
    def get_rainf(self, pattern='Rainf_WFDE5_CRU+GPCC_*_v1.1.nc'):
        """ Get rainfall
        """
        self.rainf = self.get_data(pattern)
        
    def get_snowf(self, pattern='Snowf_WFDE5_CRU+GPCC_*_v1.1.nc'):
        """ Get snowfall
        """
        self.snowf = self.get_data(pattern)
    
    def get_swdown(self, pattern='SWdown_WFDE5_CRU_*_v1.1.nc'):
        """ Get downwelling shortwave radiation
        """
        self.sw_down = self.get_data(pattern)
        
    def get_tair(self, pattern='Tair_WFDE5_CRU_*_v1.1.nc'):
        """ Get air temperature
        """
        self.t_air = self.get_data(pattern)
    
    def get_wind(self, pattern='Wind_WFDE5_CRU_*_v1.1.nc'):
        """ Get wind
        """
        self.wind = self.get_data(pattern)
        
def get_temporal_mean(rootgrp_list, var, compute=True, results_dir='results'):
    """
    """
    file_path = os.path.join(results_dir,
                             '%s_WFDE5_temporal_mean_v1.1.nc' % var)
    if compute: # Create new netCDF file
        rootgrp_out = Dataset(file_path, 'w', format='NETCDF4')
        
        wfde5_init_mon = rootgrp_list[0]
        # Create dimensions
        lat = rootgrp_out.createDimension('lat',
                                          len(wfde5_init_mon.dimensions['lat']))
        lon = rootgrp_out.createDimension('lon',
                                          len(wfde5_init_mon.dimensions['lon']))
        # Create variables
        latitudes = rootgrp_out.createVariable('lat', 'f8', ('lat'))
        longitudes = rootgrp_out.createVariable('lon', 'f8', ('lon'))
        mean_var = rootgrp_out.createVariable(var, 'f8', ('lat', 'lon'))
        
        # Add attributes
        latitudes.standard_name = wfde5_init_mon.variables['lat'].standard_name
        latitudes.units = wfde5_init_mon.variables['lat'].units
        latitudes.axis = wfde5_init_mon.variables['lat'].axis
        latitudes.long_name = wfde5_init_mon.variables['lat'].long_name
        
        longitudes.standard_name = wfde5_init_mon.variables['lon'].standard_name
        longitudes.units = wfde5_init_mon.variables['lon'].units
        longitudes.axis = wfde5_init_mon.variables['lon'].axis
        longitudes.long_name = wfde5_init_mon.variables['lon'].long_name
        
        mean_var.long_name = wfde5_init_mon.variables[var].long_name
        mean_var.standard_name = wfde5_init_mon.variables[var].standard_name
        mean_var.units = wfde5_init_mon.variables[var].units
        
        # Write data
        latitudes[:] = wfde5_init_mon.variables['lat'][:]
        longitudes[:] = wfde5_init_mon.variables['lon'][:]
        
        print('Computing WFDE5 temporal mean "%s"...' % var)
        # Initialize temporal mean numpy masked array
        wfde5_time_sum = np.ma.zeros(wfde5_init_mon.variables[var][0].shape)
        for i, wfde5_month in enumerate(rootgrp_list):
            wfde5_time_sum += np.ma.mean(wfde5_month.variables[var][:], axis=0)
        mean_var[:,:] = wfde5_time_sum / float(i + 1)
        
        rootgrp_out.close()
    
    # Open and return temporal mean dataset
    return(Dataset(file_path, 'r'))

def close_rootgrps(rootgrps):
    """ Close files in WFDE5 rootgrp list
    """
    for i, rootgrp in enumerate(rootgrps):
        rootgrp.close()

def clean_data():
    """ 1. Gather raw data
        2. Call ncks to change the time dimension to the record dimension.
    """
    args = get_args()
    
    wfde5_files = os.listdir(args.wfde5_raw_data_path)
    ncks_mk_time_rec_dmn.call_ncks(args.wfde5_raw_data_path, wfde5_files,
                                   args.wfde5_clean_data_path)

def test():
    data = WFDE5()
    
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
    
    data.get_tair()
    close_rootgrps(data.t_air)
    
    data.get_wind()
    close_rootgrps(data.wind)

def run():
    #clean_data()
    #test()
    data = WFDE5()
    data.get_tair()
    mean_temperature = get_temporal_mean(data.t_air, 'Tair', compute=True)
    close_rootgrps(data.t_air)
    
    data.get_rainf()
    mean_rainf = get_temporal_mean(data.rainf, 'Rainf', compute = True)
    close_rootgrps(data.rainf)
    
    data.get_snowf()
    mean_snowf = get_temporal_mean(data.snowf, 'Snowf', compute=True)
    close_rootgrps(data.snowf)
    
    ipdb.set_trace()
    
def main():
    run()

if __name__=='__main__':
    main()
