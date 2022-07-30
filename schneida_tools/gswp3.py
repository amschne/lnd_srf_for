#!/usr/bin/env python

""" Map near surface meteorological variables from the Global Soil Wetness
    Project Phase 3 (GSWP3) atmospheric forcing data (Kim, 2017).

    References:
    Hyungjun Kim. (2017). Global Soil Wetness Project Phase 3 Atmospheric
        Boundary Conditions (Experiment 1) [Data set]. Data Integration and
        Analysis System (DIAS). https://doi.org/10.20783/DIAS.501
"""

import os
import argparse

from netCDF4 import MFDataset
from netCDF4 import Dataset
import numpy as np
import ncks_mk_time_rec_dmn
from schneida_args import get_args

import ipdb

class GSWP3(object):
    def __init__(self, file_prefix='clmforc.GSWP3.c2011.0.5x0.5.'):
        self.file_prefix = file_prefix
        
        # Get program arguments
        args = get_args()
        self.clean_data_path = args.gswp3_clean_data_path
        self.precip_dir = args.gswp3_precip_dir
        self.solar_dir = args.gswp3_solar_dir
        self.tphwl_dir = args.gswp3_tphwl_dir
        
    def get_precip(self):
        """ Get precipitation
        """
        pattern = self.file_prefix + 'Prec.*.nc'
        files = os.path.join(self.clean_data_path, self.precip_dir, pattern)
        print('Opening %s' % files)
        self.precip_rootgrp = MFDataset(files)
        
    def get_solar(self):
        """ Get downwelling solar radiation
        """
        pattern = self.file_prefix + 'Solr.*.nc'
        files = os.path.join(self.clean_data_path, self.solar_dir, pattern)
        print('Opening %s' % files)
        self.solar_rootgrp = MFDataset(files)
        
    def get_tphwl(self):
        """ Get temperature, pressure, humidity, wind, and downwelling longwave
            radiation
        """
        pattern = self.file_prefix + 'TPQWL.*.nc'
        files = os.path.join(self.clean_data_path, self.tphwl_dir, pattern)
        print('Opening %s' % files) 
        self.tphwl_rootgrp = MFDataset(files)

def get_temporal_mean(rootgrp_in, var, compute=True, results_dir='results'):
    """ Works for var = 'TBOT'
    """
    # Get high level args
    data = GSWP3()
    
    file_path = os.path.join(results_dir,
                             data.file_prefix + var + '.temporal_mean.nc')
    if compute: # Create new netCDF file
        rootgrp_out = Dataset(file_path, 'w', format='NETCDF4')
        # Create dimensions
        lat = rootgrp_out.createDimension('lat',
                                          len(rootgrp_in.dimensions['lat']))
        lon = rootgrp_out.createDimension('lon',
                                          len(rootgrp_in.dimensions['lon']))
        # Create variables
        latitudes = rootgrp_out.createVariable('LATIXY', 'f8', ('lat', 'lon'))
        longitudes = rootgrp_out.createVariable('LONGXY', 'f8', ('lat', 'lon'))
        mean_var = rootgrp_out.createVariable(var, 'f8', ('lat', 'lon'))
        
        # Add attributes
        latitudes.long_name = rootgrp_in.variables['LATIXY'].long_name
        latitudes.units = rootgrp_in.variables['LATIXY'].units
        latitudes.mode = rootgrp_in.variables['LATIXY'].mode
        
        longitudes.long_name = rootgrp_in.variables['LONGXY'].long_name
        longitudes.units = rootgrp_in.variables['LONGXY'].units
        longitudes.mode = rootgrp_in.variables['LONGXY'].mode
        
        mean_var.long_name = rootgrp_in.variables[var].long_name
        mean_var.units = rootgrp_in.variables[var].units
        mean_var.mode = rootgrp_in.variables[var].mode
        mean_var.missing_value = rootgrp_in.variables[var].missing_value
        
        # Write data
        latitudes[:,:] = rootgrp_in.variables['LATIXY'][:]
        longitudes[:,:] = rootgrp_in.variables['LONGXY'][:]
        
        print('Computing GSWP3 temporal mean "%s"...' % var)
        mean_var[:,:] = np.ma.mean(rootgrp_in.variables[var][:], axis=0)
        
        rootgrp_out.close()
    
    # Open and return temporal mean dataset
    return(Dataset(file_path, 'r'))
        
def clean_data():
    """ 1. Gather raw data
        2. Call ncks to change the time dimension to the record dimension.
    """
    args = get_args()
    
    # Clean precipiation data
    precip_path = os.path.join(args.gswp3_raw_data_path,
                               args.gswp3_precip_dir)
    precip_files = os.listdir(precip_path)
    ncks_mk_time_rec_dmn.call_ncks(precip_path, precip_files,
                                   os.path.join(args.gswp3_clean_data_path,
                                   args.gswp3_precip_dir))
    
    # Clean solar data
    solar_path = os.path.join(args.gswp3_raw_data_path,
                              args.gswp3_solar_dir)
    solar_files = os.listdir(solar_path)
    ncks_mk_time_rec_dmn.call_ncks(solar_path, solar_files,
                                   os.path.join(args.gswp3_clean_data_path,
                                   args.gswp3_solar_dir))
    
    # Clean temperature, pressure, humidity, wind, and longwave data
    tphwl_path = os.path.join(args.gswp3_raw_data_path,
                              args.gswp3_tphwl_dir)
    tphwl_files = os.listdir(tphwl_path)
    ncks_mk_time_rec_dmn.call_ncks(tphwl_path, tphwl_files,
                                   os.path.join(args.gswp3_clean_data_path,
                                   args.gswp3_tphwl_dir))

def test():
    data = GSWP3()
    
    data.get_precip()
    data.precip_rootgrp.close()
    
    data.get_solar()
    data.solar_rootgrp.close()
    
    data.get_tphwl()
    data.tphwl_rootgrp.close()

def run():
    #clean_data()
    #test()
    data = GSWP3()
    data.get_tphwl()
    mean_temperature = get_temporal_mean(data.tphwl_rootgrp, 'TBOT', compute=True)
    data.tphwl_rootgrp.close()
    
    data.get_precip()
    mean_precip = get_temporal_mean(data.precip_rootgrp, 'PRECTmms', compute=True)
    data.precip_rootgrp.close()
    
    ipdb.set_trace()
    
def main():
    run()

if __name__=='__main__':
    main()
