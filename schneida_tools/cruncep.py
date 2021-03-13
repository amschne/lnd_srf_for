#!/usr/bin/env python

""" Map near surface meteorological variables from the Climate Research Unit /
    National Centers for Environmental Prediction (Version 7) atmospheric
    forcing data (Viovy, 2018).

    References:
    Viovy, N. (2018), CRUNCEP Version 7 - Atmospheric Forcing Data for the
        Community Land Model, https://doi.org/10.5065/PZ8F-F017, Research
        Data Archive at the National Center for Atmospheric Research,
        Computational and Information Systems Laboratory, Boulder, Colo.
        Accessed 04 Mar 2021.
"""

import os
import argparse
from netCDF4 import MFDataset

from schneida_tools import ncks_mk_time_rec_dmn
from schneida_tools.schneida_args import get_args

class CRUNCEP7(object):
    def __init__(self):
        args = get_args()
        self.clean_data_path = args.cruncep_clean_data_path
        self.precip_dir = args.cruncep_precip_dir
        self.solar_dir = args.cruncep_solar_dir
        self.tphwl_dir = args.cruncep_tphwl_dir
        
    def get_precip(self, pattern='clmforc.cruncep.V7.c2016.0.5d.Prec.*.nc'):
        """ Get precipitation
        """
        files = os.path.join(self.clean_data_path, self.precip_dir, pattern)
        print('Opening %s' % files)
        self.precip_rootgrp = MFDataset(files)
        
    def get_solar(self, pattern='clmforc.cruncep.V7.c2016.0.5d.Solr.*.nc'):
        """ Get downwelling solar radiation
        """
        files = os.path.join(self.clean_data_path, self.solar_dir, pattern)
        print('Opening %s' % files)
        self.solar_rootgrp = MFDataset(files)
        
    def get_tphwl(self, pattern='clmforc.cruncep.V7.c2016.0.5d.TPQWL.*.nc'):
        """ Get temperature, pressure, humidity, wind, and downwelling longwave
            radiation
        """
        files = os.path.join(self.clean_data_path, self.tphwl_dir, pattern)
        print('Opening %s' % files) 
        self.tphwl_rootgrp = MFDataset(files)

def clean_data():
    """ 1. Gather raw data
        2. Call ncks to change the time dimension to the record dimension.
    """
    args = get_args()
    
    # Clean precipiation data
    precip_path = os.path.join(args.cruncep_raw_data_path,
                               args.cruncep_precip_dir)
    precip_files = os.listdir(precip_path)
    ncks_mk_time_rec_dmn.call_ncks(precip_path, precip_files,
                                   os.path.join(args.cruncep_clean_data_path,
                                   args.cruncep_precip_dir))
    
    # Clean solar data
    solar_path = os.path.join(args.cruncep_raw_data_path,
                              args.cruncep_solar_dir)
    solar_files = os.listdir(solar_path)
    ncks_mk_time_rec_dmn.call_ncks(solar_path, solar_files,
                                   os.path.join(args.cruncep_clean_data_path,
                                   args.cruncep_solar_dir))
    
    # Clean temperature, pressure, humidity, wind, and longwave data
    tphwl_path = os.path.join(args.cruncep_raw_data_path,
                              args.cruncep_tphwl_dir)
    tphwl_files = os.listdir(tphwl_path)
    ncks_mk_time_rec_dmn.call_ncks(tphwl_path, tphwl_files,
                                   os.path.join(args.cruncep_clean_data_path,
                                   args.cruncep_tphwl_dir))

def test():
    data = CRUNCEP7()
    
    data.get_precip()
    data.precip_rootgrp.close()
    
    data.get_solar()
    data.solar_rootgrp.close()
    
    data.get_tphwl()
    data.tphwl_rootgrp.close()

def run():
    #clean_data()
    test()

def main():
    run()

if __name__=='__main__':
    main()
