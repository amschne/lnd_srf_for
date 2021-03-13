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
    
    def get_lwdown(self):
        """ Get downwelling longwave radiation
        """
        self.lw_down = self.get_data('LWdown_WFDE5_CRU_*_v1.1.nc')
        
    def get_psurf(self):
        """ Get surface pressure
        """
        self.p_surf = self.get_data('PSurf_WFDE5_CRU_*_v1.1.nc')
    
    def get_qair(self):
        """ Get humidity
        """
        self.q_air = self.get_data('Qair_WFDE5_CRU_*_v1.1.nc')
        
    def get_rainf(self):
        """ Get rainfall
        """
        self.rainf = self.get_data('Rainf_WFDE5_CRU+GPCC_*_v1.1.nc')
        
    def get_snowf(self):
        """ Get snowfall
        """
        self.snowf = self.get_data('Snowf_WFDE5_CRU+GPCC_*_v1.1.nc')
    
    def get_swdown(self):
        """ Get downwelling shortwave radiation
        """
        self.sw_down = self.get_data('SWdown_WFDE5_CRU_*_v1.1.nc')
        
    def get_tair(self):
        """ Get air temperature
        """
        self.t_air = self.get_data('Tair_WFDE5_CRU_*_v1.1.nc')
    
    def get_wind(self):
        """ Get wind
        """
        self.wind = self.get_data('Wind_WFDE5_CRU_*_v1.1.nc')

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

def run():
    clean_data()
    
def main():
    run()

if __name__=='__main__':
    main()
