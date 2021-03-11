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

from schneida_tools import ncks_mk_time_rec_dmn
from schneida_tools.schneida_args import get_args

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
