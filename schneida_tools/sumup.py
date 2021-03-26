#!/usr/bin/env python

""" Read data from the Surface Mass Balance and Snow on Sea Ice Working Group
    (SUMup) dataset (Montgomery et al., 2018).

    References

    Montgomery, L., Koenig, L., and Alexander, P.: The SUMup dataset: compiled
    measurements of surface mass balance components over ice sheets and sea ice
    with analysis over Greenland, Earth Syst. Sci. Data, 10, 1959–1985,
    https://doi.org/10.5194/essd-10-1959-2018, 2018.
"""

from os import path
from netCDF4 import Dataset

import ipdb

class SUMupAccumulation(object):
    def __init__(self, file_path=path.join('data_raw',
                     'SUMup_dataset_july2018_accumulation_on_land_ice.nc')):
        self.rootgrp = Dataset(file_path)

def get_accumulation():
    """
    """
    sumup_accum = SUMupAccumulation()
    return sumup_accum.rootgrp
    
def run():
    rootgrp = get_accumulation()
    ipdb.set_trace()

def main():
    run()

if __name__=='__main__':
    main()