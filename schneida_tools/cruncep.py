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

from schneida_tools import ncks_mk_time_rec_dmn

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cruncep_raw_data_path',
                        default=os.path.join('data_raw',
                          'atm_forcing.datm7.cruncep_qianFill.0.5d.v7.c160715'),
                        help='Path to top level directory where raw CRUNCEP '
                             'data is stored.')
    parser.add_argument('--cruncep_clean_data_path',
                        default=os.path.join('data_clean',
                        'atm_forcing.datm7.cruncep_qianFill.0.5d.v7.c160715'),
                        help='Path to top level directory where clean CRUNCEP '
                             'data is stored.')
    parser.add_argument('--cruncep_precip_dir', default='Precip6Hrly')
    parser.add_argument('--cruncep_solar_dir', default='Solar6Hrly')
    parser.add_argument('--cruncep_tphwl_dir', default='TPHWL6Hrly')
    
    args = parser.parse_args()
    
    return args
    
def call_ncks(input_path, input_files, output_path):
    for i, input_file in enuemrate(input_files):
        ncks_mk_time_rec_dmn.run(os.path.abspath(os.path.join(input_path,
                                                              input_file)),
                                 os.path.abspath(os.path.join(output_path,
                                                              input_file)))
def clean_data():
    """ 1. Gather raw data
        2. Call ncks to change the time dimension to the record dimension.
    """
    args = get_args()
    
    # Clean precipiation data
    precip_path = os.path.join(args.cruncep_raw_data_path,
                               args.cruncep_precip_dir)
    precip_files = os.listdir(precip_path)
    call_ncks(precip_path, precip_files,
              os.path.join(args.cruncep_clean_data_path,
                           args.cruncep_precip_dir))
    
    # Clean solar data
    solar_path = os.path.join(args.cruncep_raw_data_path,
                              args.cruncep_solar_dir)
    solar_files = os.listdir(solar_path)
    call_ncks(solar_path, solar_files,
              os.path.join(args.cruncep_clean_data_path,
                           args.cruncep_solar_dir))
    
    # Clean temperature, pressure, humidity, wind, and longwave data
    tphwl_path = os.path.join(args.cruncep_raw_data_path,
                              args.cruncep_tphwl_dir)
    tphwl_files = os.listdir(tphwl_path)
    call_ncks(tphwl_path, tphwl_files,
              os.path.join(args.cruncep_clean_data_path,
                           args.cruncep_tphwl_dir))
    
def run():
    clean_data()

def main():
    run()

if __name__=='__main__':
    main()