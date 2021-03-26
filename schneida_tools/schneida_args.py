import argparse
from os import path

def get_args():
    parser = argparse.ArgumentParser()
    
    # CRUNCEP data info
    parser.add_argument('--cruncep_raw_data_path',
                        default=path.join('data_raw',
                          'atm_forcing.datm7.cruncep_qianFill.0.5d.v7.c160715'),
                        help='Path to top level directory where raw CRUNCEP '
                             'data is stored.')
    parser.add_argument('--cruncep_clean_data_path',
                        default=path.join('data_clean',
                        'atm_forcing.datm7.cruncep_qianFill.0.5d.v7.c160715'),
                        help='Path to top level directory where clean CRUNCEP '
                             'data is stored.')
    parser.add_argument('--cruncep_precip_dir', default='Precip6Hrly')
    parser.add_argument('--cruncep_solar_dir', default='Solar6Hrly')
    parser.add_argument('--cruncep_tphwl_dir', default='TPHWL6Hrly')
    
    # WFDE5 data info
    parser.add_argument('--wfde5_raw_data_path',
                        default=path.join('data_raw', 'wfde5'),
                        help='Path to top level directory where raw WFDE5 '
                             'data is stored.')
    parser.add_argument('--wfde5_clean_data_path',
                        default=path.join('data_clean', 'wfde5'),
                        help='Path to top level directory where clean WFDE5 '
                             'data is stored.')
    
    # SUMup analysis parameters
    parser.add_argument('--sumup_start_year', type=int,
                        default=1980,
                        help='First year to include in SUMup analysis')
    parser.add_argument('--sumup_stop_year', type=int,
                        default=1990,
                        help='First year to exclude in SUMup analysis')
    
    args = parser.parse_args()
    
    return args
