import argparse
from os import path

TFRZ = 273.15 # K
GRIS_EXTENT = (-361, 361, 45, 90)
NH_EXTENT = (-361, 361, -30, 90)
AIS_EXTENT = (-361, 361, -90, -45)

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
    
    # GSWP3 data info
    parser.add_argument('--gswp3_raw_data_path',
                        default=path.join('data_raw',
                          'atm_forcing.datm7.GSWP3.0.5d.v1.c170516'),
                        help='Path to top level directory where raw GSWP3 '
                             'data is stored.')
    parser.add_argument('--gswp3_clean_data_path',
                        default=path.join('data_clean',
                        'atm_forcing.datm7.GSWP3.0.5d.v1.c170516'),
                        help='Path to top level directory where clean GSWP3 '
                             'data is stored.')
    parser.add_argument('--gswp3_precip_dir', default='Precip')
    parser.add_argument('--gswp3_solar_dir', default='Solar')
    parser.add_argument('--gswp3_tphwl_dir', default='TPHWL')
    
    # WFDE5 data info
    parser.add_argument('--wfde5_raw_data_path',
                        default=path.join('data_raw', 'wfde5'),
                        help='Path to top level directory where raw WFDE5 '
                             'data is stored.')
    parser.add_argument('--wfde5_clean_data_path',
                        default=path.join('data_clean', 'wfde5'),
                        help='Path to top level directory where clean WFDE5 '
                             'data is stored.')
                             
    # ERA5 data info
    parser.add_argument('--era5_raw_data_path',
                        default=path.join('data_raw', 'era5'),
                        help='Path to top level directory where raw ERA5 '
                             'data is stored.')
    parser.add_argument('--era5_clean_data_path',
                        default=path.join('data_clean', 'era5'),
                        help='Path to top level directory where clean ERA5 '
                             'data is stored.')
                             
    # MERRA2 data info
    parser.add_argument('--merra2_raw_data_path',
                        default=path.join('data_raw', 'merra2',
                                          'MERRA2_MONTHLY', 'M2TMNXFLX.5.12.4'),
                        help='Path to top level directory where raw MERRA-2 '
                             'data is stored.')
    parser.add_argument('--merra2_clean_data_path',
                        default=path.join('data_clean', 'merra2'),
                        help='Path to top level directory where clean MERRA-2 '
                             'data is stored.')
    
    # SUMup analysis parameters
    parser.add_argument('--sumup_start_year', type=int,
                        default=1980,
                        help='First year to include in SUMup analysis')
    parser.add_argument('--sumup_stop_year', type=int,
                        default=1990,
                        help='First year to exclude in SUMup analysis')
    
    args = parser.parse_args()
    
    
    # Add constants to args
    args.TFRZ = TFRZ
    args.GRIS_EXTENT = GRIS_EXTENT
    args.NH_EXTENT = NH_EXTENT
    args.AIS_EXTENT = AIS_EXTENT
    
    return args
