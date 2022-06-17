#!/usr/bin/env python

"""
"""

from os import path

import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib import pyplot as plt
import colorcet as cc

from schneida_tools.schneida_args import get_args
from schneida_tools import merra2
from schneida_tools import analysis_era5
from schneida_tools import coordinate_space
from schneida_tools import gris_dem
from schneida_tools import ais_dem
from schneida_tools import verify_precip_era5 as verify_precip

SURF_AIR_TEMP = 'TLML'
TOT_PREC = 'PRECTOT'

class Analysis(object):
    def __init__(self, compute_means=True, greenland=False, antarctica=False):
        self.compute_means=compute_means
        self.greenland=greenland
        self.antarctica=antarctica
        self.args = get_args()

    def compare_temperature(self, cmap='cet_CET_L8', degc_min=-7, degc_max=37):
        # Get MERRA-2 temporal mean temperature
        merra2_data = merra2.MERRA2()
        merra2_data.get_tair()
        merra2_mean_t_rootgrp = merra2.get_temporal_mean(merra2_data.t_air, SURF_AIR_TEMP,
                                                       compute=self.compute_means)
        merra2.close_rootgrps(merra2_data.t_air)

        # Convert to Celcius
        merra2_time_mean_tc = -self.args.TFRZ + merra2_mean_t_rootgrp.variables[SURF_AIR_TEMP][:]
        
        # Setup maps
        ax = analysis_era5.setup_map(greenland=self.greenland,
                                     antarctica=self.antarctica)
            
        #set_map_titles(axes)
        print('Mapping temperature data to figure...')
        
        # Map data
        longxy, latixy = np.meshgrid(merra2_mean_t_rootgrp.variables['lon'][:],
                                     merra2_mean_t_rootgrp.variables['lat'][:])

        merra2_quad_mesh = ax.pcolormesh(longxy, latixy,
                                      np.ma.clip(coordinate_space.mask_vals(longxy,
                                                           latixy,
                                                            merra2_time_mean_tc,
                                                            greenland=self.greenland,
                                                            antarctica=self.antarctica),
                                                    degc_min, degc_max),
                                         shading='nearest', cmap=cmap,
                                         vmin=degc_min, vmax=degc_max,
                                         transform=ccrs.PlateCarree())

        fig = plt.gcf()
        merra2_cbar = fig.colorbar(merra2_quad_mesh,
                                    ax=ax, orientation='vertical')
        merra2_cbar.set_label('Temperature ($^{\circ}$ C)')

        self.draw_elevation_contours(ax)
        self.merra2_mean_t_rootgrp = merra2_mean_t_rootgrp
        
        return(ax)

    def close_mean_t_rootgrps(self):
        self.merra2_mean_t_rootgrp.close()

    def draw_elevation_contours(self, ax):
        """
        Draw contours
        """
        if self.greenland:
            dem = gris_dem.GrISDEM(path.join('data_raw',
                                                'gimpdem_90m_v01.1.tif'))
            #dem.print_dataset_info()
            dem.draw_contours(ax,
                              path.join('data_clean', 'gimpdem_90m_v01.1.nc'),
                              downsample=10)
    
        elif self.antarctica:
            # Add elevation contours
            dem = ais_dem.AisDem(path.join('data_raw', 'krigged_dem_nsidc.bin'))
            dem.setup_map(path.join('data_clean', "krigged_dem_nsidc.nc"),
                          path.join('data_clean', "krigged_dem_errormap_nsidc.nc"),
                          new_map=False)
            
            dem.draw_contours(ax, path.join('data_clean',
                                               'krigged_dem_nsidc.nc'),
                              path.join('data_clean',
                                           'krigged_dem_errormap_nsidc.nc'))
                
        else:
            coordinate_space.draw_elevation_contours(ax)

def run():
    # Greenland analysis
    greenland_analysis = Analysis(compute_means=False,
                                  greenland=True)
    # Temperature
    ax = greenland_analysis.compare_temperature(degc_min=-50, degc_max=0,
                                           #cmap='cet_CET_L7',
                                           #cmap='cet_CET_CBL2'
                                           cmap='cet_CET_L17')
    # Set the figure title
    plt.suptitle('Greenland mean 1980 to 1990 surface air temperature')
    savefig_name = path.join('results', 'greenland_tair_merra2.png')
    print('Writing temperature maps to %s' % savefig_name)
    plt.savefig(savefig_name, dpi=300)
    # Close figure and files
    plt.close()
    greenland_analysis.close_mean_t_rootgrps()
    
    # Antarctica analysis
    antarctica_analysis = Analysis(compute_means=False,
                                   antarctica=True)
    # Temperature
    ax = antarctica_analysis.compare_temperature(degc_min=-50, degc_max=0,
                                                   #cmap='cet_CET_L7'
                                                   #cmap='cet_CET_CBL2'
                                                   cmap='cet_CET_L17'
                                                  )
    # Set the figure title
    plt.suptitle('Antarctica mean 1980 to 1990 surface air temperature')
    savefig_name = path.join('results', 'antarctica_tair_merra2.png')
    print('Writing temperature maps to %s' % savefig_name)
    plt.savefig(savefig_name, dpi=300)
    # Close figure and files
    plt.close()
    antarctica_analysis.close_mean_t_rootgrps()
    
    # Northern hemisphere analysis
    northern_hemisphere_analysis = Analysis(compute_means=False)
    ax = northern_hemisphere_analysis.compare_temperature()
    # Set the figure title
    plt.suptitle('Northern Hemisphere mean 1980 to 1990 surface air temperature')
    savefig_name = path.join('results', 'nh_tair_merra2.png')
    print('Writing temperature maps to %s' % savefig_name)
    plt.savefig(savefig_name, dpi=300)
    # Close figure and files
    plt.close()
    northern_hemisphere_analysis.close_mean_t_rootgrps()

def main():
    run()

if __name__=='__main__':
    main()