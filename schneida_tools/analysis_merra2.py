#!/usr/bin/env python

"""
"""

from os import path

import numpy as np
from cartopy import util
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib import pyplot as plt
import colorcet as cc
from mpi4py import MPI

'''
from schneida_tools.schneida_args import get_args
from schneida_tools import merra2
from schneida_tools import analysis_era5
from schneida_tools import coordinate_space
from schneida_tools import gris_dem
from schneida_tools import ais_dem
from schneida_tools import verify_precip_era5 as verify_precip
'''

from schneida_args import get_args
import merra2
import analysis_era5
import coordinate_space
import gris_dem
import ais_dem
import verify_precip_merra2 as verify_precip

import ipdb

class Analysis(object):
    def __init__(self, compute_means=True, greenland=False, antarctica=False):
        self.compute_means=compute_means
        self.greenland=greenland
        self.antarctica=antarctica
        self.args = get_args()

    def compare_temperature(self, cmap='cet_CET_C3_r', degc_min=-7, degc_max=37):
        # Get MERRA-2 temporal mean temperature
        merra2_data = merra2.MERRA2()
        merra2_data.get_tair()
        merra2_mean_t_rootgrp = merra2.get_temporal_mean(merra2_data.t_air, merra2.SURF_AIR_TEMP,
                                                       compute=self.compute_means)
        merra2.close_rootgrps(merra2_data.t_air)

        # Convert to Celcius
        merra2_time_mean_tc = -self.args.TFRZ + merra2_mean_t_rootgrp.variables[merra2.SURF_AIR_TEMP][:]
        
        # Setup maps
        ax = analysis_era5.setup_map(greenland=self.greenland,
                                     antarctica=self.antarctica)
            
        #set_map_titles(axes)
        print('Mapping temperature data to figure...')
        
        # Map data
        longxy, latixy = np.meshgrid(merra2_mean_t_rootgrp.variables['lon'][:],
                                     merra2_mean_t_rootgrp.variables['lat'][:])
        if False:
            merra2_quad_mesh = ax.pcolormesh(longxy, latixy,
                                      np.ma.clip(coordinate_space.mask_vals(longxy,
                                                           latixy,
                                                            merra2_time_mean_tc,
                                                            greenland=self.greenland,
                                                            antarctica=self.antarctica),
                                                    degc_min, degc_max),
                                         shading='nearest',
                                         #shading='gouraud',
                                         cmap=cmap,
                                         vmin=degc_min, vmax=degc_max,
                                         transform=ccrs.PlateCarree())
        else:
            # Add cyclic value to arrays
            longxy=util.add_cyclic_point(longxy)
            for i in range(longxy.shape[0]):
                if True:#longxy[i,-1] < 0:
                    # check for negative longitude on right most side,
                    # change to positive
                    longxy[i, -1] +=360
                
            latixy=util.add_cyclic_point(latixy)
            merra2_time_mean_tc=util.add_cyclic_point(merra2_time_mean_tc)
            merra2_quad_mesh = ax.contourf(longxy, latixy,
                                        np.ma.clip(coordinate_space.mask_vals(longxy,
                                                             latixy,
                                                              merra2_time_mean_tc,
                                                              greenland=self.greenland,
                                                              antarctica=self.antarctica),
                                                      degc_min, degc_max),
                                        levels=degc_max-degc_min,
                                        cmap=cmap,
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
    def compare_precip(self,
                       #cmap='cet_CET_L6_r',
                       cmap='cet_CET_L7_r', cm_per_year_min=0,
                       cm_per_year_max=180, ax=None, elevation_levels=500):
        """
        """
        # Get MERRA-2 temporal mean precipitation
        merra2_data = merra2.MERRA2()
        merra2_data.get_precip()
        merra2_mean_precip_rootgrp = merra2.get_temporal_mean(merra2_data.precip,
                                                          merra2.TOT_PREC,
                                                          compute=self.compute_means)
        merra2.close_rootgrps(merra2_data.precip)
        # Integrate mean precip, convert from mm / s to cm / yr.
        merra2_time_mean_precip = (60.*60.*24.*365.25 *
                    merra2_mean_precip_rootgrp.variables[merra2.TOT_PREC][:]) / 10.
        # Setup maps
        if ax is None:
            ax = analysis_era5.setup_map(greenland=self.greenland,
                                         antarctica=self.antarctica)
        # Map data
        print('Mapping precipitation data to figure...')
        longxy, latixy = np.meshgrid(merra2_mean_precip_rootgrp.variables['lon'][:],
                                     merra2_mean_precip_rootgrp.variables['lat'][:])
        if False:
            merra2_quad_mesh = ax.pcolormesh(longxy, latixy,
                                         np.ma.clip(coordinate_space.mask_vals(longxy, latixy,
                                                              merra2_time_mean_precip,
                                                              greenland=self.greenland,
                                                              antarctica=self.antarctica),
                                                    cm_per_year_min, cm_per_year_max),
                                         shading='nearest', cmap=cmap,
                                         vmin=cm_per_year_min, vmax=cm_per_year_max,
                                         transform=ccrs.PlateCarree())
        else:
            # Add cyclic value to arrays
            longxy=util.add_cyclic_point(longxy)
            for i in range(longxy.shape[0]):
                if True:#longxy[i,-1] < 0:
                    # check for negative longitude on right most side,
                    # change to positive
                    longxy[i, -1] +=360
                
            latixy=util.add_cyclic_point(latixy)
            merra2_time_mean_precip=util.add_cyclic_point(merra2_time_mean_precip)
            merra2_quad_mesh = ax.contourf(longxy, latixy,
                                     np.ma.clip(coordinate_space.mask_vals(
                                                          longxy,
                                                          latixy,
                                                          merra2_time_mean_precip,
                                                          greenland=self.greenland,
                                                          antarctica=self.antarctica),
                                                cm_per_year_min, cm_per_year_max),
                                        #extend='both',
                                        levels=int((cm_per_year_max-cm_per_year_min)/5),
                                        cmap=cmap,
                                        vmin=cm_per_year_min, vmax=cm_per_year_max,
                                        transform=ccrs.PlateCarree())

        '''
        # Colorbar
        fig = plt.gcf()
        merra2_cbar = fig.colorbar(merra2_quad_mesh,
                            ax=ax, orientation='vertical')
        merra2_cbar.set_label('precipitation (cm w.eq. yr$^{-1}$)')
        '''
        self.draw_elevation_contours(ax, levels_interval=elevation_levels)

        self.merra2_mean_precip_rootgrp = merra2_mean_precip_rootgrp
        
        self.precip_cmap = cmap
        self.precip_cm_per_year_min = cm_per_year_min
        self.precip_cm_per_year_max = cm_per_year_max
        
        ax.set_title(' MERRA-2')
        return ax
        
    def close_mean_precip_rootgrps(self):
        '''
        '''
        self.merra2_mean_precip_rootgrp.close()

    def draw_elevation_contours(self, ax, levels_interval=500):
        """
        Draw contours
        """
        if self.greenland:
            dem = gris_dem.GrISDEM(path.join('data_raw',
                                                'gimpdem_90m_v01.1.tif'))
            #dem.print_dataset_info()
            dem.draw_contours(ax,
                              path.join('data_clean', 'gimpdem_90m_v01.1.nc'),
                              downsample=10, levels_interval=levels_interval)
    
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

def run(debug=False):
    if debug:
        rank=0
    else:
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
    # Greenland analysis
    greenland_analysis = Analysis(compute_means=False,
                                      greenland=True)
    antarctica_analysis = Analysis(compute_means=False,
                                       antarctica=True)
    northern_hemisphere_analysis = Analysis(compute_means=False)
    if rank ==0:
        # Temperature
        ax = greenland_analysis.compare_temperature(degc_min=-50, degc_max=0,
                                               #cmap='cet_CET_L7',
                                               cmap='cet_CET_CBL3_r'
                                               #cmap='cet_CET_L12'
                                               )
        # Set the figure title
        plt.suptitle('Greenland mean 1980 to 1990 surface air temperature')
        savefig_name = path.join('results', 'greenland_tair_merra2.png')
        print('Writing temperature maps to %s' % savefig_name)
        plt.savefig(savefig_name, dpi=600)
        # Close figure and files
        plt.close()
        greenland_analysis.close_mean_t_rootgrps()
    
    if rank==1:
        # Precipitation
        (sumup_gris, sumup_ais) = verify_precip.grid_sumup2merra2()
        comm.send(sumup_ais, dest=3)
        ax = greenland_analysis.compare_precip(cm_per_year_min=0., cm_per_year_max=150.)
        ax.scatter(sumup_gris[1], sumup_gris[0], s=sumup_gris[3], c=sumup_gris[2],
                        cmap=greenland_analysis.precip_cmap, vmin=greenland_analysis.precip_cm_per_year_min,
                        vmax=greenland_analysis.precip_cm_per_year_max, edgecolors='black',
                        linewidths=0.5,
                        transform=ccrs.PlateCarree())
        # Set the figure title
        plt.suptitle('Greenland mean 1980 to 1990 precipitation')
        savefig_name = path.join('results', 'greenland_precip_merra2.png')
        print('Writing precipitation maps to %s' % savefig_name)
        plt.savefig(savefig_name, dpi=600)
        # Close figure and files
        plt.close()
        greenland_analysis.close_mean_precip_rootgrps()

    if rank==2:
        # Temperature
        ax = antarctica_analysis.compare_temperature(degc_min=-50, degc_max=0,
                                                       #cmap='cet_CET_L7'
                                                       cmap='cet_CET_CBL3_r'
                                                       #cmap='cet_CET_L12'
                                                      )
        # Set the figure title
        plt.suptitle('Antarctica mean 1980 to 1990 surface air temperature')
        savefig_name = path.join('results', 'antarctica_tair_merra2.png')
        print('Writing temperature maps to %s' % savefig_name)
        plt.savefig(savefig_name, dpi=600)
        # Close figure and files
        plt.close()
        antarctica_analysis.close_mean_t_rootgrps()
    
    if rank==3:
        # Precipitation
        sumup_ais = comm.recv(source=1)
        ax = antarctica_analysis.compare_precip(cm_per_year_min=0, cm_per_year_max=150)
        ax.scatter(sumup_ais[1], sumup_ais[0], s=sumup_ais[3], c=sumup_ais[2],
                        cmap=antarctica_analysis.precip_cmap, vmin=antarctica_analysis.precip_cm_per_year_min,
                        vmax=antarctica_analysis.precip_cm_per_year_max, edgecolors='black',
                        linewidths=0.5,
                        transform=ccrs.PlateCarree())
        
         # Set the figure title
        plt.suptitle('Antarctica mean 1980 to 1990 precipitation')
        savefig_name = path.join('results', 'antarctica_precip_merra2.png')
        print('Writing precipitation maps to %s' % savefig_name)
        plt.savefig(savefig_name, dpi=600)
        # Close figure and files
        plt.close()
        antarctica_analysis.close_mean_precip_rootgrps()

    if rank==4:
        # Northern hemisphere analysis
        ax = northern_hemisphere_analysis.compare_temperature()
        # Set the figure title
        plt.suptitle('Northern Hemisphere mean 1980 to 1990 surface air temperature')
        savefig_name = path.join('results', 'nh_tair_merra2.png')
        print('Writing temperature maps to %s' % savefig_name)
        plt.savefig(savefig_name, dpi=600)
        # Close figure and files
        plt.close()
        northern_hemisphere_analysis.close_mean_t_rootgrps()
        
    if rank==5:
        # Precipitation
        ax = northern_hemisphere_analysis.compare_precip(cmap='cet_CET_L1_r')
        # Set the figure title
        plt.suptitle('Northern Hemisphere mean 1980 to 1990 precipitation')
        savefig_name = path.join('results', 'nh_precip_merra2.png')
        print('Writing precipitation maps to %s' % savefig_name)
        plt.savefig(savefig_name, dpi=600)
        # Close figure and files
        plt.close()
        northern_hemisphere_analysis.close_mean_precip_rootgrps()
        
def main():
    run()

if __name__=='__main__':
    main()
