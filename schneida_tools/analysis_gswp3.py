#!/usr/bin/env python

from os import path

import numpy as np
from cartopy import util
import cartopy.feature as cfeature
import cartopy.crs as ccrs
from matplotlib import pyplot as plt
import colorcet as cc
from mpi4py import MPI

from schneida_args import get_args
import gswp3
import wfde5
import coordinate_space
import gris_dem
import ais_dem
import verify_precip
import analysis_era5

import ipdb

def set_map_titles(axes):
    axes[0].set_title('GSWP3')
    axes[1].set_title('WFDE5')
    axes[2].set_title('GSWP3 - WFDE5')
    axes[3].set_title('((GSWP3 - WFDE5)\n/ WFDE5) ' + r'$\times$ 100')

class Analysis(object):
    def __init__(self, compute_means=True, greenland=False, antarctica=False):
        self.compute_means=compute_means
        self.greenland=greenland
        self.antarctica=antarctica
        self.args = get_args()
        
    def compare_temperature(self, cmap='cet_CET_C3_r', degc_min=-7, degc_max=37):
        """
        """
        # Get GSWP3 temporal mean temperature
        gswp3_data = gswp3.GSWP3()
        gswp3_data.get_tphwl()
        gswp3_mean_t_rootgrp = gswp3.get_temporal_mean(gswp3_data.tphwl_rootgrp,
                                                           'TBOT',
                                                           compute=self.compute_means)
        gswp3_data.tphwl_rootgrp.close()
        # Convert from K to degrees C
        gswp3_time_mean_tc = -self.args.TFRZ + gswp3_mean_t_rootgrp.variables['TBOT'][:]
    
        # Get WFDE5 temporal mean temperature
        wfde5_data = wfde5.WFDE5()
        wfde5_data.get_tair()
        wfde5_mean_t_rootgrp = wfde5.get_temporal_mean(wfde5_data.t_air, 'Tair',
                                                       compute=self.compute_means)
        wfde5.close_rootgrps(wfde5_data.t_air)
        # Convert to Celcius and shift WFDE5 data to CRUNCEP grid
        wfde5_time_mean_tc = np.roll(-self.args.TFRZ + wfde5_mean_t_rootgrp.variables['Tair'][:],
                                     360, axis=1)
        # Calculate difference
        print('Computing temperature differences...')
        time_mean_tc_diffs = gswp3_time_mean_tc - wfde5_time_mean_tc
        time_mean_tc_diffs_rel = time_mean_tc_diffs / (wfde5_time_mean_tc + self.args.TFRZ)
    
        # Setup maps
        if True:
            axes = coordinate_space.four_map_horizontal_comparison(greenland=self.greenland,
                                                               antarctica=self.antarctica)
        else:
            ax = analysis_era5.setup_map(greenland=self.greenland,
                                         antarctica=self.antarctica)
        set_map_titles(axes)
        print('Mapping temperature data to figure...')
        # Map data
        longxy = gswp3_mean_t_rootgrp.variables['LONGXY'][:]
        latixy = gswp3_mean_t_rootgrp.variables['LATIXY'][:]
    
        gswp3_quad_mesh = axes[0].pcolormesh(longxy.data, latixy.data,
                                           np.ma.clip(coordinate_space.mask_vals(longxy,
                                                                latixy,
                                                                gswp3_time_mean_tc,
                                                                greenland=self.greenland,
                                                                antarctica=self.antarctica),
                                                      degc_min, degc_max),
                                           shading='nearest', cmap=cmap, 
                                           vmin=degc_min, vmax=degc_max,
                                           transform=ccrs.PlateCarree())
    
        wfde5_quad_mesh = axes[1].pcolormesh(longxy.data, latixy.data,
                                         np.ma.clip(coordinate_space.mask_vals(longxy,
                                                              latixy,
                                                              wfde5_time_mean_tc,
                                                              greenland=self.greenland,
                                                              antarctica=self.antarctica),
                                                    degc_min, degc_max),
                                         shading='nearest', cmap=cmap,
                                         vmin=degc_min, vmax=degc_max,
                                         transform=ccrs.PlateCarree())
    
        diffs_quad_mesh = axes[2].pcolormesh(longxy.data, latixy.data,
                                         np.ma.clip(coordinate_space.mask_vals(longxy,
                                                              latixy,
                                                              time_mean_tc_diffs,
                                                              greenland=self.greenland,
                                                              antarctica=self.antarctica),
                                                    -10, 10),
                                         shading='nearest', cmap='cet_CET_D1',
                                         vmin=-10, vmax=10,
                                         transform=ccrs.PlateCarree())
                    
        rel_diffs_quad_mesh = axes[3].pcolormesh(longxy.data, latixy.data,
                                             np.ma.clip(100. *
                                                        coordinate_space.mask_vals(longxy,
                                                                  latixy,
                                                                  time_mean_tc_diffs_rel,
                                                                  greenland=self.greenland,
                                                                  antarctica=self.antarctica),
                                                        -5, 5),
                                             shading='nearest', cmap='cet_CET_D1',
                                             vmin=-5, vmax=5,
                                             transform=ccrs.PlateCarree())
        # Colorbar
        fig = plt.gcf()
        gswp3_cbar = fig.colorbar(gswp3_quad_mesh,
                                    ax=axes[0:2], orientation='horizontal')
        gswp3_cbar.set_label('Temperature ($^{\circ}$ C)')

        '''
        wfde5_cbar = fig.colorbar(wfde5_quad_mesh,
                                    ax=axes[1], orientation='horizontal')
        wfde5_cbar.set_label('Temperature ($^{\circ}$ C)')
        '''    

        diffs_cbar = fig.colorbar(diffs_quad_mesh,
                                  ax=axes[2], orientation='horizontal')
        diffs_cbar.set_label('Difference (K)')

        rel_cbar = fig.colorbar(rel_diffs_quad_mesh,
                                  ax=axes[3], orientation='horizontal')
        rel_cbar.set_label(r'Difference ($\%$)')
    
        self.draw_elevation_contours(axes)
        
        self.gswp3_mean_t_rootgrp = gswp3_mean_t_rootgrp
        self.wfde5_mean_t_rootgrp = wfde5_mean_t_rootgrp
        
        return(axes)
        
    def close_mean_t_rootgrps(self):
        self.gswp3_mean_t_rootgrp.close()
        self.wfde5_mean_t_rootgrp.close()
    
    def draw_elevation_contours(self, axes):
        """
        Draw contours
        """
        if self.greenland:
            dem = gris_dem.GrISDEM(path.join('data_raw',
                                                'gimpdem_90m_v01.1.tif'))
            #dem.print_dataset_info()
            for i, ax in enumerate(axes):
                # Add elevation contours
                dem.draw_contours(ax,
                              path.join('data_clean', 'gimpdem_90m_v01.1.nc'),
                              downsample=10)
    
        elif self.antarctica:
            # Add elevation contours
            dem = ais_dem.AisDem(path.join('data_raw', 'krigged_dem_nsidc.bin'))
            dem.setup_map(path.join('data_clean', "krigged_dem_nsidc.nc"),
                          path.join('data_clean', "krigged_dem_errormap_nsidc.nc"),
                          new_map=False)
            
            for i, ax in enumerate(axes):
                #Add elevation contours
                dem.draw_contours(ax, path.join('data_clean',
                                               'krigged_dem_nsidc.nc'),
                              path.join('data_clean',
                                           'krigged_dem_errormap_nsidc.nc'))
                
        else:
            for i, ax in enumerate(axes):
                # Add elevation contours
                coordinate_space.draw_elevation_contours(ax)
    
    def compare_precip(self,
                       #cmap='cet_CET_L6_r',
                       cmap='cet_CET_L7_r', cm_per_year_min=0,
                       cm_per_year_max=180):
        """
        """
        # Get GSWP3 temporal mean precipitation
        gswp3_data = gswp3.GSWP3()
        gswp3_data.get_precip()
        gswp3_mean_precip_rootgrp = gswp3.get_temporal_mean(gswp3_data.precip_rootgrp,
                                                                'PRECTmms',
                                                              compute=self.compute_means)
        gswp3_data.precip_rootgrp.close()
        # mm H2O / sec -> cm H2O / yr.
        gswp3_time_mean_precip = (60. * 60. * 24. * 365.25 *
                            gswp3_mean_precip_rootgrp.variables['PRECTmms'][:]) / 10.
    
        # Get WFDE5 temporal mean precipitation
        wfde5_data = wfde5.WFDE5()
        # WFDE5 rainfall
        wfde5_data.get_rainf()
        wfde5_mean_rainf_rootgrp = wfde5.get_temporal_mean(wfde5_data.rainf,
                                                           'Rainf',
                                                           compute=self.compute_means)
        wfde5.close_rootgrps(wfde5_data.rainf)
        # WFDE5 snowfall
        wfde5_data.get_snowf()
        wfde5_mean_snowf_rootgrp = wfde5.get_temporal_mean(wfde5_data.snowf,
                                                           'Snowf',
                                                           compute=self.compute_means)
        wfde5.close_rootgrps(wfde5_data.snowf)
        # Integrate total precip, convert to cm / yr., and shift WFDE5 data
        # toCRUNCEP grid
        wfde5_time_mean_precip = np.roll((60.* 60. * 24. * 365.25 *
                                (wfde5_mean_rainf_rootgrp.variables['Rainf'][:] +
                                 wfde5_mean_snowf_rootgrp.variables['Snowf'][:])) / 10.,
                                 360, axis=1)
        # Calculate difference
        print('Computing precipitation differences...')
        time_mean_precip_diffs = gswp3_time_mean_precip - wfde5_time_mean_precip
        time_mean_precip_diffs_rel = time_mean_precip_diffs / wfde5_time_mean_precip
    
        # Setup maps
        if True:
            axes = coordinate_space.four_map_horizontal_comparison(greenland=self.greenland,
                                                               antarctica=self.antarctica)
            set_map_titles(axes)
        else:
            ax_wfde5 = analysis_era5.setup_map(greenland=self.greenland, antarctica=self.antarctica)
            ax_gs
        # Map data
        print('Mapping precipitation data to figure...')
        longxy = gswp3_mean_precip_rootgrp.variables['LONGXY'][:]
        latixy = gswp3_mean_precip_rootgrp.variables['LATIXY'][:]
        # Add cyclic value to arrays
        longxy=util.add_cyclic_point(longxy)
        for i in range(longxy.shape[0]):
            if True:#longxy[i,-1] < 0:
                # check for negative longitude on right most side,
                # change to positive
                longxy[i, -1] +=360
        latixy=util.add_cyclic_point(latixy)
        era5_time_mean_precip=util.add_cyclic_point(era5_time_mean_precip)
        
        gswp3_quad_mesh = axes[0].contourf(longxy.data, latixy.data,
                                           np.ma.clip(coordinate_space.mask_vals(longxy, latixy,
                                                                gswp3_time_mean_precip,
                                                                greenland=self.greenland,
                                                                antarctica=self.antarctica),
                                                      cm_per_year_min, cm_per_year_max),
                                                levels=int((cm_per_year_max-cm_per_year_min)/5.),
                                                #extend='both',
                                                cmap=cmap,
                                                vmin=cm_per_year_min, vmax=cm_per_year_max,
                                                transform=ccrs.PlateCarree())
                                                
        wfde5_quad_mesh = axes[1].contourf(longxy.data, latixy.data,
                                         np.ma.clip(coordinate_space.mask_vals(longxy, latixy,
                                                              wfde5_time_mean_precip,
                                                              greenland=self.greenland,
                                                              antarctica=self.antarctica),
                                                    cm_per_year_min, cm_per_year_max),
                                levels=int((cm_per_year_max-cm_per_year_min)/5.),
                                #extend='both',
                                cmap=cmap,
                                vmin=cm_per_year_min, vmax=cm_per_year_max,
                                transform=ccrs.PlateCarree())
        '''
        gswp3_quad_mesh = axes[0].pcolormesh(longxy.data, latixy.data,
                                           np.ma.clip(coordinate_space.mask_vals(longxy, latixy,
                                                                gswp3_time_mean_precip,
                                                                greenland=self.greenland,
                                                                antarctica=self.antarctica),
                                                      cm_per_year_min, cm_per_year_max),
                                           shading='nearest', cmap=cmap,
                                           vmin=cm_per_year_min, vmax=cm_per_year_max,
                                           transform=ccrs.PlateCarree())
                                       
        wfde5_quad_mesh = axes[1].pcolormesh(longxy.data, latixy.data,
                                         np.ma.clip(coordinate_space.mask_vals(longxy, latixy,
                                                              wfde5_time_mean_precip,
                                                              greenland=self.greenland,
                                                              antarctica=self.antarctica),
                                                    cm_per_year_min, cm_per_year_max),
                                         shading='nearest', cmap=cmap,
                                         vmin=cm_per_year_min, vmax=cm_per_year_max,
                                         transform=ccrs.PlateCarree())
        '''
        diffs_quad_mesh = axes[2].pcolormesh(longxy.data, latixy.data,
                                         np.ma.clip(coordinate_space.mask_vals(longxy, latixy,
                                                              time_mean_precip_diffs,
                                                              greenland=self.greenland,
                                                              antarctica=self.antarctica),
                                                    -50, 50),
                                         shading='nearest', cmap='cet_CET_D10',
                                         vmin=-50, vmax=50,
                                         transform=ccrs.PlateCarree())
                                     
        rel_diffs_quad_mesh = axes[3].pcolormesh(longxy.data, latixy.data,
                                             np.ma.clip(100. *
                                                        coordinate_space.mask_vals(longxy, latixy,
                                                                  time_mean_precip_diffs_rel,
                                                                  greenland=self.greenland,
                                                                  antarctica=self.antarctica),
                                                        -50, 50),
                                             shading='nearest', cmap='cet_CET_D10',
                                             vmin=-50, vmax=50,
                                             transform=ccrs.PlateCarree())
        # Colorbar
        fig = plt.gcf()
        gswp3_cbar = fig.colorbar(gswp3_quad_mesh,
                            ax=axes[0:2], orientation='vertical')
        gswp3_cbar.set_label('precipitation rate (cm w.eq. yr$^{-1}$)')
        '''
        wfde5_cbar = fig.colorbar(wfde5_quad_mesh,
                            ax=axes[1], orientation='horizontal')
        wfde5_cbar.set_label('Precipitation (cm H$_2$O yr$^{-1}$)')
        '''
        diffs_cbar = fig.colorbar(diffs_quad_mesh,
                            ax=axes[2], orientation='vertical')
        diffs_cbar.set_label('Difference (cm w.eq. yr$^{-1}$)')

        rel_cbar = fig.colorbar(rel_diffs_quad_mesh,
                            ax=axes[3], orientation='vertical')
        rel_cbar.set_label(r'Difference ($\%$)')

        self.draw_elevation_contours(axes)

        self.wfde5_mean_rainf_rootgrp = wfde5_mean_rainf_rootgrp
        self.wfde5_mean_snowf_rootgrp = wfde5_mean_snowf_rootgrp
        self.gswp3_mean_precip_rootgrp = gswp3_mean_precip_rootgrp
        self.precip_cmap = cmap
        self.precip_cm_per_year_min = cm_per_year_min
        self.precip_cm_per_year_max = cm_per_year_max
        
        return axes
        
    def close_mean_precip_rootgrps(self):
        self.wfde5_mean_rainf_rootgrp.close()
        self.wfde5_mean_snowf_rootgrp.close()
        self.gswp3_mean_precip_rootgrp.close()
        
def test():
    gswp3.test()
    wfde5.test()

def run(debug=False):
    '''
    plt.style.use(path.join('schneida_tools', 'gmd_movie_frame.mplstyle'))
    plt.style.use('uci_darkblue')
    '''
    plt.style.use('agu_online_poster_presentation')
    #plt.style.use('agu_full')
    #plt.style.use('hofmann_talk')
    
    if debug:
        rank=0
    else:
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
    
    greenland_analysis = Analysis(compute_means=False,
                                  greenland=True)
    antarctica_analysis = Analysis(compute_means=False,
                                   antarctica=True)
    northern_hemisphere_analysis = Analysis(compute_means=False)
    
    # Greenland analysis
    if rank==0:
        # Temperature
        axes = greenland_analysis.compare_temperature(degc_min=-50, degc_max=0,
                                               #cmap='cet_CET_L7',
                                                cmap='cet_CET_CBL3_r')
        # Set the figure title
        plt.suptitle('Greenland mean 1980 to 1990 surface air temperature')
        savefig_name = path.join('results', 'greenland_tair_gswp3_vs_wfde5.png')
        print('Writing temperature maps to %s' % savefig_name)
        plt.savefig(savefig_name, dpi=600)
        # Close figure and files
        plt.close()
        greenland_analysis.close_mean_t_rootgrps()

    if rank==1:
        # Precipitation
        (sumup_gris, sumup_ais) = verify_precip.grid_sumup2wfde5()
        comm.send(sumup_ais, dest=3)
        axes = greenland_analysis.compare_precip(cm_per_year_min=0, cm_per_year_max=150)
        # Get and plot SUMup locations
        axes[0].scatter(sumup_gris[1], sumup_gris[0], s=sumup_gris[3], c=sumup_gris[2],
                        cmap=greenland_analysis.precip_cmap, vmin=greenland_analysis.precip_cm_per_year_min,
                        vmax=greenland_analysis.precip_cm_per_year_max, edgecolors='black',
                        linewidths=0.5,
                        transform=ccrs.PlateCarree())
        axes[1].scatter(sumup_gris[1], sumup_gris[0], s=sumup_gris[3], c=sumup_gris[2],
                        cmap=greenland_analysis.precip_cmap, vmin=greenland_analysis.precip_cm_per_year_min,
                        vmax=greenland_analysis.precip_cm_per_year_max, edgecolors='black',
                        linewidths=0.5,
                        transform=ccrs.PlateCarree())
    
        # Set the figure title
        plt.suptitle('Greenland mean 1980 to 1990 precipitation')
        savefig_name = path.join('results', 'greenland_precip_gswp3_vs_wfde5.png')
        print('Writing precipitation maps to %s' % savefig_name)
        plt.savefig(savefig_name, dpi=600)
        # Close figure and files
        plt.close()
        greenland_analysis.close_mean_precip_rootgrps()
    
    # Antarctica analysis
    if rank==2:
        # Temperature
        axes = antarctica_analysis.compare_temperature(degc_min=-50, degc_max=0,
                                                       #cmap='cet_CET_L7'
                                                       cmap='cet_CET_CBL2'
                                                      )
        # Set the figure title
        plt.suptitle('Antarctica mean 1980 to 1990 surface air temperature')
        savefig_name = path.join('results', 'antarctica_tair_gswp3_vs_wfde5.png')
        print('Writing temperature maps to %s' % savefig_name)
        plt.savefig(savefig_name, dpi=600)
        # Close figure and files
        plt.close()
        antarctica_analysis.close_mean_t_rootgrps()
    
    if rank==3:
        # Precipitation
        sumup_ais = comm.recv(source=1)
        axes = antarctica_analysis.compare_precip(cm_per_year_min=0, cm_per_year_max=150)
        # Get and plot SUMup locations
        axes[0].scatter(sumup_ais[1], sumup_ais[0], s=sumup_ais[3], c=sumup_ais[2],
                        cmap=antarctica_analysis.precip_cmap, vmin=antarctica_analysis.precip_cm_per_year_min,
                        vmax=antarctica_analysis.precip_cm_per_year_max, edgecolors='black',
                        linewidths=0.5,
                        transform=ccrs.PlateCarree())
        axes[1].scatter(sumup_ais[1], sumup_ais[0], s=sumup_ais[3], c=sumup_ais[2],
                        cmap=antarctica_analysis.precip_cmap, vmin=antarctica_analysis.precip_cm_per_year_min,
                        vmax=antarctica_analysis.precip_cm_per_year_max, edgecolors='black',
                        linewidths=0.5,
                        transform=ccrs.PlateCarree())
    
        # Set the figure title
        plt.suptitle('Antarctica mean 1980 to 1990 precipitation')
        savefig_name = path.join('results', 'antarctica_precip_gswp3_vs_wfde5.png')
        plt.savefig(savefig_name, dpi=600)
        # Close figure and files
        plt.close()
        antarctica_analysis.close_mean_precip_rootgrps()
    
    if rank==4:
        # Northern hemisphere analysis
        axes = northern_hemisphere_analysis.compare_temperature()
        # Set the figure title
        plt.suptitle('Northern Hemisphere mean 1980 to 1990 surface air temperature')
        savefig_name = path.join('results', 'nh_tair_gswp3_vs_wfde5.png')
        print('Writing temperature maps to %s' % savefig_name)
        plt.savefig(savefig_name, dpi=600)
        # Close figure and files
        plt.close()
        northern_hemisphere_analysis.close_mean_t_rootgrps()
    
    if rank==5:
        # Precipitation
        axes = northern_hemisphere_analysis.compare_precip()
        # Set the figure title
        plt.suptitle('Northern Hemisphere mean 1980 to 1990 precipitation')
        savefig_name = path.join('results', 'nh_precip_gswp3_vs_wfde5.png')
        plt.savefig(savefig_name, dpi=600)
        # Close figure and files
        plt.close()
        northern_hemisphere_analysis.close_mean_precip_rootgrps()
    
def main():
    run()

if __name__=='__main__':
    main()
