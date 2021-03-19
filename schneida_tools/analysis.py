#!/usr/bin/env python

from os import path

import numpy as np
import cartopy.crs as ccrs
from matplotlib import pyplot as plt
import colorcet as cc

from schneida_tools import cruncep
from schneida_tools import wfde5
from schneida_tools import coordinate_space
from schneida_tools import gris_dem

import ipdb

TFRZ = 273.15 # K

def set_map_titles(axes):
    axes[0].set_title('CRUNCEP7')
    axes[1].set_title('WFDE5')
    axes[2].set_title('CRUNCEP7 - WFDE5')
    axes[3].set_title('((CRUNCEP7 - WFDE5)\n/ WFDE5) ' + r'$\times$ 100')

def compare_temperature(compute_means=True, cmap='cet_CET_L8', vmin=-7, vmax=37,
                        greenland=False):
    """
    """
    # Get CRUNCEP temporal mean temperature
    cruncep_data = cruncep.CRUNCEP7()
    cruncep_data.get_tphwl()
    cruncep_mean_t_rootgrp = cruncep.get_temporal_mean(cruncep_data.tphwl_rootgrp,
                                                       'TBOT',
                                                       compute=compute_means)
    cruncep_data.tphwl_rootgrp.close()
    # Convert from K to degrees C
    cruncep_time_mean_tc = -TFRZ + cruncep_mean_t_rootgrp.variables['TBOT'][:]
    
    # Get WFDE5 temporal mean temperature
    wfde5_data = wfde5.WFDE5()
    wfde5_data.get_tair()
    wfde5_mean_t_rootgrp = wfde5.get_temporal_mean(wfde5_data.t_air, 'Tair',
                                                   compute=compute_means)
    wfde5.close_rootgrps(wfde5_data.t_air)
    # Convert to Celcius and shift WFDE5 data to CRUNCEP grid
    wfde5_time_mean_tc = np.roll(-TFRZ + wfde5_mean_t_rootgrp.variables['Tair'][:],
                                 360, axis=1)
    # Calculate difference
    print('Computing temperature differences...')
    time_mean_tc_diffs = cruncep_time_mean_tc - wfde5_time_mean_tc
    time_mean_tc_diffs_rel = time_mean_tc_diffs / (wfde5_time_mean_tc + TFRZ)
    
    # Setup maps
    print('Mapping temperature data to figure...')
    if greenland:
        axes = coordinate_space.nh_horizontal_comparison(map_lon_min=-55,
                                                         map_lon_max=-29,
                                                         map_lat_min=59,
                                                         map_lat_max=84,
                                                         map_lat_0=71.4,
                                                         map_lon_0=-42.1)
    else:
        axes = coordinate_space.nh_horizontal_comparison()
    
    set_map_titles(axes)

    # Map data
    cruncep_quad_mesh = axes[0].pcolor(
                    cruncep_mean_t_rootgrp.variables['LONGXY'][:],
                    cruncep_mean_t_rootgrp.variables['LATIXY'][:],
                    np.ma.clip(cruncep_time_mean_tc, vmin, vmax),
                    shading='nearest',
                    cmap=cmap, vmin=vmin, vmax=vmax,
                    transform=ccrs.PlateCarree())
    wfde5_quad_mesh = axes[1].pcolor(
                    cruncep_mean_t_rootgrp.variables['LONGXY'][:],
                    cruncep_mean_t_rootgrp.variables['LATIXY'][:],
                    np.ma.clip(wfde5_time_mean_tc, vmin, vmax),
                    shading='nearest',
                    cmap=cmap, vmin=vmin, vmax=vmax,
                    transform=ccrs.PlateCarree())
    diffs_quad_mesh = axes[2].pcolor(
                    cruncep_mean_t_rootgrp.variables['LONGXY'][:],
                    cruncep_mean_t_rootgrp.variables['LATIXY'][:],
                    np.ma.clip(time_mean_tc_diffs, -2, 2),
                    shading='nearest',
                    cmap='cet_CET_D4', vmin=-2, vmax=2,
                    transform=ccrs.PlateCarree())
                    
    rel_diffs_quad_mesh = axes[3].pcolor(
                    cruncep_mean_t_rootgrp.variables['LONGXY'][:],
                    cruncep_mean_t_rootgrp.variables['LATIXY'][:],
                    np.ma.clip(100. * time_mean_tc_diffs_rel, -0.5, 0.5),
                    shading='nearest',
                    cmap='cet_CET_D4', vmin=-0.5, vmax=0.5,
                    transform=ccrs.PlateCarree())

    if greenland:
        dem = GrISDEM(path.join('data_raw','gimpdem_90m_v01.1.tif'))
        for i, ax in enumerate(axes):
            # Add evelvation contours
            dem.draw_contours(ax, path.join('data_clean', "gimpdem_90m_v01.1.nc"))
            
        # Colorbar
        fig = plt.gcf()
        cruncep_cbar = fig.colorbar(cruncep_quad_mesh,
                                    ax=axes[0], orientation='vertical')
        cruncep_cbar.set_label('Temperature ($^{\circ}$ C)')
    
        wfde5_cbar = fig.colorbar(wfde5_quad_mesh,
                                    ax=axes[1], orientation='vertical')
        wfde5_cbar.set_label('Temperature ($^{\circ}$ C)')
    
        diffs_cbar = fig.colorbar(diffs_quad_mesh,
                                  ax=axes[2], orientation='vertical')
        diffs_cbar.set_label('Difference (K)')
    
        rel_cbar = fig.colorbar(rel_diffs_quad_mesh,
                                  ax=axes[3], orientation='vertical')
        rel_cbar.set_label(r'Difference ($\%$)')
    
        # Set the figure title
        plt.suptitle('Greenland mean 1980 to 1990 surface air temperature')
    
        # Save results
        print('Writing results')
        plt.savefig(path.join('results', 'greenland_tair_cruncep_vs_wfde5_.png'),
                    dpi=300)
    
    else:
        for i, ax in enumerate(axes):
            # Add evelvation contours
            coordinate_space.draw_elevation_contours(ax)

        # Colorbar
        fig = plt.gcf()
        cruncep_cbar = fig.colorbar(cruncep_quad_mesh,
                                    ax=axes[0:2], orientation='horizontal')
        cruncep_cbar.set_label('Temperature ($^{\circ}$ C)')
    
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
    

        # Set the figure title
        plt.suptitle('Northern Hemisphere mean 1980 to 1990 surface air temperature')
    
        # Save results
        print('Writing results')
        plt.savefig(path.join('results', 'nh_tair_cruncep_vs_wfde5.png'),
                    dpi=300)
    
    # Close figure and files
    plt.close()
    cruncep_mean_t_rootgrp.close()
    wfde5_mean_t_rootgrp.close()
    
def compare_precip(compute_means=True, cmap='cet_CET_L6', vmin=0, vmax=180,
                   Greenland=False):
    """
    """
    # Get CRUNCEP temporal mean precipitation
    cruncep_data = cruncep.CRUNCEP7()
    cruncep_data.get_precip()
    cruncep_mean_precip_rootgrp = cruncep.get_temporal_mean(cruncep_data.precip_rootgrp,
                                                            'PRECTmms',
                                                          compute=compute_means)
    cruncep_data.precip_rootgrp.close()
    # mm H2O / sec -> cm H2O / yr.
    cruncep_time_mean_precip = (60. * 60. * 24. * 365.25 *
                        cruncep_mean_precip_rootgrp.variables['PRECTmms'][:]) / 10.
    
    # Get WFDE5 temporal mean precipitation
    wfde5_data = wfde5.WFDE5()
    # WFDE5 rainfall
    wfde5_data.get_rainf()
    wfde5_mean_rainf_rootgrp = wfde5.get_temporal_mean(wfde5_data.rainf,
                                                       'Rainf',
                                                       compute=compute_means)
    wfde5.close_rootgrps(wfde5_data.rainf)
    # WFDE5 snowfall
    wfde5_data.get_snowf()
    wfde5_mean_snowf_rootgrp = wfde5.get_temporal_mean(wfde5_data.snowf,
                                                       'Snowf',
                                                       compute=compute_means)
    wfde5.close_rootgrps(wfde5_data.snowf)
    # Integrate total precip, convert to cm / yr., and shift WFDE5 data
    # toCRUNCEP grid
    wfde5_time_mean_precip = np.roll((60.* 60. * 24. * 365.25 *
                            (wfde5_mean_rainf_rootgrp.variables['Rainf'][:] +
                             wfde5_mean_snowf_rootgrp.variables['Snowf'][:])) / 10.,
                             360, axis=1)
    # Calculate difference
    print('Computing precipitation differences...')
    time_mean_precip_diffs = cruncep_time_mean_precip - wfde5_time_mean_precip
    time_mean_precip_diffs_rel = time_mean_precip_diffs / wfde5_time_mean_precip
    
    # Setup maps
    if greenland:
        axes = coordinate_space.nh_horizontal_comparison(map_lon_min=-55,
                                                         map_lon_max=-29,
                                                         map_lat_min=59,
                                                         map_lat_max=84,
                                                         map_lat_0=71.4,
                                                         map_lon_0=-42.1)
    else:
        axes = coordinate_space.nh_horizontal_comparison()
    
    set_map_titles(axes)

    # Map data
    print('Mapping precipitation data to figure...')
    cruncep_quad_mesh = axes[0].pcolor(
                    cruncep_mean_precip_rootgrp.variables['LONGXY'][:],
                    cruncep_mean_precip_rootgrp.variables['LATIXY'][:],
                    np.ma.clip(cruncep_time_mean_precip, vmin, vmax),
                    shading='nearest',
                    cmap=cmap, vmin=vmin, vmax=vmax,
                    transform=ccrs.PlateCarree())
    wfde5_quad_mesh = axes[1].pcolor(
                    cruncep_mean_precip_rootgrp.variables['LONGXY'][:],
                    cruncep_mean_precip_rootgrp.variables['LATIXY'][:],
                    np.ma.clip(wfde5_time_mean_precip, vmin, vmax),
                    shading='nearest',
                    cmap=cmap, vmin=vmin, vmax=vmax,
                    transform=ccrs.PlateCarree())
    diffs_quad_mesh = axes[2].pcolor(
                    cruncep_mean_precip_rootgrp.variables['LONGXY'][:],
                    cruncep_mean_precip_rootgrp.variables['LATIXY'][:],
                    np.ma.clip(time_mean_precip_diffs, -50, 50),
                    shading='nearest',
                    cmap='cet_CET_D6_r', vmin=-50, vmax=50,
                    transform=ccrs.PlateCarree())
    rel_diffs_quad_mesh = axes[3].pcolor(
                    cruncep_mean_precip_rootgrp.variables['LONGXY'][:],
                    cruncep_mean_precip_rootgrp.variables['LATIXY'][:],
                    np.ma.clip(100.*time_mean_precip_diffs_rel, -50, 50),
                    shading='nearest',
                    cmap='cet_CET_D6_r', vmin=-50, vmax=50,
                    transform=ccrs.PlateCarree())
                    
    if greenland:
        # Add evelvation contours
        dem = GrISDEM(path.join('data_raw','gimpdem_90m_v01.1.tif'))
        for i, ax in enumerate(axes):
            # Add evelvation contours
            dem.draw_contours(ax, path.join('data_clean', "gimpdem_90m_v01.1.nc"))

        # Colorbar
        fig = plt.gcf()
        cruncep_cbar = fig.colorbar(cruncep_quad_mesh,
                            ax=axes[0], orientation='vertical')
        cruncep_cbar.set_label('Precipitation (cm H$_2$O yr$^{-1}$)')
        
        wfde5_cbar = fig.colorbar(wfde5_quad_mesh,
                            ax=axes[1], orientation='vertical')
        wfde5_cbar.set_label('Precipitation (cm H$_2$O yr$^{-1}$)')
        
        diffs_cbar = fig.colorbar(diffs_quad_mesh,
                            ax=axes[2], orientation='vertical')
        diffs_cbar.set_label('Difference (cm H$_2$O yr$^{-1}$)')
    
        rel_cbar = fig.colorbar(rel_diffs_quad_mesh,
                            ax=axes[3], orientation='vertical')
        rel_cbar.set_label(r'Difference ($\%$)')

        # Set the figure title
        plt.suptitle('Greenland mean 1980 to 1990 precipitation')
    
        # Save results
        print('Writing results')
        plt.savefig(path.join('results', 'greenland_precip_cruncep_vs_wfde5.png'),
                    dpi=300)
    
    else:
        # Add evelvation contours
        for i, ax in enumerate(axes):
            coordinate_space.draw_elevation_contours(ax)

        # Colorbar
        fig = plt.gcf()
        cruncep_cbar = fig.colorbar(cruncep_quad_mesh,
                            ax=axes[0:2], orientation='horizontal')
        cruncep_cbar.set_label('Precipitation rate (cm H$_2$O yr$^{-1}$)')
        '''
        wfde5_cbar = fig.colorbar(wfde5_quad_mesh,
                            ax=axes[1], orientation='horizontal')
        wfde5_cbar.set_label('Precipitation (cm H$_2$O yr$^{-1}$)')
        '''
        diffs_cbar = fig.colorbar(diffs_quad_mesh,
                            ax=axes[2], orientation='horizontal')
        diffs_cbar.set_label('Difference (cm H$_2$O yr$^{-1}$)')
    
        rel_cbar = fig.colorbar(rel_diffs_quad_mesh,
                            ax=axes[3], orientation='horizontal')
        rel_cbar.set_label(r'Difference ($\%$)')

        # Set the figure title
        plt.suptitle('Northern Hemisphere mean 1980 to 1990 precipitation')
    
        # Save results
        print('Writing results')
        plt.savefig(path.join('results', 'nh_precip_cruncep_vs_wfde5.png'), dpi=300)
    
    # Close figure and files
    plt.close()
    wfde5_mean_rainf_rootgrp.close()
    wfde5_mean_snowf_rootgrp.close()
    cruncep_mean_precip_rootgrp.close()
    
def test():
    cruncep.test()
    wfde5.test()

def run():
    #compare_temperature(compute_means=False)
    compare_temperature(compute_means=False, greenland = True, vmin=-32, vmax=6)
    #compare_precip(compute_means=False)
    compare_temperature(compute_means=False, greenland = True, vmin=0, vmax=150)

def main():
    run()

if __name__=='__main__':
    main()
