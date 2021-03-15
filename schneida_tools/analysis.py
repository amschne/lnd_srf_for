#!/usr/bin/env python

from os import path

import numpy as np
import cartopy.crs as ccrs
from matplotlib import pyplot as plt
import colorcet as cc

from schneida_tools import cruncep
from schneida_tools import wfde5
from schneida_tools import coordinate_space

import ipdb

TFRZ = 273.15 # K

def set_map_titles(axes):
    axes[0].set_title('CRUNCEP7')
    axes[1].set_title('WFDE5')
    axes[2].set_title(r'(CRUNCEP7 - WFDE5) $\times$ 10')

def compare_temperature(sample_step=1, cmap='cet_CET_D4', vmin=-30, vmax=30):
    cruncep_data = cruncep.CRUNCEP7()
    cruncep_data.get_tphwl()
    wfde5_data = wfde5.WFDE5()
    wfde5_data.get_tair()
    #ipdb.set_trace()
    # Calculate temporal means
    print('Computing temporal means...')
    cruncep_time_mean_tc = -TFRZ + np.ma.mean(
                    cruncep_data.tphwl_rootgrp.variables['TBOT'][::sample_step],
                     axis=0)
    
    wfde5_init_tair = wfde5_data.t_air[0]
    wfde5_time_mean_tc = np.ma.zeros(wfde5_init_tair.variables['Tair'][0].shape)
    file_counter = 0.
    for i, wfde5_month in enumerate(wfde5_data.t_air):
        month_mean = np.ma.mean(wfde5_month.variables['Tair'][::sample_step],
                                axis=0)
        wfde5_time_mean_tc += month_mean
        file_counter += 1.
    # Average, convert to Celcius, and shift WFDE5 data to CRUNCEP grid
    wfde5_time_mean_tc = np.roll(-TFRZ + (wfde5_time_mean_tc / file_counter),
                                 360, axis=1)
    
    # Calculate difference
    print('Computing differences...')
    time_mean_tc_diffs = cruncep_time_mean_tc - wfde5_time_mean_tc
    
    # Setup maps
    axes = coordinate_space.nh_horizontal_comparison()
    set_map_titles(axes)

    # Map data
    print('Mapping data to figure...')
    cruncep_quad_mesh = axes[0].pcolor(
                    cruncep_data.tphwl_rootgrp.variables['LONGXY'][:],
                    cruncep_data.tphwl_rootgrp.variables['LATIXY'][:],
                    np.ma.clip(cruncep_time_mean_tc, vmin, vmax),
                    shading='nearest',
                    cmap=cmap, vmin=vmin, vmax=vmax,
                    transform=ccrs.PlateCarree())
    wfde5_quad_mesh = axes[1].pcolor(
                    cruncep_data.tphwl_rootgrp.variables['LONGXY'][:],
                    cruncep_data.tphwl_rootgrp.variables['LATIXY'][:],
                    np.ma.clip(wfde5_time_mean_tc, vmin, vmax),
                    shading='nearest',
                    cmap=cmap, vmin=vmin, vmax=vmax,
                    transform=ccrs.PlateCarree())
    diffs_quad_mesh = axes[2].pcolor(
                    cruncep_data.tphwl_rootgrp.variables['LONGXY'][:],
                    cruncep_data.tphwl_rootgrp.variables['LATIXY'][:],
                    np.ma.clip(10*time_mean_tc_diffs, vmin, vmax),
                    shading='nearest',
                    cmap=cmap, vmin=vmin, vmax=vmax,
                    transform=ccrs.PlateCarree())

    # Add evelvation contours
    for i, ax in enumerate(axes):
        coordinate_space.draw_elevation_contours(ax)

    # Colorbar
    fig = plt.gcf()
    cbar = fig.colorbar(cruncep_quad_mesh,
                        ax=axes[:], orientation='horizontal')
    cbar.set_label('Temperature ($^{\circ}$ C)')

    # Set the figure title
    plt.suptitle('Northern Hemisphere mean 1980 to 1990 surface air temperature')
    
    # Save results
    print('Writing results')
    plt.savefig(path.join('results', 'tair_cruncep_vs_wfde5.png'), dpi=300)
    
    # Close figure and files
    plt.close()
    wfde5.close_rootgrps(wfde5_data.t_air)
    cruncep.tphwl_rootgrp.close()
    
def compare_precip(sample_step=1, cmap='cet_CET_D6_r', vmin=-150, vmax=150):
    cruncep_data = cruncep.CRUNCEP7()
    cruncep_data.get_precip()
    wfde5_data = wfde5.WFDE5()
    wfde5_data.get_rainf()
    wfde5_data.get_snowf()
    #ipdb.set_trace()
    # Calculate temporal means
    print('Computing temporal means...')
    # mm H2O / sec -> cm H2O / yr.
    cruncep_time_mean_precip = (60. * 60. * 24. * 365.25 * np.ma.mean(
                    cruncep_data.precip_rootgrp.variables['PRECTmms'][::sample_step],
                     axis=0) / 10.)
    
    # WFDE5 rainfall
    wfde5_init_rainf = wfde5_data.rainf[0]
    wfde5_time_mean_rainf = np.ma.zeros(wfde5_init_rainf.variables['Rainf'][0].shape)
    file_counter = 0.
    for i, wfde5_month in enumerate(wfde5_data.rainf):
        month_mean = np.ma.mean(wfde5_month.variables['Rainf'][::sample_step],
                                axis=0)
        wfde5_time_mean_rainf += month_mean
        file_counter += 1.
    wfde5_time_mean_rainf = wfde5_time_mean_rainf / file_counter
    
    # WFDE5 snowfall
    wfde5_init_snowf = wfde5_data.snowf[0]
    wfde5_time_mean_snowf = np.ma.zeros(wfde5_init_snowf.variables['Snowf'][0].shape)
    file_counter = 0.
    for i, wfde5_month in enumerate(wfde5_data.snowf):
        month_mean = np.ma.mean(wfde5_month.variables['Snowf'][::sample_step],
                                axis=0)
        wfde5_time_mean_snowf += month_mean
        file_counter += 1.
    wfde5_time_mean_snowf = wfde5_time_mean_snowf / file_counter
    
    # Integrate total precip, convert to cm / yr., and shift WFDE5 data
    # toCRUNCEP grid
    wfde5_time_mean_precip = np.roll((60.* 60. * 24. * 365.25 *
                                     (wfde5_time_mean_rainf +
                                      wfde5_time_mean_snowf)) / 10.,
                                     360, axis=1)
    # Calculate difference
    print('Computing differences...')
    time_mean_precip_diffs = cruncep_time_mean_precip - wfde5_time_mean_precip
    
    # Setup maps
    axes = coordinate_space.nh_horizontal_comparison()
    set_map_titles(axes)

    # Map data
    print('Mapping data to figure...')
    cruncep_quad_mesh = axes[0].pcolor(
                    cruncep_data.precip_rootgrp.variables['LONGXY'][:],
                    cruncep_data.precip_rootgrp.variables['LATIXY'][:],
                    np.ma.clip(cruncep_time_mean_precip, vmin, vmax),
                    shading='nearest',
                    cmap=cmap, vmin=vmin, vmax=vmax,
                    transform=ccrs.PlateCarree())
    wfde5_quad_mesh = axes[1].pcolor(
                    cruncep_data.precip_rootgrp.variables['LONGXY'][:],
                    cruncep_data.precip_rootgrp.variables['LATIXY'][:],
                    np.ma.clip(wfde5_time_mean_precip, vmin, vmax),
                    shading='nearest',
                    cmap=cmap, vmin=vmin, vmax=vmax,
                    transform=ccrs.PlateCarree())
    diffs_quad_mesh = axes[2].pcolor(
                    cruncep_data.precip_rootgrp.variables['LONGXY'][:],
                    cruncep_data.precip_rootgrp.variables['LATIXY'][:],
                    np.ma.clip(10*time_mean_precip_diffs, vmin, vmax),
                    shading='nearest',
                    cmap=cmap, vmin=vmin, vmax=vmax,
                    transform=ccrs.PlateCarree())

    # Add evelvation contours
    for i, ax in enumerate(axes):
        coordinate_space.draw_elevation_contours(ax)

    # Colorbar
    fig = plt.gcf()
    cbar = fig.colorbar(cruncep_quad_mesh,
                        ax=axes[:], orientation='horizontal')
    cbar.set_label('Precipitation rate (cm H$_2$O yr$^{-1})')

    # Set the figure title
    plt.suptitle('Northern Hemisphere mean 1980 to 1990 precipitation')
    
    # Save results
    print('Writing results')
    plt.savefig(path.join('results', 'precip_cruncep_vs_wfde5.png'), dpi=300)
    
    # Close figure and files
    plt.close()
    wfde5.close_rootgrps(wfde5_data.precip)
    cruncep.precip_rootgrp.close()
    
def test():
    cruncep.test()
    wfde5.test()

def run():
    #compare_temperature()
    compare_precip()

def main():
    run()

if __name__=='__main__':
    main()
