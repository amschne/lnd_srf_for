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
    axes[2].set_title('CRUNCEP7 - WFDE5')

def compare_temperature(sample_step=1, vmin=-30, vmax=30):
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
    file_counter = 0
    for i, wfde5_month in enumerate(wfde5_data.t_air):
        month_mean = np.ma.mean(wfde5_month.variables['Tair'][::sample_step],
                                axis=0)
        wfde5_time_mean_tc += month_mean
        file_counter += 1
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
                    cmap='cet_CET_D4', vmin=vmin, vmax=vmax,
                    transform=ccrs.PlateCarree())
    wfde5_quad_mesh = axes[1].pcolor(
                    cruncep_data.tphwl_rootgrp.variables['LONGXY'][:],
                    cruncep_data.tphwl_rootgrp.variables['LATIXY'][:],
                    np.ma.clip(wfde5_time_mean_tc, vmin, vmax),
                    shading='nearest',
                    cmap='cet_CET_D4', vmin=vmin, vmax=vmax,
                    transform=ccrs.PlateCarree())
    diffs_quad_mesh = axes[2].pcolor(
                    cruncep_data.tphwl_rootgrp.variables['LONGXY'][:],
                    cruncep_data.tphwl_rootgrp.variables['LATIXY'][:],
                    np.ma.clip(time_mean_tc_diffs, vmin, vmax),
                    shading='nearest',
                    cmap='cet_CET_D4', vmin=vmin, vmax=vmax,
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
                    
def test():
    cruncep.test()
    wfde5.test()

def run():
    compare_temperature()

def main():
    run()

if __name__=='__main__':
    main()
