#!/usr/bin/env python

from os import path

import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib import pyplot as plt
import colorcet as cc

#from schneida_tools import gswp3
#from schneida_tools import wfde5
from schneida_tools import era5
from schneida_tools import coordinate_space
from schneida_tools import gris_dem
from schneida_tools import ais_dem
from schneida_tools import verify_precip_era5 as verify_precip

import ipdb

TFRZ = 273.15 # K
GRIS_EXTENT = (0, 360, 45, 90)
NH_EXTENT = (0, 360, -30, 90)
AIS_EXTENT = (0, 360, -90, -45)

'''
def set_map_titles(axes):
    axes[0].set_title('GSWP3')
    axes[1].set_title('WFDE5')
    axes[2].set_title('GSWP3 - WFDE5')
    axes[3].set_title('((GSWP3 - WFDE5)\n/ WFDE5) ' + r'$\times$ 100')
'''

def setup_map(greenland=False, antarctica=False,
              lon_lines = np.arange(-180, 180, 30),
              lat_lines = np.arange(-90, 90, 30)):
    """ Uses the Lambert Azimuthal Equal Area map projection
    """
    if greenland:
        map_lon_min = -55
        map_lon_max=-29
        map_lat_min=59
        map_lat_max=84
        map_lat_0=71.4
        map_lon_0=-42.1
    elif antarctica:
        map_lon_min = -180
        map_lon_max=180
        map_lat_min=-90
        map_lat_max=-65
        map_lat_0=-90
        map_lon_0=0
    else: # Northern hempiphere
        map_lon_min=-180
        map_lon_max=180
        map_lat_min=15
        map_lat_max=90
        map_lat_0=90.0
        map_lon_0=-80
        
    map_proj = ccrs.LambertAzimuthalEqualArea(central_longitude=map_lon_0,
                                              central_latitude=map_lat_0,
                                              false_easting=0.0,
                                              false_northing=0.0)
                             
    #plt.style.use(path.join('schneida_tools', 'gmd_movie_frame.mplstyle'))
    
    nrows=1
    ncols=1
    geo_ax = plt.subplot(nrows, ncols, 1, projection=map_proj)

    geo_ax.set_extent((map_lon_min, map_lon_max, map_lat_min, map_lat_max),
                      crs=ccrs.PlateCarree())
    geo_ax.add_feature(cfeature.LAND, color='#C6BEB5')
    geo_ax.add_feature(cfeature.OCEAN, color='#C6BEB5')
    if greenland:
        geo_ax.hlines(-9*10**5, 4*10**5, 4*10**5 + 250*10**3, colors='black', lw=2)
        geo_ax.text(3*10**5, -8.5*10**5, '250 km')
    elif antarctica:
        geo_ax.hlines(-2.7*10**6, -1.5*10**6, -1.5*10**6 + 1000*10**3, colors='black', lw=2)
        geo_ax.text(-1.5*10**6, -2.6*10**6, '1000 km')
        
    return geo_ax

def mask_vals(longxy, latixy, var_arr, greenland=False, antarctica=False):
    """ Mask areas outside map domain
    """
    if greenland:
        extent = GRIS_EXTENT
    elif antarctica:
        extent = AIS_EXTENT 
    else:
        extent = NH_EXTENT
        
    var_arr = np.ma.masked_where(longxy < extent[0], var_arr)
    var_arr = np.ma.masked_where(longxy > extent[1], var_arr)
    var_arr = np.ma.masked_where(latixy < extent[2], var_arr)
    var_arr = np.ma.masked_where(latixy > extent[3], var_arr)
    
    return var_arr

class Analysis(object):
    def __init__(self, compute_means=True, greenland=False, antarctica=False):
        self.compute_means=compute_means
        self.greenland=greenland
        self.antarctica=antarctica
        
    def compare_temperature(self, cmap='cet_CET_L8', degc_min=-7, degc_max=37):
        """
        """
        '''
        # Get GSWP3 temporal mean temperature
        gswp3_data = gswp3.GSWP3()
        gswp3_data.get_tphwl()
        gswp3_mean_t_rootgrp = gswp3.get_temporal_mean(gswp3_data.tphwl_rootgrp,
                                                           'TBOT',
                                                           compute=self.compute_means)
        gswp3_data.tphwl_rootgrp.close()
        # Convert from K to degrees C
        gswp3_time_mean_tc = -TFRZ + gswp3_mean_t_rootgrp.variables['TBOT'][:]
    
        '''
        # Get ERA5 temporal mean temperature
        era5_data = era5.ERA5()
        era5_data.get_tair()
        era5_mean_t_rootgrp = era5.get_temporal_mean(era5_data.t_air, 't2m',
                                                       compute=self.compute_means)
        era5.close_rootgrps(era5_data.t_air)
        if False:
            # Convert to Celcius and shift ERA5 data to CRUNCEP grid
            era5_time_mean_tc = np.roll(-TFRZ + era5_mean_t_rootgrp.variables['t2m'][:],
                                         360, axis=1)
        else:
            # just convert to Celcius
            era5_time_mean_tc = -TFRZ + era5_mean_t_rootgrp.variables['t2m'][:]
        
        '''
        # Calculate difference
        print('Computing temperature differences...')
        time_mean_tc_diffs = gswp3_time_mean_tc - era5_time_mean_tc
        time_mean_tc_diffs_rel = time_mean_tc_diffs / (era5_time_mean_tc + TFRZ)
        '''
    
        # Setup maps
        #axes = coordinate_space.four_map_horizontal_comparison(greenland=self.greenland,
         #                                                      antarctica=self.antarctica)
        ax = setup_map(greenland=self.greenland, antarctica=self.antarctica)
            
        #set_map_titles(axes)
        print('Mapping temperature data to figure...')
        
        # Map data
        longxy, latixy = np.meshgrid(era5_mean_t_rootgrp.variables['lon'][:],
                                     era5_mean_t_rootgrp.variables['lat'][:])
        '''
        longxy = gswp3_mean_t_rootgrp.variables['LONGXY'][:]
        latixy = gswp3_mean_t_rootgrp.variables['LATIXY'][:]
    
        gswp3_quad_mesh = axes[0].pcolormesh(longxy.data, latixy.data,
                                           np.ma.clip(mask_vals(longxy,
                                                                latixy,
                                                                gswp3_time_mean_tc,
                                                                greenland=self.greenland,
                                                                antarctica=self.antarctica),
                                                      degc_min, degc_max),
                                           shading='nearest', cmap=cmap, 
                                           vmin=degc_min, vmax=degc_max,
                                           transform=ccrs.PlateCarree())
        '''
        era5_quad_mesh = ax.pcolormesh(longxy, latixy,
                                      np.ma.clip(mask_vals(longxy,
                                                           latixy,
                                                            era5_time_mean_tc,
                                                            greenland=self.greenland,
                                                            antarctica=self.antarctica),
                                                    degc_min, degc_max),
                                         shading='nearest', cmap=cmap,
                                         vmin=degc_min, vmax=degc_max,
                                         transform=ccrs.PlateCarree())
        '''
        diffs_quad_mesh = axes[2].pcolormesh(longxy.data, latixy.data,
                                         np.ma.clip(mask_vals(longxy,
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
                                                        mask_vals(longxy,
                                                                  latixy,
                                                                  time_mean_tc_diffs_rel,
                                                                  greenland=self.greenland,
                                                                  antarctica=self.antarctica),
                                                        -5, 5),
                                             shading='nearest', cmap='cet_CET_D1',
                                             vmin=-5, vmax=5,
                                             transform=ccrs.PlateCarree())
        
        # Colorbar
        '''
        fig = plt.gcf()
        '''
        gswp3_cbar = fig.colorbar(gswp3_quad_mesh,
                                    ax=axes[0:2], orientation='horizontal')
        gswp3_cbar.set_label('Temperature ($^{\circ}$ C)')

        '''    
        era5_cbar = fig.colorbar(era5_quad_mesh,
                                    ax=ax, orientation='vertical')
        era5_cbar.set_label('Temperature ($^{\circ}$ C)')
        '''    
        
        diffs_cbar = fig.colorbar(diffs_quad_mesh,
                                  ax=axes[2], orientation='horizontal')
        diffs_cbar.set_label('Difference (K)')

        rel_cbar = fig.colorbar(rel_diffs_quad_mesh,
                                  ax=axes[3], orientation='horizontal')
        rel_cbar.set_label(r'Difference ($\%$)')
        '''
        # Draw contours
        if self.greenland:
            dem = gris_dem.GrisDEM()
            '''
            for i, ax in enumerate(axes):
                # Add elevation contours
                dem.draw_contours(ax)
            '''
            dem.draw_contours(ax)
        elif self.antarctica:
            # Add elevation contours
            dem = ais_dem.AisDem()
            '''
            for i, ax in enumerate(axes):
                # Add elevation contours
                dem.draw_contours(ax)
            '''
            dem.draw_contours(ax) 
        
        else:
            '''
            for i, ax in enumerate(axes):
                # Add elevation contours
                coordinate_space.draw_elevation_contours(ax)
            '''
            coordinate_space.draw_elevation_contours(ax)
        #self.gswp3_mean_t_rootgrp = gswp3_mean_t_rootgrp
        
        self.era5_mean_t_rootgrp = era5_mean_t_rootgrp
        
        #return(axes)
        return(ax)
        
    def close_mean_t_rootgrps(self):
        #self.gswp3_mean_t_rootgrp.close()
        self.era5_mean_t_rootgrp.close()
    
    def compare_precip(self,
                       #cmap='cet_CET_L6_r',
                       cmap='cet_CET_L7_r', cm_per_year_min=0,
                       cm_per_year_max=180):
        """
        """
        '''
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
    
        '''
        # Get ERA5 temporal mean precipitation
        era5_data = era5.ERA5()
        # WFDE5 rainfall
        era5_data.get_precip()
        era5_mean_precip_rootgrp = era5.get_temporal_mean(era5_data.precip,
                                                           'tp',
                                                           compute=self.compute_means)
        era5.close_rootgrps(era5_data.precip)
        '''
        # WFDE5 snowfall
        wfde5_data.get_snowf()
        wfde5_mean_snowf_rootgrp = wfde5.get_temporal_mean(wfde5_data.snowf,
                                                           'Snowf',
                                                           compute=self.compute_means)
        wfde5.close_rootgrps(wfde5_data.snowf)
        '''
        if False:
            # Integrate total precip, convert to cm / yr., and shift WFDE5 data
            # to CRUNCEP grid
            era5_time_mean_precip = np.roll((60.* 60. * 24. * 365.25 *
                                (era5_mean_rainf_rootgrp.variables['Rainf'][:] +
                                 era5_mean_snowf_rootgrp.variables['Snowf'][:])) / 10.,
                                 360, axis=1)
        else:
            era5_time_mean_precip = ((365. * 10 + 3) *
                    era5_mean_precip_rootgrp.variables['tp'][:]) / 0.1
        # Calculate difference
        '''
        print('Computing precipitation differences...')
        time_mean_precip_diffs = gswp3_time_mean_precip - wfde5_time_mean_precip
        time_mean_precip_diffs_rel = time_mean_precip_diffs / wfde5_time_mean_precip
        
        # Setup maps
        axes = coordinate_space.four_map_horizontal_comparison(greenland=self.greenland,
                                                               antarctica=self.antarctica)
        set_map_titles(axes)
        '''
        ax = setup_map(greenland=self.greenland, antarctica=self.antarctica)
        # Map data
        print('Mapping precipitation data to figure...')
        longxy, latixy = np.meshgrid(era5_mean_precip_rootgrp.variables['lon'][:],
                                     era5_mean_precip_rootgrp.variables['lat'][:])
        #longxy = gswp3_mean_precip_rootgrp.variables['LONGXY'][:]
        #latixy = gswp3_mean_precip_rootgrp.variables['LATIXY'][:]
        '''
        gswp3_quad_mesh = axes[0].pcolormesh(longxy.data, latixy.data,
                                           np.ma.clip(mask_vals(longxy, latixy,
                                                                gswp3_time_mean_precip,
                                                                greenland=self.greenland,
                                                                antarctica=self.antarctica),
                                                      cm_per_year_min, cm_per_year_max),
                                           shading='nearest', cmap=cmap,
                                           vmin=cm_per_year_min, vmax=cm_per_year_max,
                                           transform=ccrs.PlateCarree())
        '''                               
        era5_quad_mesh = ax.pcolormesh(longxy, latixy,
                                         np.ma.clip(mask_vals(longxy, latixy,
                                                              era5_time_mean_precip,
                                                              greenland=self.greenland,
                                                              antarctica=self.antarctica),
                                                    cm_per_year_min, cm_per_year_max),
                                         shading='nearest', cmap=cmap,
                                         vmin=cm_per_year_min, vmax=cm_per_year_max,
                                         transform=ccrs.PlateCarree())
        '''
        diffs_quad_mesh = axes[2].pcolormesh(longxy.data, latixy.data,
                                         np.ma.clip(mask_vals(longxy, latixy,
                                                              time_mean_precip_diffs,
                                                              greenland=self.greenland,
                                                              antarctica=self.antarctica),
                                                    -50, 50),
                                         shading='nearest', cmap='cet_CET_D10',
                                         vmin=-50, vmax=50,
                                         transform=ccrs.PlateCarree())
                                     
        rel_diffs_quad_mesh = axes[3].pcolormesh(longxy.data, latixy.data,
                                             np.ma.clip(100. *
                                                        mask_vals(longxy, latixy,
                                                                  time_mean_precip_diffs_rel,
                                                                  greenland=self.greenland,
                                                                  antarctica=self.antarctica),
                                                        -50, 50),
                                             shading='nearest', cmap='cet_CET_D10',
                                             vmin=-50, vmax=50,
                                             transform=ccrs.PlateCarree())
        '''
        # Colorbar
        fig = plt.gcf()
        '''
        gswp3_cbar = fig.colorbar(gswp3_quad_mesh,
                            ax=axes[0:2], orientation='horizontal')
        gswp3_cbar.set_label('Precipitation rate (cm H$_2$O yr$^{-1}$)')
        '''
        era5_cbar = fig.colorbar(era5_quad_mesh,
                            ax=ax, orientation='vertical')
        era5_cbar.set_label('Precipitation (cm H$_2$O yr$^{-1}$)')
        '''
        diffs_cbar = fig.colorbar(diffs_quad_mesh,
                            ax=axes[2], orientation='horizontal')
        diffs_cbar.set_label('Difference (cm H$_2$O yr$^{-1}$)')

        rel_cbar = fig.colorbar(rel_diffs_quad_mesh,
                            ax=axes[3], orientation='horizontal')
        rel_cbar.set_label(r'Difference ($\%$)')
        '''
        # Save results
        if self.greenland:
            # Add elevation contours
            print('Drawing Greenland ')
            dem = gris_dem.GrisDEM()
            '''
            for i, ax in enumerate(axes):
                # Add evelvation contours
                dem.draw_contours(ax)
            '''
            dem.draw_contours(ax)
            
        elif self.antarctica:
            # Add elevation contours
            dem = ais_dem.AisDem()
            '''
            for i, ax in enumerate(axes):
                # Add elevation contours
                dem.draw_contours(ax)
            '''
            dem.draw_contours(ax)

        else:
            # Add elevation contours
            '''
            for i, ax in enumerate(axes):
                coordinate_space.draw_elevation_contours(ax)
            '''
            coordinate_space.draw_elevation_contours(ax)

        self.era5_mean_precip_rootgrp = era5_mean_precip_rootgrp
        #self.wfde5_mean_snowf_rootgrp = wfde5_mean_snowf_rootgrp
        #self.gswp3_mean_precip_rootgrp = gswp3_mean_precip_rootgrp
        self.precip_cmap = cmap
        self.precip_cm_per_year_min = cm_per_year_min
        self.precip_cm_per_year_max = cm_per_year_max
        
        return ax
        
    def close_mean_precip_rootgrps(self):
        '''
        self.wfde5_mean_rainf_rootgrp.close()
        self.wfde5_mean_snowf_rootgrp.close()
        self.gswp3_mean_precip_rootgrp.close()
        '''
        self.era5_mean_precip_rootgrp.close()
def test():
    gswp3.test()
    wfde5.test()

def run():
    '''
    plt.style.use(path.join('schneida_tools', 'gmd_movie_frame.mplstyle'))
    plt.style.use('uci_darkblue')
    '''
    plt.style.use('agu_online_poster_presentation')
    #plt.style.use('uci_blue')
    # Greenland analysis
    greenland_analysis = Analysis(compute_means=False,
                                  greenland=True)
    # Temperature
    ax = greenland_analysis.compare_temperature(degc_min=-50, degc_max=0,
                                           #cmap='cet_CET_L7',
                                            cmap='cet_CET_CBL2')
    # Set the figure title
    plt.suptitle('Greenland mean 1980 to 1990 surface air temperature')
    savefig_name = path.join('results', 'greenland_tair_era5.png')
    print('Writing temperature maps to %s' % savefig_name)
    plt.savefig(savefig_name, dpi=300)
    # Close figure and files
    plt.close()
    greenland_analysis.close_mean_t_rootgrps()

    # Precipitation
    (sumup_gris, sumup_ais) = verify_precip.grid_sumup2era5()
    ax = greenland_analysis.compare_precip(cm_per_year_min=0, cm_per_year_max=150)
    # Get and plot SUMup locations
    ax.scatter(sumup_gris[1], sumup_gris[0], s=sumup_gris[3], c=sumup_gris[2],
                    cmap=greenland_analysis.precip_cmap, vmin=greenland_analysis.precip_cm_per_year_min,
                    vmax=greenland_analysis.precip_cm_per_year_max, edgecolors='white',
                    linewidths=0.5,
                    transform=ccrs.PlateCarree())
    '''
    axes[1].scatter(sumup_gris[1], sumup_gris[0], s=sumup_gris[3], c=sumup_gris[2],
                    cmap=greenland_analysis.precip_cmap, vmin=greenland_analysis.precip_cm_per_year_min,
                    vmax=greenland_analysis.precip_cm_per_year_max, edgecolors='white',
                    linewidths=0.5,
                    transform=ccrs.PlateCarree())
    '''
    # Set the figure title
    plt.suptitle('Greenland mean 1980 to 1990 precipitation')
    savefig_name = path.join('results', 'greenland_precip_era5.png')
    print('Writing precipitation maps to %s' % savefig_name)
    plt.savefig(savefig_name, dpi=300)
    # Close figure and files
    plt.close()
    greenland_analysis.close_mean_precip_rootgrps()
    
    # Antarctica analysis
    antarctica_analysis = Analysis(compute_means=False,
                                   antarctica=True)
    # Temperature
    ax = antarctica_analysis.compare_temperature(degc_min=-50, degc_max=0,
                                                   #cmap='cet_CET_L7'
                                                   cmap='cet_CET_CBL2'
                                                  )
    # Set the figure title
    plt.suptitle('Antarctica mean 1980 to 1990 surface air temperature')
    savefig_name = path.join('results', 'antarctica_tair_era5.png')
    print('Writing temperature maps to %s' % savefig_name)
    plt.savefig(savefig_name, dpi=300)
    # Close figure and files
    plt.close()
    antarctica_analysis.close_mean_t_rootgrps()
    
    # Precipitation
    ax = antarctica_analysis.compare_precip(cm_per_year_min=0, cm_per_year_max=150)
    # Get and plot SUMup locations
    ax.scatter(sumup_ais[1], sumup_ais[0], s=sumup_ais[3], c=sumup_ais[2],
                    cmap=antarctica_analysis.precip_cmap, vmin=antarctica_analysis.precip_cm_per_year_min,
                    vmax=antarctica_analysis.precip_cm_per_year_max, edgecolors='white',
                    linewidths=0.5,
                    transform=ccrs.PlateCarree())
    '''
    axes[1].scatter(sumup_ais[1], sumup_ais[0], s=sumup_ais[3], c=sumup_ais[2],
                    cmap=antarctica_analysis.precip_cmap, vmin=antarctica_analysis.precip_cm_per_year_min,
                    vmax=antarctica_analysis.precip_cm_per_year_max, edgecolors='white',
                    linewidths=0.5,
                    transform=ccrs.PlateCarree())
    '''
    # Set the figure title
    plt.suptitle('Antarctica mean 1980 to 1990 precipitation')
    savefig_name = path.join('results', 'antarctica_precip_era5.png')
    plt.savefig(savefig_name, dpi=300)
    # Close figure and files
    plt.close()
    antarctica_analysis.close_mean_precip_rootgrps()
    
    # Northern hemisphere analysis
    northern_hemisphere_analysis = Analysis(compute_means=False)
    ax = northern_hemisphere_analysis.compare_temperature()
    # Set the figure title
    plt.suptitle('Northern Hemisphere mean 1980 to 1990 surface air temperature')
    savefig_name = path.join('results', 'nh_tair_era5.png')
    print('Writing temperature maps to %s' % savefig_name)
    plt.savefig(savefig_name, dpi=300)
    # Close figure and files
    plt.close()
    northern_hemisphere_analysis.close_mean_t_rootgrps()
    
    # Precipitation
    ax = northern_hemisphere_analysis.compare_precip()
    # Set the figure title
    plt.suptitle('Northern Hemisphere mean 1980 to 1990 precipitation')
    savefig_name = path.join('results', 'nh_precip_era5.png')
    plt.savefig(savefig_name, dpi=300)
    # Close figure and files
    plt.close()
    northern_hemisphere_analysis.close_mean_precip_rootgrps()
    
def main():
    run()

if __name__=='__main__':
    main()
