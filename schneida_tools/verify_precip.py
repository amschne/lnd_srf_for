#!/usr/bin/env python

""" Verify Greenland ice sheet precipitation data with measurements
    from the Surface Mass Balance and Snow on Sea Ice Working Group (SUMup)
    dataset (Montgomery et al., 2018).

    References

    Montgomery, L., Koenig, L., and Alexander, P.: The SUMup dataset: compiled
        measurements of surface mass balance components over ice sheets and sea ice
        with analysis over Greenland, Earth Syst. Sci. Data, 10, 1959–1985,
        https://doi.org/10.5194/essd-10-1959-2018, 2018.
"""

from os import path

import numpy as np
from scipy import stats
from matplotlib import pyplot as plt

from schneida_args import get_args
import sumup
import wfde5
import cruncep
import gswp3

import ipdb

def get_wfde5_temporal_means():
    """
    """
    # Snowfall
    try:
        wfde5_mean_snowf_rootgrp = wfde5.get_temporal_mean(list(),
                                                           'Snowf',
                                                           compute=False)
    except FileNotFoundError:
        wfde5_data = wfde5.WFDE5()
        wfde5_data.get_snowf()
        wfde5_mean_snowf_rootgrp = wfde5.get_temporal_mean(wfde5_data.snowf,
                                                           'Snowf')
        wfde5.close_rootgrps(wfde5_data.snowf)
        
    # Rainfall
    try:
        wfde5_mean_rainf_rootgrp = wfde5.get_temporal_mean(list(),
                                                           'Rainf',
                                                           compute=False)
    except FileNotFoundError:
        wfde5_data = wfde5.WFDE5()
        wfde5_data.get_rainf()
        wfde5_mean_rainf_rootgrp = wfde5.get_temporal_mean(wfde5_data.rainf,
                                                           'Rainf')
        wfde5.close_rootgrps(wfde5_data.rainf)
        
    return(wfde5_mean_snowf_rootgrp, wfde5_mean_rainf_rootgrp)
    
def get_cruncep_temporal_means():
    """
    """
    try:
        cruncep_mean_precip_rootgrp = cruncep.get_temporal_mean(None,
                                                                'PRECTmms',
                                                                compute=False)
    except FileNotFoundError:
        cruncep_data = cruncep.CRUNCEP7()
        cruncep_data.get_precip()
        cruncep_mean_precip_rootgrp = cruncep.get_temporal_mean(cruncep_data.precip_rootgrp,
                                                                'PRECTmms')
        cruncep_data.precip_rootgrp.close()
        
    return cruncep_mean_precip_rootgrp
    
def get_gswp3_temporal_means():
    """
    """
    try:
        gswp3_mean_precip_rootgrp = gswp3.get_temporal_mean(None,
                                                                'PRECTmms',
                                                                compute=False)
    except FileNotFoundError:
        gswp3_data = gswp3.GSWP3()
        gswp3_data.get_precip()
        gswp3_mean_precip_rootgrp = gswp3.get_temporal_mean(gswp3_data.precip_rootgrp,
                                                                'PRECTmms')
        gswp3_data.precip_rootgrp.close()
        
    return gswp3_mean_precip_rootgrp
    
def grid_sumup2wfde5(xlim=140, ylim=140, axes=None, savefig=True, sublimation_data=None):
    """ Loop through measurements and filter out:
        1. Measurements outside time period of analysis
        2. All measurements that are not from ice cores
    """
    args = get_args()
    # Get reanalysis data
    wfde5_mean_snowf_rootgrp, wfde5_mean_rainf_rootgrp = get_wfde5_temporal_means()
    # Integrate and convert from mm per day to cm per year
    wfde5_mean_precip = (60.* 60. * 24. * 365.25 *
                         (wfde5_mean_snowf_rootgrp.variables['Snowf'][:] +
                          wfde5_mean_rainf_rootgrp['Rainf'][:])) / 10.
    
    cruncep_mean_precip_rootgrp = get_cruncep_temporal_means()
    # Convert from mm per day to cm per year and shift cruncep data to WFDE5 grid
    cruncep_mean_precip = np.roll((60. * 60. * 24. * 365.25 *
                    cruncep_mean_precip_rootgrp.variables['PRECTmms'][:]) / 10.,
                                     -360, axis=1)
                                     
    gswp3_mean_precip_rootgrp = get_gswp3_temporal_means()
    # Convert from mm per day to cm per year and shift gswp3 data to WFDE5 grid
    gswp3_mean_precip = np.roll((60. * 60. * 24. * 365.25 *
                    gswp3_mean_precip_rootgrp.variables['PRECTmms'][:]) / 10.,
                                     -360, axis=1)
    # Get SUMup data
    sumup_rootgrp = sumup.get_accumulation()
    
    print('Clustering SUMup dataset for measurements valid from %d to %d'
          % (args.sumup_start_year, args.sumup_stop_year))
    valid_wfde5_lat_idx = dict()
    valid_wfde5_lon_idx = dict()
    valid_sumup_lat = dict()
    valid_sumup_lon = dict()
    valid_sumup_accum = dict()
    valid_sumup_error = dict()
    valid_sumup_elev = dict()
    valid_sumup_idxs = list()
    grid_sumup_lat = sumup_rootgrp.variables['Latitude'][:]
    grid_sumup_lon = sumup_rootgrp.variables['Longitude'][:]
    for i, year in enumerate(sumup_rootgrp.variables['Year'][:]):
        method = sumup_rootgrp.variables['Method'][i]
        if year >= args.sumup_start_year and year < args.sumup_stop_year and method==1.:
            # Valid measurement! Find the nearest grid point
            wfde5_lat_idx = np.argmin(np.abs(wfde5_mean_snowf_rootgrp.variables['lat'][:] -
                                             sumup_rootgrp.variables['Latitude'][i]))
            wfde5_lon_idx = np.argmin(np.abs(wfde5_mean_snowf_rootgrp.variables['lon'][:] -
                                             sumup_rootgrp.variables['Longitude'][i]))
            wfde5_lat = wfde5_mean_snowf_rootgrp.variables['lat'][wfde5_lat_idx]
            wfde5_lon = wfde5_mean_snowf_rootgrp.variables['lon'][wfde5_lon_idx]
            
            # Initialize SUMup dictionarys
            key = '%s,%s' % (str(wfde5_lat), str(wfde5_lon))
            
            # Store valid WFDE5 indicies
            valid_wfde5_lat_idx[key] = wfde5_lat_idx
            valid_wfde5_lon_idx[key] = wfde5_lon_idx
            
            valid_sumup_lat[key] = list()
            valid_sumup_lon[key] = list()
            valid_sumup_accum[key] = list()
            valid_sumup_error[key] = list()
            valid_sumup_elev[key] = list()
            
            # Adjust sumup grid coordinate arrays
            grid_sumup_lat[i] = wfde5_lat
            grid_sumup_lon[i] = wfde5_lon
            
            # Store valid SUMup index
            valid_sumup_idxs.append(i)
    
    print('Organizing valid SUMup data into lists...')
    for j, sumup_i in enumerate(valid_sumup_idxs):
        key = '%s,%s' % (str(grid_sumup_lat[sumup_i]),
                         str(grid_sumup_lon[sumup_i]))
        valid_sumup_accum[key].append(float(sumup_rootgrp.variables['Accumulation'][sumup_i]))
        valid_sumup_error[key].append(float(sumup_rootgrp.variables['Error'][sumup_i]))
        valid_sumup_elev[key].append(float(sumup_rootgrp.variables['Elevation'][sumup_i]))
        valid_sumup_lat[key].append(float(sumup_rootgrp.variables['Latitude'][sumup_i]))
        valid_sumup_lon[key].append(float(sumup_rootgrp.variables['Longitude'][sumup_i]))
        if False:
            print(float(sumup_rootgrp.variables['Latitude'][sumup_i]),
                  float(sumup_rootgrp.variables['Longitude'][sumup_i]),
                  float(sumup_rootgrp.variables['Accumulation'][sumup_i]))

    # Setup sample lists
    print('Separating results for Greenland and Antarctica')
    lat_gris_sample = list()
    lon_gris_sample = list()
    gris_n_samples = list()
    ais_n_samples = list()
    lat_ais_sample = list()
    lon_ais_sample = list()
    wfde5_mean_precip_gris_sample = list()
    wfde5_mean_precip_ais_sample = list()
    cruncep_mean_precip_gris_sample = list()
    cruncep_mean_precip_ais_sample = list()
    gswp3_mean_precip_gris_sample = list()
    gswp3_mean_precip_ais_sample = list()
    sumup_mean_accum_gris = list()
    sumup_mean_accum_ais = list()
    sumup_mad_accum_gris = list()
    sumup_mad_accum_ais = list()
    for key, accum in valid_sumup_accum.items():
        sumup_median_lat = np.ma.median(valid_sumup_lat[key])
        sumup_median_lon = np.ma.median(valid_sumup_lon[key])
        
        if sublimation_data is None:
            net_vapor_flux = 0.
        else:
            net_vapor_flux = sublimation_data.get_net_vapor_flux(sumup_median_lat,
                                                                 sumup_median_lon)
        
        sumup_mean_accum = np.ma.median(np.array(accum) - net_vapor_flux)
        sumup_mad_accum = stats.median_abs_deviation(np.array(accum) - net_vapor_flux,
                                                     nan_policy='omit')
        if sumup_mean_accum >= 0: # valid accumulation rate
            wfde5_lat_idx = valid_wfde5_lat_idx[key]
            wfde5_lon_idx = valid_wfde5_lon_idx[key]
            if wfde5_mean_snowf_rootgrp.variables['lat'][wfde5_lat_idx] > 0 and wfde5_mean_snowf_rootgrp.variables['lon'][wfde5_lon_idx] < 0:
                # Greenland
                if False:
                    # Store WFDE5 grid points
                    lat_gris_sample.append(wfde5_mean_snowf_rootgrp.variables['lat'][wfde5_lat_idx])
                    lon_gris_sample.append(wfde5_mean_snowf_rootgrp.variables['lon'][wfde5_lon_idx])
                else:
                    # Store median SUMup locations
                    lat_gris_sample.append(sumup_median_lat)
                    lon_gris_sample.append(sumup_median_lon)
                
                gris_n_samples.append(len(valid_sumup_lat[key]))
                # Convert sumup accumulation rate from m per year to cm per year
                sumup_mean_accum_gris.append(100. * sumup_mean_accum)
                sumup_mad_accum_gris.append(100. * sumup_mad_accum)
                wfde5_mean_precip_gris_sample.append(wfde5_mean_precip[wfde5_lat_idx,
                                                                       wfde5_lon_idx])
                cruncep_mean_precip_gris_sample.append(cruncep_mean_precip[wfde5_lat_idx,
                                                                           wfde5_lon_idx])
                gswp3_mean_precip_gris_sample.append(gswp3_mean_precip[wfde5_lat_idx,
                                                                           wfde5_lon_idx])
            elif wfde5_mean_snowf_rootgrp.variables['lat'][wfde5_lat_idx] < 0:
                # Antarctica
                if False:
                    # Store WFDE5 grid points
                    lat_ais_sample.append(wfde5_mean_snowf_rootgrp.variables['lat'][wfde5_lat_idx])
                    lon_ais_sample.append(wfde5_mean_snowf_rootgrp.variables['lon'][wfde5_lon_idx])
                else:
                    # Store median SUMup locations
                    lat_ais_sample.append(sumup_median_lat)
                    lon_ais_sample.append(sumup_median_lon)
                
                ais_n_samples.append(len(valid_sumup_lat[key]))
                # Convert sumup accumulation rate from m per year to cm per year
                sumup_mean_accum_ais.append(100. * sumup_mean_accum)
                sumup_mad_accum_ais.append(100. * sumup_mad_accum)
                wfde5_mean_precip_ais_sample.append(wfde5_mean_precip[wfde5_lat_idx,
                                                                      wfde5_lon_idx])
                cruncep_mean_precip_ais_sample.append(cruncep_mean_precip[wfde5_lat_idx,
                                                                          wfde5_lon_idx])
                gswp3_mean_precip_ais_sample.append(gswp3_mean_precip[wfde5_lat_idx,
                                                                          wfde5_lon_idx])
    
    wfde5_gris_errors = np.array(wfde5_mean_precip_gris_sample) - np.array(sumup_mean_accum_gris)
    wfde5_ais_errors = np.array(wfde5_mean_precip_ais_sample) - np.array(sumup_mean_accum_ais)
    cruncep_gris_errors = np.array(cruncep_mean_precip_gris_sample) - np.array(sumup_mean_accum_gris)
    cruncep_ais_errors = np.array(cruncep_mean_precip_ais_sample) - np.array(sumup_mean_accum_ais)
    gswp3_gris_errors = np.array(gswp3_mean_precip_gris_sample) - np.array(sumup_mean_accum_gris)
    gswp3_ais_errors = np.array(gswp3_mean_precip_ais_sample) - np.array(sumup_mean_accum_ais)
    
    gris_sample_matrix = np.array([sumup_mean_accum_gris,
                                  wfde5_mean_precip_gris_sample,
                                  cruncep_mean_precip_gris_sample,
                                  gswp3_mean_precip_gris_sample])
    ais_sample_matrix = np.array([sumup_mean_accum_ais,
                                  wfde5_mean_precip_ais_sample,
                                  cruncep_mean_precip_ais_sample,
                                  gswp3_mean_precip_ais_sample])
    
    gris_covariances, gris_correlations = covariance(gris_sample_matrix)
    ais_covariances, ais_correlations = covariance(ais_sample_matrix)
    
    if False:
        p_array_cruncep = compute_joint_distribution(sumup_mean_accum_gris,
                                                     sumup_mean_accum_ais,
                                                     cruncep_mean_precip_gris_sample,
                                                     cruncep_mean_precip_ais_sample,
                                                     xlim=xlim, ylim=ylim)
        p_array_wfde5 = compute_joint_distribution(sumup_mean_accum_gris,
                                                   sumup_mean_accum_ais,
                                                   wfde5_mean_precip_gris_sample,
                                                   wfde5_mean_precip_ais_sample,
                                                   xlim=xlim, ylim=ylim)
                                               
        p_array_gswp3 = compute_joint_distribution(sumup_mean_accum_gris,
                                                   sumup_mean_accum_ais,
                                                   gswp3_mean_precip_gris_sample,
                                                   gswp3_mean_precip_ais_sample,
                                                   xlim=xlim, ylim=ylim)
    
    print('Taylor diagram results (CRUNCEP7, GSWP3, WFDE5)')
    print('-----------------------------------------------')
    print('Greenland ice sheet:')
    print('CRUNCEP7 arccos(r) = %r radians' % np.arccos(gris_correlations[0,2]))
    print('std_cruncep; std_sumup = %r; %r cm/yr' %
                    (np.std(cruncep_mean_precip_gris_sample),
                     np.std(sumup_mean_accum_gris)))
    print('CRUNCEP7 RMSE = %r cm/yr' % np.sqrt(np.mean(cruncep_gris_errors**2)))

    print('GSWP3 arccos(r) = %r radians' % np.arccos(gris_correlations[0,3]))
    print('std_gswp3; std_sumup = %r; %r cm/yr' %
                     (np.std(gswp3_mean_precip_gris_sample),
                      np.std(sumup_mean_accum_gris)))
    print('GSWP3 RMSE = %r cm/yr' % np.sqrt(np.mean(gswp3_gris_errors**2)))
        
    print('WFDE5 arccos(r) = %r radians' % np.arccos(gris_correlations[0,1]))
    print('std_wfde5; std_sumup = %r; %r cm/yr' %
                      (np.std(wfde5_mean_precip_gris_sample),
                       np.std(sumup_mean_accum_gris)))
    print('WFDE5 RMSE = %r cm/yr' % np.sqrt(np.mean(wfde5_gris_errors**2)))
    
    print('Antarctic ice sheet:')
    print('CRUNCEP7 arccos(r) = %r radians' % np.arccos(ais_correlations[0,2]))
    print('std_cruncep; std_sumup = %r; %r cm/yr' %
                    (np.std(cruncep_mean_precip_ais_sample),
                     np.std(sumup_mean_accum_ais)))
    print('CRUNCEP7 RMSE = %r cm/yr' % np.sqrt(np.mean(cruncep_ais_errors**2)))

    print('GSWP3 arccos(r) = %r radians' % np.arccos(ais_correlations[0,3]))
    print('std_gswp3; std_sumup = %r; %r cm/yr' %
                    (np.std(gswp3_mean_precip_ais_sample),
                     np.std(sumup_mean_accum_ais)))
    print('GSWP3 RMSE = %r cm/yr' % np.sqrt(np.mean(gswp3_ais_errors**2)))
        
    print('WFDE5 arccos(r) = %r radians' % np.arccos(ais_correlations[0,1]))
    print('std_wfde5; std_sumup = %r; %r cm/yr' %
                   (np.std(wfde5_mean_precip_ais_sample),
                    np.std(sumup_mean_accum_ais)))
    print('WFDE5 RMSE = %r cm/yr' % np.sqrt(np.mean(wfde5_ais_errors**2)))
    # Scatter data
    if axes is None:
        fig, axes = plt.subplots(nrows=5, ncols=2, sharex=True, sharey=False)
    #fig.suptitle('Mean 1980 to 1990 precipitation reanalyses vs. accumulation measurements (ice cores)')
    #axes[0,0].set_title('Greenland')
    axes[2,0].set_ylabel('CRUNCEP: precipitation\n(cm w.eq. yr$^{-1}$)')
    #axes[1,0].set_title('Antarctica')
    #axes[1,0].set_ylabel('CRUNCEP precipitation (cm w.eq. yr$^{-1}$)')
    
    #axes[0,1].set_title('Greenland')
    axes[3,0].set_ylabel('GSWP3: precipitation\n(cm w.eq. yr$^{-1}$)')
    #axes[1,1].set_title('Antarctica')
    #axes[1,1].set_ylabel('GSWP3 precipitation (cm w.eq. yr$^{-1}$)')
    
    #axes[0,2].set_title('Greenland')
    axes[1,0].set_ylabel('WFDE5: precipitation\n(cm w.eq. yr$^{-1}$)')
    #axes[1,2].set_title('Antarctica')
    #axes[1,2].set_ylabel('WFDE5 precipitation (cm w.eq. yr$^{-1}$)')
    
    #axes[1,0].set_xlabel('SUMup accumulation rate (cm w.eq. yr$^{-1}$)')
    #axes[1,1].set_xlabel('SUMup accumulation rate (cm w.eq. yr$^{-1}$)')
    #axes[1,2].set_xlabel('SUMup accumulation rate (cm w.eq. yr$^{-1}$)')
    
    axes[2,0].errorbar(sumup_mean_accum_gris, cruncep_mean_precip_gris_sample,
                       xerr=sumup_mad_accum_gris,
                       #yerr=get_yerrs(cruncep_gris_errors),
                       #fmt='x', color='#0064A4', ls='')
                   color='black', marker="o",
                   markeredgecolor='None',
                       #(27./255.,61./255.,109./255.,0.1),
                   alpha=0.25, ls='None'
                   )
    '''
    p_gris_cruncep = axes[0,0].contourf(p_array_cruncep['gris'][0], p_array_cruncep['gris'][1],
                                p_array_cruncep['gris'][2], levels=100,
                                cmap='cet_linear_worb_100_25_c53',
                                vmin=0)
    '''# gris_correlations[0,2]
    axes[2,0].text(xlim/1.5 - 5,ylim/20.-2,
          '$n$ = %d\n$r^2$ = %s\nMAE = %s cm yr$^{-1}$\nbias = %s cm yr$^{-1}$'
                   % (len(cruncep_gris_errors),
                      np.around(gris_correlations[0,2]**2, decimals=4),
                      np.around(np.mean(np.abs(cruncep_gris_errors)),
                                decimals=2),
                      np.around(np.mean(cruncep_gris_errors),
                                         decimals=2),
                        ))
                       
    
    axes[2,1].errorbar(sumup_mean_accum_ais, cruncep_mean_precip_ais_sample,
                       xerr=sumup_mad_accum_ais,
                       #yerr=get_yerrs(cruncep_ais_errors),
                   #fmt='x', color='#0064A4', ls='')
                   color='#f78d2d', marker="o",
                   markeredgecolor='None',
                   #(27./255.,61./255.,109./255.,0.1),
                   alpha=0.25, ls='None')
    '''
    
    p_ais_cruncep = axes[1,0].contourf(p_array_cruncep['ais'][0], p_array_cruncep['ais'][1],
                                p_array_cruncep['ais'][2], levels=100,
                                cmap='cet_linear_worb_100_25_c53',
                                vmin=0)
    
    '''#ais_correlaitons[0,2]
    axes[2,1].text(xlim/28., ylim - ylim/2.5,
          '$n$ = %d\n$r^2$ = %s\nMAE = %s cm yr$^{-1}$\nbias = %s cm yr$^{-1}$'
                   % (len(cruncep_ais_errors),
                      np.around(ais_correlations[0,2]**2, decimals=4),
                      np.around(np.mean(np.abs(cruncep_ais_errors)),
                                decimals=2),
                      np.around(np.mean(cruncep_ais_errors),
                                         decimals=2),
                        ))
                               
    axes[3,0].errorbar(sumup_mean_accum_gris, gswp3_mean_precip_gris_sample,
                       xerr=sumup_mad_accum_gris,
                       #yerr=get_yerrs(gswp3_gris_errors),
                       #fmt='x', color='#0064A4', ls='')
                       color='#555759', marker="o",
                   markeredgecolor='None',
                   #(27./255.,61./255.,109./255.,0.1),
                   alpha=0.25, ls='None')
    '''
    p_gris_gswp3 = axes[0,1].contourf(p_array_gswp3['gris'][0], p_array_gswp3['gris'][1],
                                p_array_gswp3['gris'][2], levels=100,
                                cmap='cet_linear_worb_100_25_c53',
                                vmin=0)
                                
    '''# gris_correlations[0,3]
    axes[3,0].text(xlim/1.5 - 5, ylim/20.-2,
          '$n$ = %d\n$r^2$ = %s\nMAE = %s cm yr$^{-1}$\nbias = %s cm yr$^{-1}$'
                   % (len(gswp3_gris_errors),
                      np.around(gris_correlations[0,3]**2, decimals=4),
                      np.around(np.mean(np.abs(gswp3_gris_errors)),
                                decimals=2),
                      np.around(np.mean(gswp3_gris_errors),
                                         decimals=2),
                        ))
                       
    
    axes[3,1].errorbar(sumup_mean_accum_ais, gswp3_mean_precip_ais_sample,
                       xerr=sumup_mad_accum_ais,
                       #yerr=get_yerrs(gswp3_ais_errors),
                       #fmt='x', color='#0064A4', ls='')
                   color='#6aa2b8', marker="o",
                   markeredgecolor='None',
                   #(27./255.,61./255.,109./255.,0.1),
                   alpha=0.25, ls='None')
    '''
    p_ais_gswp3 = axes[1,1].contourf(p_array_gswp3['ais'][0], p_array_gswp3['ais'][1],
                                p_array_gswp3['ais'][2], levels=100,
                                cmap='cet_linear_worb_100_25_c53',
                                vmin=0)
    '''# ais_correlations[0,3]
    axes[3,1].text(xlim/28., ylim - ylim/2.5,
          '$n$ = %d\n$r^2$ = %s\nMAE = %s cm yr$^{-1}$\nbias = %s cm yr$^{-1}$'
                   % (len(gswp3_ais_errors),
                      np.around(ais_correlations[0,3]**2, decimals=4),
                      np.around(np.mean(np.abs(gswp3_ais_errors)),
                                decimals=2),
                      np.around(np.mean(gswp3_ais_errors),
                                         decimals=2),
                        ))
    
    axes[1,0].errorbar(sumup_mean_accum_gris, wfde5_mean_precip_gris_sample,
                   xerr=sumup_mad_accum_gris,
                       #yerr=get_yerrs(wfde5_gris_errors),
                       #fmt='x', color='#0064A4', ls='')
                             color='#1b3d6d', marker="o",
                             markeredgecolor='None',
                             #(27./255.,61./255.,109./255.,0.1),
                             alpha=0.25, ls='None')
    '''

    p_gris_wfde5 = axes[0,2].contourf(p_array_wfde5['gris'][0], p_array_wfde5['gris'][1],
                                      p_array_wfde5['gris'][2], levels=100,
                                      cmap='cet_linear_worb_100_25_c53',
                                      vmin=0)
    
    '''# gris_correlations[0,1]
    axes[1,0].text(xlim/1.5 - 5,ylim/20.-2,
          '$n$ = %d\n$r^2$ = %s\nMAE = %s cm yr$^{-1}$\nbias = %s cm yr$^{-1}$'
                   % (len(wfde5_gris_errors),
                      np.around(gris_correlations[0,1]**2, decimals=4),
                      np.around(np.mean(np.abs(wfde5_gris_errors)),
                                decimals=2),
                      np.around(np.mean(wfde5_gris_errors),
                                         decimals=2),
                        ))
                       
    
    axes[1,1].errorbar(sumup_mean_accum_ais, wfde5_mean_precip_ais_sample,
                       xerr=sumup_mad_accum_ais,
                       #yerr=get_yerrs(wfde5_ais_errors),
                       #fmt='x', color='#0064A4', ls='')
                             color='#f7eb5f', marker="o",
                             markeredgecolor='None',
                             #(27./255.,61./255.,109./255.,0.1),
                             alpha=0.25, ls='None')
    '''
    p_ais_wfde5 = axes[1,2].contourf(p_array_wfde5['ais'][0], p_array_wfde5['ais'][1],
                                      p_array_wfde5['ais'][2], levels=100,
                                      cmap='cet_linear_worb_100_25_c53',
                                      vmin=0)
    
    '''#ais_correlations[0,1]
    axes[1,1].text(xlim/28., ylim - ylim/2.5,
          '$n$ = %d\n$r^2$ = %s\nMAE = %s cm yr$^{-1}$\nbias = %s cm yr$^{-1}$'
                   % (len(wfde5_ais_errors),
                      np.around(ais_correlations[0,1]**2, decimals=4),
                      np.around(np.mean(np.abs(wfde5_ais_errors)),
                                decimals=2),
                      np.around(np.mean(wfde5_ais_errors),
                                         decimals=2),
                        ))
    '''
    for i in range(2):
        for j in range(3):
            axes[i,j].set_ylim(0, ylim)
            axes[i,j].set_xlim(0, xlim)
            axes[i,j].plot([0,xlim], [0,ylim], color=(0,0,0,0.2), marker='None')
            axes[i,j].set_xticks(np.arange(0,xlim+1,20))
            axes[i,j].set_xticks(np.arange(0,xlim,5), minor=True)
            axes[i,j].set_yticks(np.arange(0,ylim+1,20))
            axes[i,j].set_yticks(np.arange(0,ylim,5), minor=True)
            axes[i,j].tick_params(axis='both', which='both', top=True, right=True)
    '''
    if savefig:
        plt.savefig(path.join('results', 'p_cruncep_wfde5_sumup_gris_ais_precip.png'), dpi=600)
    plt.close()
    
    return((lat_gris_sample, lon_gris_sample, sumup_mean_accum_gris, gris_n_samples),
           (lat_ais_sample, lon_ais_sample, sumup_mean_accum_ais, ais_n_samples),
            axes)
    
def compute_joint_distribution(sumup_mean_accum_gris, sumup_mean_accum_ais,
                               ra_mean_precip_gris_sample, ra_mean_precip_ais_sample,
                               xlim=140, ylim=140):
    # Compute join distribution of precipitation renanalyses and SUMup
    #.accumulation rates
    h_gris, xedges_gris, yedges_gris, image_gris= plt.hist2d(sumup_mean_accum_gris,
                                    ra_mean_precip_gris_sample, bins=int(xlim/10.),
                                    range=[[0,xlim],[0,ylim]], density=True,
                                    cmap='cet_linear_worb_100_25_c53')
    h_ais, xedges_ais, yedges_ais, image_ais = plt.hist2d(sumup_mean_accum_ais,
                                     ra_mean_precip_ais_sample, bins=int(xlim/10.),
                                     range=[[0,xlim],[0,ylim]], density=True,
                                     cmap='cet_linear_worb_100_25_c53')
    plt.close()
    
    xvals_gris = xedges_gris[:-1] + np.diff(xedges_gris)/2.
    yvals_gris = yedges_gris[:-1] + np.diff(yedges_gris)/2.
    
    xvals_ais = xedges_ais[:-1] + np.diff(xedges_ais)/2.
    yvals_ais = yedges_ais[:-1] + np.diff(yedges_ais)/2.
    
    return({'gris':[xvals_gris, yvals_gris, h_gris.T],
            'ais':[xvals_ais, yvals_ais, h_ais.T]})

def covariance(sample_matrix):
    sample_means = sample_matrix.mean(axis=1)
    sample_std = sample_matrix.std(axis=1, ddof=1)
    sample_covariance = np.empty((sample_matrix.shape[0],
                                    sample_matrix.shape[0]))
                                    
    sample_correlation = np.empty((sample_matrix.shape[0],
                                    sample_matrix.shape[0]))
                                    
    for i in range(sample_matrix.shape[0]):
        for j in range(sample_matrix.shape[0]):
            #print('Computing covariance (%d, %d)' % (i, j))
            sample_covariance[i,j] = (sample_matrix[i] -
                                      sample_means[i]) @ (sample_matrix[j] -
                                                          sample_means[j])
            
            sample_correlation[i,j] = ((sample_matrix[i] - sample_means[i]) @
                                       (sample_matrix[j] - sample_means[j]) / 
                                       (sample_std[i] * sample_std[j]))
                                         
                                         
    sample_covariance = (1. / (sample_matrix.shape[1] - 1)) * sample_covariance
    sample_correlation = (1. / (sample_matrix.shape[1] - 1)) * sample_correlation

    return(sample_covariance, sample_correlation)

def get_yerrs(err_arr):
    yerr = np.zeros((2, err_arr.size))
    for i, err in enumerate(err_arr):
        if err > 0:
            yerr[0, i] = err
        elif err < 0:
            yerr[1, i] = -err
            
    return yerr

def test():
    get_wfde5_temporal_means()
    get_cruncep_temporal_means()
    get_gswp3_temporal_means()

def run():
    #test()
    #plt.style.use(path.join('schneida_tools', 'gmd_movie_frame.mplstyle'))
    #plt.style.use('hofmann')
    plt.style.use('agu_online_poster_presentation')
    #plt.style.use('uci_darkblue')
    #plt.style.use('agu_half_horizontal')
    (sumup_gris, sumup_ais) = grid_sumup2wfde5()
    ipdb.set_trace()

def main():
    run()

if __name__=='__main__':
    main()
