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
from matplotlib import pyplot as plt

from schneida_tools.schneida_args import get_args
from schneida_tools import sumup
from schneida_tools import wfde5
from schneida_tools import cruncep
from schneida_tools import gswp3

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
        gswp3_data = gswp3.CRUNCEP7()
        gswp3_data.get_precip()
        gswp3_mean_precip_rootgrp = gswp3.get_temporal_mean(gswp3_data.precip_rootgrp,
                                                                'PRECTmms')
        gswp3_data.precip_rootgrp.close()
        
    return gswp3_mean_precip_rootgrp
    
def grid_sumup2wfde5(xlim=150, ylim=150):
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
    print('Separating results for Greenland and Antartica')
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
    for key, accum in valid_sumup_accum.items():
        sumup_mean_accum = np.ma.mean(accum)
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
                    sumup_median_lat = np.ma.median(valid_sumup_lat[key])
                    sumup_median_lon = np.ma.median(valid_sumup_lon[key])
                    lat_gris_sample.append(sumup_median_lat)
                    lon_gris_sample.append(sumup_median_lon)
                
                gris_n_samples.append(len(valid_sumup_lat[key]))
                # Convert sumup accumulation rate from m per year to cm per year
                sumup_mean_accum_gris.append(100. * sumup_mean_accum)
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
                    sumup_median_lat = np.ma.median(valid_sumup_lat[key])
                    sumup_median_lon = np.ma.median(valid_sumup_lon[key])
                    lat_ais_sample.append(sumup_median_lat)
                    lon_ais_sample.append(sumup_median_lon)
                
                ais_n_samples.append(len(valid_sumup_lat[key]))
                # Convert sumup accumulation rate from m per year to cm per year
                sumup_mean_accum_ais.append(100. * sumup_mean_accum)
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
                                  gswp3_mean_precip_ais_sample]])
    
    gris_covariances, gris_correlations = covariance(gris_sample_matrix)
    ais_covariances, ais_correlations = covariance(ais_sample_matrix)
    
    # Scatter data
    fig, axes = plt.subplots(nrows=2, ncols=3, sharex=True, sharey=True)
    fig.suptitle('Mean 1980 to 1990 precipitation reanalyses vs. accumulation measurements (ice cores)')
    axes[0,0].set_title('CRUNCEP7')
    axes[0,1].set_title('GSWP3')
    axes[0,2].set_title('WFDE5')
    axes[0,0].set_ylabel('Greenland ice sheet\nprecipitation rate (cm H$_2$O yr$^{-1}$)')
    axes[1,0].set_ylabel('Antarctic ice sheet\nprecipitation rate (cm H$_2$O yr$^{-1}$)')
    axes[1,0].set_xlabel('SUMup accumulation rate (cm H$_2$O yr$^{-1}$)')
    axes[1,1].set_xlabel('SUMup accumulation rate (cm H$_2$O yr$^{-1}$)')
    axes[1,2].set_xlabel('SUMup accumulation rate (cm H$_2$O yr$^{-1}$)')
    
    axes[0,0].set_ylim(0, ylim)
    axes[0,0].set_xlim(0, xlim)
    
    for i in range(2):
        for j in range(3):
            axes[i,j].set_ylim(0, ylim)
            axes[i,j].set_xlim(0, xlim)
            axes[i,j].plot([0,xlim], [0,ylim], color='#555759')
            axes[i,j].set_xticks(np.arange(0,xlim+1,50))
            axes[i,j].set_xticks(np.arange(0,xlim,10), minor=True)
            axes[i,j].set_yticks(np.arange(0,ylim+1,50))
            axes[i,j].set_yticks(np.arange(0,ylim,10), minor=True)
            axes[i,j].tick_params(axis='both', which='both', top=True, right=True)
            
    axes[0,0].errorbar(sumup_mean_accum_gris, cruncep_mean_precip_gris_sample,
                       yerr=get_yerrs(cruncep_gris_errors),
                       fmt='x', color='#0064A4', ls='')
    axes[0,0].text(0.5*xlim, 5, 'R$^2$ = %s\nMAE = %s cm yr$^{-1}$\n$n$ = %d'
                            % (np.around(gris_correlations[0,2]**2, decimals=4),
                               np.around(np.mean(np.abs(cruncep_gris_errors)),
                                         decimals=2),
                               len(cruncep_gris_errors)))
                       
    axes[1,0].errorbar(sumup_mean_accum_ais, cruncep_mean_precip_ais_sample,
                       yerr=get_yerrs(cruncep_ais_errors),
                       fmt='x', color='#0064A4', ls='')
    axes[1,0].text(0.5*xlim, 5, 'R$^2$ = %s\nMAE = %s cm yr$^{-1}$\n$n$ = %d'
                            % (np.around(ais_correlations[0,2]**2, decimals=4),
                               np.around(np.mean(np.abs(cruncep_ais_errors)),
                                         decimals=2),
                               len(cruncep_ais_errors)))
                               
    axes[0,1].errorbar(sumup_mean_accum_gris, gswp3_mean_precip_gris_sample,
                       yerr=get_yerrs(gswp3_gris_errors),
                       fmt='x', color='#0064A4', ls='')
    axes[0,1].text(0.5*xlim, 5, 'R$^2$ = %s\nMAE = %s cm yr$^{-1}$\n$n$ = %d'
                            % (np.around(gris_correlations[0,3]**2, decimals=4),
                               np.around(np.mean(np.abs(gswp3_gris_errors)),
                                         decimals=2),
                               len(gswp3_gris_errors)))
                       
    axes[1,1].errorbar(sumup_mean_accum_ais, gswp3_mean_precip_ais_sample,
                       yerr=get_yerrs(gswp3_ais_errors),
                       fmt='x', color='#0064A4', ls='')
    axes[1,1].text(0.5*xlim, 5, 'R$^2$ = %s\nMAE = %s cm yr$^{-1}$\n$n$ = %d'
                            % (np.around(ais_correlations[0,3]**2, decimals=4),
                               np.around(np.mean(np.abs(gswp3_ais_errors)),
                                         decimals=2),
                               len(gswp3_ais_errors)))
                       
    axes[0,2].errorbar(sumup_mean_accum_gris, wfde5_mean_precip_gris_sample,
                       yerr=get_yerrs(wfde5_gris_errors),
                       fmt='x', color='#0064A4', ls='')
    axes[0,2].text(0.5*xlim, 5, 'R$^2$ = %s\nMAE = %s cm yr$^{-1}$\n$n$ = %d'
                            % (np.around(gris_correlations[0,1]**2, decimals=4),
                               np.around(np.mean(np.abs(wfde5_gris_errors)),
                                         decimals=2),
                                len(wfde5_gris_errors)))
                       
    axes[1,2].errorbar(sumup_mean_accum_ais, wfde5_mean_precip_ais_sample,
                       yerr=get_yerrs(wfde5_ais_errors),
                       fmt='x', color='#0064A4', ls='')
    axes[1,2].text(0.5*xlim, 5, 'R$^2$ = %s\nMAE = %s cm yr$^{-1}$\n$n$ = %d'
                            % (np.around(ais_correlations[0,1]**2, decimals=4),
                               np.around(np.mean(np.abs(wfde5_ais_errors)),
                                         decimals=2),
                               len(wfde5_ais_errors)))
    plt.savefig(path.join('results', 'cruncep_wfde5_sumup_gris_ais_precip.pdf'))
    plt.close()
    
    return((lat_gris_sample, lon_gris_sample, sumup_mean_accum_gris, gris_n_samples),
           (lat_ais_sample, lon_ais_sample, sumup_mean_accum_ais, ais_n_samples))
    
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
