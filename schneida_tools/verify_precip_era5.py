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
from schneida_tools import era5
#from schneida_tools import cruncep
#from schneida_tools import gswp3
import colorcet as cc

import ipdb
    
def get_era5_temporal_means():
    """
    """
    try:
        era5_mean_precip_rootgrp = era5.get_temporal_mean(None,
                                                          'tp',
                                                          compute=False)
    except FileNotFoundError:
        era5_data = era5.ERA5()
        era5_data.get_precip()
        era5_mean_precip_rootgrp = era5.get_temporal_mean(era5_data.precip_rootgrp,
                                                                'tp')
        gswp3_data.precip_rootgrp.close()
        
    return era5_mean_precip_rootgrp
    
def grid_sumup2era5(xlim=140, ylim=140):
    """ Loop through measurements and filter out:
        1. Measurements outside time period of analysis
        2. All measurements that are not from ice cores
    """
    args = get_args()
    # Get reanalysis data
    era5_mean_precip_rootgrp = get_era5_temporal_means()
    if False:
        # Convert from mm per day to cm per year and shift cruncep data to WFDE5 grid
        cruncep_mean_precip = np.roll((60. * 60. * 24. * 365.25 *
                    cruncep_mean_precip_rootgrp.variables['PRECTmms'][:]) / 10.,
                                     -360, axis=1)
    else: # convert from m per day to cm per year
        era5_mean_precip = ((365. * 10 + 3) *
                    era5_mean_precip_rootgrp.variables['tp'][:]) / 0.1
    # Get SUMup data
    sumup_rootgrp = sumup.get_accumulation()
    
    print('Clustering SUMup dataset for measurements valid from %d to %d'
          % (args.sumup_start_year, args.sumup_stop_year))
    valid_era5_lat_idx = dict()
    valid_era5_lon_idx = dict()
    valid_sumup_lat = dict()
    valid_sumup_lon = dict()
    valid_sumup_accum = dict()
    valid_sumup_error = dict()
    valid_sumup_elev = dict()
    valid_sumup_idxs = list()
    grid_sumup_lat = sumup_rootgrp.variables['Latitude'][:]
    grid_sumup_lon = sumup_rootgrp.variables['Longitude'][:]
    if True:
        # apply longitude correction for [0, 360]
        print('Shifting sumup longitude coordinates from [-180, 180] to '
              '[0, 360]')
        grid_sumup_lon = np.where(grid_sumup_lon > 0, grid_sumup_lon,
                                  grid_sumup_lon + 360)
    
    for i, year in enumerate(sumup_rootgrp.variables['Year'][:]):
        method = sumup_rootgrp.variables['Method'][i]
        if year >= args.sumup_start_year and year < args.sumup_stop_year and method==1.:
            # Valid measurement! Find the nearest grid point
            era5_lat_idx = np.argmin(np.abs(era5_mean_precip_rootgrp.variables['lat'][:] -
                                             grid_sumup_lat[i]))
            era5_lon_idx = np.argmin(np.abs(era5_mean_precip_rootgrp.variables['lon'][:] -
                                             grid_sumup_lon[i]))
            era5_lat = era5_mean_precip_rootgrp.variables['lat'][era5_lat_idx]
            era5_lon = era5_mean_precip_rootgrp.variables['lon'][era5_lon_idx]
            
            # Initialize SUMup dictionarys
            key = '%s,%s' % (str(era5_lat), str(era5_lon))
            
            # Store valid WFDE5 indicies
            valid_era5_lat_idx[key] = era5_lat_idx
            valid_era5_lon_idx[key] = era5_lon_idx
            
            valid_sumup_lat[key] = list()
            valid_sumup_lon[key] = list()
            valid_sumup_accum[key] = list()
            valid_sumup_error[key] = list()
            valid_sumup_elev[key] = list()
            
            # Adjust sumup grid coordinate arrays
            grid_sumup_lat[i] = era5_lat
            grid_sumup_lon[i] = era5_lon
            
            # Store valid SUMup index
            valid_sumup_idxs.append(i)
    
    print('Organizing valid SUMup data into lists...')
    for j, sumup_i in enumerate(valid_sumup_idxs):
        key = '%s,%s' % (str(grid_sumup_lat[sumup_i]),
                         str(grid_sumup_lon[sumup_i]))
        valid_sumup_accum[key].append(float(sumup_rootgrp.variables['Accumulation'][sumup_i]))
        valid_sumup_error[key].append(float(sumup_rootgrp.variables['Error'][sumup_i]))
        valid_sumup_elev[key].append(float(sumup_rootgrp.variables['Elevation'][sumup_i]))
        valid_sumup_lat[key].append(float(grid_sumup_lat[sumup_i]))
        valid_sumup_lon[key].append(float(grid_sumup_lon[sumup_i]))
        if False:
            print(float(grid_sumup_lat[sumup_i]),
                  float(grid_sumup_lon[sumup_i]),
                  float(sumup_rootgrp.variables['Accumulation'][sumup_i]))

    # Setup sample lists
    print('Separating results for Greenland and Antarctica')
    lat_gris_sample = list()
    lon_gris_sample = list()
    gris_n_samples = list()
    ais_n_samples = list()
    lat_ais_sample = list()
    lon_ais_sample = list()
    era5_mean_precip_gris_sample = list()
    era5_mean_precip_ais_sample = list()
    sumup_mean_accum_gris = list()
    sumup_mean_accum_ais = list()
    for key, accum in valid_sumup_accum.items():
        sumup_mean_accum = np.ma.mean(accum)
        if sumup_mean_accum >= 0: # valid accumulation rate
            era5_lat_idx = valid_era5_lat_idx[key]
            era5_lon_idx = valid_era5_lon_idx[key]
            if era5_mean_precip_rootgrp.variables['lat'][era5_lat_idx] > 0 and era5_mean_precip_rootgrp.variables['lon'][era5_lon_idx] > 180:
                # Greenland
                if False:
                    # Store WFDE5 grid points
                    lat_gris_sample.append(era5_mean_precip_rootgrp.variables['lat'][era5_lat_idx])
                    lon_gris_sample.append(era5_mean_precip_rootgrp.variables['lon'][era5_lon_idx])
                else:
                    # Store median SUMup locations
                    sumup_median_lat = np.ma.median(valid_sumup_lat[key])
                    sumup_median_lon = np.ma.median(valid_sumup_lon[key])
                    lat_gris_sample.append(sumup_median_lat)
                    lon_gris_sample.append(sumup_median_lon)
                
                gris_n_samples.append(len(valid_sumup_lat[key]))
                # Convert sumup accumulation rate from m per year to cm per year
                sumup_mean_accum_gris.append(100. * sumup_mean_accum)
                '''
                era5_mean_precip_gris_sample.append(era5_mean_precip[era5_lat_idx,
                                                                       era5_lon_idx])
                '''
                era5_mean_precip_gris_sample.append(era5_mean_precip[era5_lat_idx,
                                                                           era5_lon_idx])
                '''
                gswp3_mean_precip_gris_sample.append(gswp3_mean_precip[wfde5_lat_idx,
                                                                           wfde5_lon_idx])
                '''
            elif era5_mean_precip_rootgrp.variables['lat'][era5_lat_idx] < 0:
                # Antarctica
                if False:
                    # Store WFDE5 grid points
                    lat_ais_sample.append(era5_mean_precip_rootgrp.variables['lat'][era5_lat_idx])
                    lon_ais_sample.append(era5_mean_precip_rootgrp.variables['lon'][era5_lon_idx])
                else:
                    # Store median SUMup locations
                    sumup_median_lat = np.ma.median(valid_sumup_lat[key])
                    sumup_median_lon = np.ma.median(valid_sumup_lon[key])
                    lat_ais_sample.append(sumup_median_lat)
                    lon_ais_sample.append(sumup_median_lon)
                
                ais_n_samples.append(len(valid_sumup_lat[key]))
                # Convert sumup accumulation rate from m per year to cm per year
                sumup_mean_accum_ais.append(100. * sumup_mean_accum)
                '''
                era5_mean_precip_ais_sample.append(era5_mean_precip[era5_lat_idx,
                                                                      era5_lon_idx])
                '''
                era5_mean_precip_ais_sample.append(era5_mean_precip[era5_lat_idx,
                                                                          era5_lon_idx])
                '''
                gswp3_mean_precip_ais_sample.append(gswp3_mean_precip[era5_lat_idx,
                                                                          era5_lon_idx])
                '''
    #wfde5_gris_errors = np.array(wfde5_mean_precip_gris_sample) - np.array(sumup_mean_accum_gris)
    #wfde5_ais_errors = np.array(wfde5_mean_precip_ais_sample) - np.array(sumup_mean_accum_ais)
    era5_gris_errors = np.array(era5_mean_precip_gris_sample) - np.array(sumup_mean_accum_gris)
    era5_ais_errors = np.array(era5_mean_precip_ais_sample) - np.array(sumup_mean_accum_ais)
    #gswp3_gris_errors = np.array(gswp3_mean_precip_gris_sample) - np.array(sumup_mean_accum_gris)
    #gswp3_ais_errors = np.array(gswp3_mean_precip_ais_sample) - np.array(sumup_mean_accum_ais)
    
    gris_sample_matrix = np.array([sumup_mean_accum_gris,
                                  #wfde5_mean_precip_gris_sample,
                                  #era_mean_precip_gris_sample,
                                  era5_mean_precip_gris_sample])
    ais_sample_matrix = np.array([sumup_mean_accum_ais,
                                  #wfde5_mean_precip_ais_sample,
                                  #cruncep_mean_precip_ais_sample,
                                  era5_mean_precip_ais_sample])
    
    gris_covariances, gris_correlations = covariance(gris_sample_matrix)
    ais_covariances, ais_correlations = covariance(ais_sample_matrix)
    
    if True:
        #fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True)
        #axes[0].set_title('Greenland')
        #axes[1].set_title('Antarctica')
        h_gris, xedges_gris, yedges_gris, image_gris= plt.hist2d(sumup_mean_accum_gris,
                                    era5_mean_precip_gris_sample, bins=int(xlim/10.),
                                    range=[[0,xlim],[0,ylim]], density=True,
                                    cmap='cet_linear_worb_100_25_c53')
        h_ais, xedges_ais, yedges_ais, image_ais = plt.hist2d(sumup_mean_accum_ais,
                                     era5_mean_precip_ais_sample, bins=int(xlim/10.),
                                     range=[[0,xlim],[0,ylim]], density=True,
                                     cmap='cet_linear_worb_100_25_c53')
        #plt.savefig(path.join('results', 'hist_era5_sumup_gris_ais_precip.png'), dpi=600)
        plt.close()
        xvals_gris = xedges_gris[:-1] + np.diff(xedges_gris)/2.
        yvals_gris = yedges_gris[:-1] + np.diff(yedges_gris)/2.
    
        xvals_ais = xedges_ais[:-1] + np.diff(xedges_ais)/2.
        yvals_ais = yedges_ais[:-1] + np.diff(yedges_ais)/2.
         
        print('Taylor diagram results (ERA5)')
        print('-----------------------------')
        print('Greenland ice sheet:')
        print('arccos(r) = %r radians' % np.arccos(gris_correlations[0,1]))
        print('std_era5 / std_sumup = %r' % (np.std(era5_mean_precip_gris_sample)/
                                             np.std(sumup_mean_accum_gris)))
        print('RMSE = %r cm/yr' % np.sqrt(np.mean(era5_gris_errors**2)))
        print('Antarctic ice sheet:')
        print('arccos(r) = %r radians' % np.arccos(ais_correlations[0,1]))
        print('std_era5 / std_sumup = %r' % (np.std(era5_mean_precip_ais_sample)/
                                             np.std(sumup_mean_accum_ais)))
        print('RMSE = %r cm/yr' % np.sqrt(np.mean(era5_ais_errors**2)))
    # Scatter data
    fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True)
    #fig.suptitle('Mean 1980 to 1990 precipitation reanalyses vs. accumulation measurements (ice cores)')
    #plt.style.use('agu_half_vertical')
    axes[0].set_title('Greenland')
    axes[1].set_title('Antarctica')
    #axes[0,1].set_title('GSWP3')
    #axes[0,2].set_title('WFDE5')
    axes[0].set_ylabel('ERA5 precipitation (cm w.eq. yr$^{-1}$)')
    axes[1].set_ylabel('ERA5 precipitation (cm w.eq. yr$^{-1}$)')
    axes[1].set_xlabel('SUMup accumulation rate (cm w.eq. yr$^{-1}$)')
    #axes[1,1].set_xlabel('SUMup accumulation rate (cm H$_2$O yr$^{-1}$)')
    #axes[1,2].set_xlabel('SUMup accumulation rate (cm H$_2$O yr$^{-1}$)')
    
    axes[0].set_ylim(0, ylim)
    axes[0].set_xlim(0, xlim)
    
    '''        
    axes[0].errorbar(sumup_mean_accum_gris, era5_mean_precip_gris_sample,
                       yerr=get_yerrs(era5_gris_errors),
                       fmt='x',
                       color='#0064A4',
                       #color='white',
                       ls='')
    
    p_gris = axes[0].contourf(xvals_gris, yvals_gris, h_gris.T, levels=100, cmap='cet_linear_worb_100_25_c53',
                     vmin=0)
    '''
    p_gris = axes[0].plot(sumup_mean_accum_gris, era5_mean_precip_gris_sample,
                             color='#1b3d6d', marker="s",
                             markeredgecolor='None',
                             #(27./255.,61./255.,109./255.,0.1),
                             alpha=0.2, ls='None')
    axes[0].text(xlim/28., ylim - ylim/3., '$n$ = %d\n$r^2$ = %s\nMAE = %s cm yr$^{-1}$\nbias = %s cm yr$^{-1}$'
                            % (len(era5_gris_errors),
                               np.around(gris_correlations[0,1]**2, decimals=4),
                               np.around(np.mean(np.abs(era5_gris_errors)),
                                         decimals=2),
                               np.around(np.mean(era5_gris_errors),
                                         decimals=2),
                               ))
    '''                   
    axes[1].errorbar(sumup_mean_accum_ais, era5_mean_precip_ais_sample,
                       yerr=get_yerrs(era5_ais_errors),
                       fmt='x',
                       color='#0064A4',
                       #color='white',
                       ls='')
    
    p_ais = axes[1].contourf(xvals_ais, yvals_ais, h_ais.T, levels=100, cmap='cet_linear_worb_100_25_c53',
                     vmin=0)
    '''
    p_ais = axes[1].plot(sumup_mean_accum_ais, era5_mean_precip_ais_sample,
                         color='#1b3d6d', marker="s",
                         markeredgecolor='None',#(27./255.,61./255.,109/255.,0.1)
                         alpha=0.2, ls='None')
    axes[1].text(xlim/28., ylim - ylim/3., '$n$ = %d\n$r^2$ = %s\nMAE = %s cm yr$^{-1}$\nbias = %s cm yr$^{-1}$'
                            % (len(era5_ais_errors),
                               np.around(ais_correlations[0,1]**2, decimals=4),
                               np.around(np.mean(np.abs(era5_ais_errors)),
                                         decimals=2),
                               np.around(np.mean(era5_ais_errors),
                                         decimals=2),
                               ))
    fig.colorbar(image_gris, ax=axes[0])
    fig.colorbar(image_ais, ax=axes[1])
    for i in range(2):
        #for j in range(2):
        axes[i].set_ylim(0, ylim)
        axes[i].set_xlim(0, xlim)
        axes[i].plot([0,xlim], [0,ylim], color=(0,0,0,0.2),marker='None')
        axes[i].set_xticks(np.arange(0,xlim+1,20))
        axes[i].set_xticks(np.arange(0,xlim,5), minor=True)
        axes[i].set_yticks(np.arange(0,ylim+1,20))
        axes[i].set_yticks(np.arange(0,ylim,5), minor=True)
        axes[i].tick_params(axis='both', which='both', top=True, right=True)
    '''                           
    axes[0,1].errorbar(sumup_mean_accum_gris, gswp3_mean_precip_gris_sample,
                       yerr=get_yerrs(gswp3_gris_errors),
                       fmt='x', color='#0064A4', ls='')
    axes[0,1].text(0.5*xlim-13, 5, 'R$^2$ = %s\nMAE = %s cm yr$^{-1}$\n$n$ = %d'
                            % (np.around(gris_correlations[0,3]**2, decimals=4),
                               np.around(np.mean(np.abs(gswp3_gris_errors)),
                                         decimals=2),
                               len(gswp3_gris_errors)))
                       
    axes[1,1].errorbar(sumup_mean_accum_ais, gswp3_mean_precip_ais_sample,
                       yerr=get_yerrs(gswp3_ais_errors),
                       fmt='x', color='#0064A4', ls='')
    axes[1,1].text(0.5*xlim-13, 5, 'R$^2$ = %s\nMAE = %s cm yr$^{-1}$\n$n$ = %d'
                            % (np.around(ais_correlations[0,3]**2, decimals=4),
                               np.around(np.mean(np.abs(gswp3_ais_errors)),
                                         decimals=2),
                               len(gswp3_ais_errors)))
                       
    axes[0,2].errorbar(sumup_mean_accum_gris, wfde5_mean_precip_gris_sample,
                       yerr=get_yerrs(wfde5_gris_errors),
                       fmt='x', color='#0064A4', ls='')
    axes[0,2].text(0.5*xlim-13, 5, 'R$^2$ = %s\nMAE = %s cm yr$^{-1}$\n$n$ = %d'
                            % (np.around(gris_correlations[0,1]**2, decimals=4),
                               np.around(np.mean(np.abs(wfde5_gris_errors)),
                                         decimals=2),
                                len(wfde5_gris_errors)))
                       
    axes[1,2].errorbar(sumup_mean_accum_ais, wfde5_mean_precip_ais_sample,
                       yerr=get_yerrs(wfde5_ais_errors),
                       fmt='x', color='#0064A4', ls='')
    axes[1,2].text(0.5*xlim-13, 5, 'R$^2$ = %s\nMAE = %s cm yr$^{-1}$\n$n$ = %d'
                            % (np.around(ais_correlations[0,1]**2, decimals=4),
                               np.around(np.mean(np.abs(wfde5_ais_errors)),
                                         decimals=2),
                               len(wfde5_ais_errors)))
    '''
    plt.savefig(path.join('results', 'p_era5_sumup_gris_ais_precip.png'), dpi=600)
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
    get_era5_temporal_means()
    #get_cruncep_temporal_means()
    #get_gswp3_temporal_means()

def run():
    #test()
    #plt.style.use(path.join('schneida_tools', 'gmd_movie_frame.mplstyle'))
    #plt.style.use('hofmann')
    plt.style.use('agu_online_poster_presentation')
    #plt.style.use('uci_darkblue')
    #plt.style.use('agu_half_horizontal')
    #plt.style.use('uci_blue')
    (sumup_gris, sumup_ais) = grid_sumup2era5()
    ipdb.set_trace()

def main():
    run()

if __name__=='__main__':
    main()
