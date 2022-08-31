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
import merra2

import ipdb

def get_merra2_temporal_means():
    """
    """
    try:
        merra2_mean_precip_rootgrp = merra2.get_temporal_mean(None,
                                                          merra2.TOT_PREC,
                                                          compute=False)
    except FileNotFoundError:
        merra2_data = merra2.MERRA2()
        merra2_data.get_precip()
        merra2_mean_precip_rootgrp = merra2.get_temporal_mean(merra2_data.precip_rootgrp,
                                                                merra2.TOT_PREC)
        
    return merra2_mean_precip_rootgrp
    
def grid_sumup2merra2(xlim=150, ylim=150, axes=None, sublimation_data=None,
                      closefig=False):
    """ Loop through measurements and filter out:
        1. Measurements outside time period of analysis
        2. All measurements that are not from ice cores
    """
    args = get_args()
    # Get reanalysis data
    merra2_mean_precip_rootgrp = get_merra2_temporal_means()
    # Convert from mm per s to cm per year
    merra2_mean_precip = (60. * 60. * 24. * 365. *
                    merra2_mean_precip_rootgrp.variables[merra2.TOT_PREC][:]) / 10.
    # Get SUMup data
    sumup_rootgrp = sumup.get_accumulation()
    
    print('Clustering SUMup dataset for measurements valid from %d to %d'
          % (args.sumup_start_year, args.sumup_stop_year))
    valid_merra2_lat_idx = dict()
    valid_merra2_lon_idx = dict()
    valid_sumup_lat = dict()
    valid_sumup_lon = dict()
    valid_sumup_accum = dict()
    valid_sumup_error = dict()
    valid_sumup_elev = dict()
    valid_sumup_idxs = list()
    grid_sumup_lat = sumup_rootgrp.variables['Latitude'][:]
    grid_sumup_lon = sumup_rootgrp.variables['Longitude'][:]
    if False:
        # apply longitude correction for [0, 360]
        print('Shifting sumup longitude coordinates from [-180, 180] to '
              '[0, 360]')
        grid_sumup_lon = np.where(grid_sumup_lon > 0, grid_sumup_lon,
                                  grid_sumup_lon + 360)
    
    for i, year in enumerate(sumup_rootgrp.variables['Year'][:]):
        method = sumup_rootgrp.variables['Method'][i]
        if year >= args.sumup_start_year and year < args.sumup_stop_year and method==1.:
            # Valid measurement! Find the nearest grid point
            merra2_lat_idx = np.argmin(np.abs(merra2_mean_precip_rootgrp.variables['lat'][:] -
                                             grid_sumup_lat[i]))
            merra2_lon_idx = np.argmin(np.abs(merra2_mean_precip_rootgrp.variables['lon'][:] -
                                             grid_sumup_lon[i]))
            merra2_lat = merra2_mean_precip_rootgrp.variables['lat'][merra2_lat_idx]
            merra2_lon = merra2_mean_precip_rootgrp.variables['lon'][merra2_lon_idx]
            
            # Initialize SUMup dictionarys
            key = '%s,%s' % (str(merra2_lat), str(merra2_lon))
            
            # Store valid MERRA-2 indicies
            valid_merra2_lat_idx[key] = merra2_lat_idx
            valid_merra2_lon_idx[key] = merra2_lon_idx
            
            valid_sumup_lat[key] = list()
            valid_sumup_lon[key] = list()
            valid_sumup_accum[key] = list()
            valid_sumup_error[key] = list()
            valid_sumup_elev[key] = list()
            
            # Adjust sumup grid coordinate arrays
            grid_sumup_lat[i] = merra2_lat
            grid_sumup_lon[i] = merra2_lon
            
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
    merra2_mean_precip_gris_sample = list()
    merra2_mean_precip_ais_sample = list()
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
        
        sumup_mad_accum = stats.median_abs_deviation(np.array(accum) - net_vapor_flux,
                                                     nan_policy='omit')
        sumup_mean_accum = np.ma.median(np.array(accum) - net_vapor_flux)
        if sumup_mean_accum >= 0: # valid accumulation rate
            merra2_lat_idx = valid_merra2_lat_idx[key]
            merra2_lon_idx = valid_merra2_lon_idx[key]
            if merra2_mean_precip_rootgrp.variables['lat'][merra2_lat_idx] > 0:
                # Greenland
                if False:
                    # Store MERRA-2 grid points
                    lat_gris_sample.append(merra2_mean_precip_rootgrp.variables['lat'][merra2_lat_idx])
                    lon_gris_sample.append(merra2_mean_precip_rootgrp.variables['lon'][merra2_lon_idx])
                else:
                    # Store median SUMup locations
                    lat_gris_sample.append(sumup_median_lat)
                    lon_gris_sample.append(sumup_median_lon)
                
                gris_n_samples.append(len(valid_sumup_lat[key]))
                # Convert sumup accumulation rate from m per year to cm per year
                sumup_mean_accum_gris.append(100. * sumup_mean_accum)
                sumup_mad_accum_gris.append(100. * sumup_mad_accum)
                merra2_mean_precip_gris_sample.append(merra2_mean_precip[merra2_lat_idx,
                                                                           merra2_lon_idx])

            elif merra2_mean_precip_rootgrp.variables['lat'][merra2_lat_idx] < 0:
                # Antarctica
                if False:
                    # Store MERRA-2 grid points
                    lat_ais_sample.append(merra2_mean_precip_rootgrp.variables['lat'][merra2_lat_idx])
                    lon_ais_sample.append(merra2_mean_precip_rootgrp.variables['lon'][merra2_lon_idx])
                else:
                    # Store median SUMup locations
                    lat_ais_sample.append(sumup_median_lat)
                    lon_ais_sample.append(sumup_median_lon)
                
                ais_n_samples.append(len(valid_sumup_lat[key]))
                # Convert sumup accumulation rate from m per year to cm per year
                sumup_mean_accum_ais.append(100. * sumup_mean_accum)
                sumup_mad_accum_ais.append(100. * sumup_mad_accum)
                merra2_mean_precip_ais_sample.append(merra2_mean_precip[merra2_lat_idx,
                                                                          merra2_lon_idx])

    merra2_gris_errors = np.array(merra2_mean_precip_gris_sample) - np.array(sumup_mean_accum_gris)
    merra2_ais_errors = np.array(merra2_mean_precip_ais_sample) - np.array(sumup_mean_accum_ais)
    
    gris_sample_matrix = np.array([sumup_mean_accum_gris,
                                  #wfde5_mean_precip_gris_sample,
                                  #era_mean_precip_gris_sample,
                                  merra2_mean_precip_gris_sample])
    ais_sample_matrix = np.array([sumup_mean_accum_ais,
                                  #wfde5_mean_precip_ais_sample,
                                  #cruncep_mean_precip_ais_sample,
                                  merra2_mean_precip_ais_sample])
    
    gris_covariances, gris_correlations = covariance(gris_sample_matrix)
    ais_covariances, ais_correlations = covariance(ais_sample_matrix)
    
    print('Taylor diagram results (MERRA-2)')
    print('--------------------------------')
    print('Greenland ice sheet:')
    print('arccos(r) = %r radians' % np.arccos(gris_correlations[0,1]))
    print('std_merra2; std_sumup = %r; %r cm/yr' %
                    (np.std(merra2_mean_precip_gris_sample),
                     np.std(sumup_mean_accum_gris)))
    print('RMSE = %r cm/yr' % np.sqrt(np.mean(merra2_gris_errors**2)))
    
    print('Antarctic ice sheet:')
    print('arccos(r) = %r radians' % np.arccos(ais_correlations[0,1]))
    print('std_merra2; std_sumup = %r; %r cm/yr' %
                      (np.std(merra2_mean_precip_ais_sample),
                       np.std(sumup_mean_accum_ais)))
    print('RMSE = %r cm/yr' % np.sqrt(np.mean(merra2_ais_errors**2)))
    
    # Scatter data
    if axes is None:
        fig, axes = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)
        fig.suptitle('Mean 1980 to 1990 precipitation reanalyses vs. accumulation measurements (ice cores)')
        axes[0, 0].set_title('MERRA-2')
        #axes[0,1].set_title('GSWP3')
        #axes[0,2].set_title('WFDE5')
    
    axes[-1, 0].set_ylabel('MERRA-2: precipitation\n(cm w.e. yr$^{-1}$)')
    #axes[1, 0].set_ylabel('Antarctic ice sheet\nprecipitation rate (cm w.e. yr$^{-1}$)')
    #axes[1, 0].set_xlabel('SUMup accumulation rate (cm w.e. yr$^{-1}$)')
    #axes[1,1].set_xlabel('SUMup accumulation rate (cm H$_2$O yr$^{-1}$)')
    #axes[1,2].set_xlabel('SUMup accumulation rate (cm H$_2$O yr$^{-1}$)')
    
    #axes[0, 0].set_ylim(0, ylim)
    #axes[0, 0].set_xlim(0, xlim)
    
    '''
    for i in range(2):
        for j in range(2):
            axes[i,j].set_ylim(0, ylim)
            axes[i,j].set_xlim(0, xlim)
            axes[i,j].plot([0,xlim], [0,ylim], color='#555759')
            axes[i,j].set_xticks(np.arange(0,xlim+1,50))
            axes[i,j].set_xticks(np.arange(0,xlim,10), minor=True)
            axes[i,j].set_yticks(np.arange(0,ylim+1,50))
            axes[i,j].set_yticks(np.arange(0,ylim,10), minor=True)
            axes[i,j].tick_params(axis='both', which='both', top=True, right=True)
            
    '''
    axes[-1,0].errorbar(sumup_mean_accum_gris, merra2_mean_precip_gris_sample,
                       xerr = sumup_mad_accum_gris,
                       #yerr=get_yerrs(merra2_gris_errors),
                       #fmt='x',
                       color='#7ab800',
                       marker="o", markeredgecolor='None',
                       alpha=0.25, ls='None')
    axes[-1,0].text(xlim/28., ylim - ylim/2.5, '$n$ = %d\nr$^2$ = %s\nMAE = %s cm yr$^{-1}$\nbias = %s cm yr$^{-1}$'
                            % (len(merra2_gris_errors),
                               np.around(gris_correlations[0,1]**2, decimals=4),
                               np.around(np.mean(np.abs(merra2_gris_errors)),
                                         decimals=2),
                               np.around(np.mean(merra2_gris_errors),
                                         decimals=2)
                               ))
                       
    axes[-1,1].errorbar(sumup_mean_accum_ais, merra2_mean_precip_ais_sample,
                        xerr=sumup_mad_accum_ais,
                        color='#c6beb5', marker="o",
                        markeredgecolor='None', alpha=0.25, ls='None',
                       #yerr=get_yerrs(merra2_ais_errors),
                       #fmt='x',
                       #color='#0064A4',
                       )
    axes[-1,1].text(xlim/28., ylim-ylim/2.5, '$n$ = %d\nr$^2$ = %s\nMAE = %s cm yr$^{-1}$\nbias = %s cm yr$^{-1}$'
                            % (len(merra2_ais_errors),
                               np.around(ais_correlations[0,1]**2, decimals=4),
                               np.around(np.mean(np.abs(merra2_ais_errors)),
                                         decimals=2),
                               np.around(np.mean(merra2_ais_errors),
                                          decimals=2)
                                ))
    if closefig:
        #plt.savefig(path.join('results', 'merra2_sumup_gris_ais_precip.png'), dpi=600)
        plt.close()
    
    return((lat_gris_sample, lon_gris_sample, sumup_mean_accum_gris, gris_n_samples),
           (lat_ais_sample, lon_ais_sample, sumup_mean_accum_ais, ais_n_samples),
           axes)
    
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
    get_merra2_temporal_means()

def run():
    #test()
    #plt.style.use(path.join('schneida_tools', 'gmd_movie_frame.mplstyle'))
    #plt.style.use('hofmann')
    #plt.style.use('agu_online_poster_presentation')
    #plt.style.use('uci_darkblue')
    #plt.style.use('agu_half_horizontal')
    #plt.style.use('uci_blue')
    (sumup_gris, sumup_ais) = grid_sumup2merra2()

def main():
    run()

if __name__=='__main__':
    main()
