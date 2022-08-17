#!/usr/bin/env python

""" MAR version 3.5.2
"""

from schneida_args import get_args

DAYS_PER_MONTH = np.array([31., 28.25, 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.])

class MARModelDataset(object):
    def __init__(self, verbose=True, plot_ssmb_maps=False,
                 plot_anomalies=False,
                 monthly_cmap=colors.ListedColormap(cc.CET_C10s).reversed()):
        
        args = get_args()
        data_stream = self.get_data(args.mar3gr_raw_data_path,
                                    args.sumup_start_year,
                                    args.sumup_stop_year,
                                    mar_forcing=args.mar3gr_forcing)
        
        #gis_smb = list()
        #total_gis_smb = list()
        monthly_gis_smb = list()
        monthly_gis_precip = list()
        monthly_gis_melt = list()
        monthly_gis_sub = list()
        #monthly_gis_ssmb = list()
        gis_ssmb_time_idx = 0
        year_idx = 0.
        for i, fi in enumerate(data_stream):
            f = Dataset(fi)
            if verbose:
                print(f)
                print(f.variables.keys())
            self.lons = f.variables['LON'][:] + 360
            self.lats = f.variables['LAT'][:]
            '''
            mask = f.variables['MSK_MAR']
            self.topo_5 = f.variables['SRF_bam01'][:]
            self.topo_mar = f.variables['SRF_MAR'][:]
            area = f.variables['AREA']
            smb = f.variables['SMBcorr'] # mmWE per month
            '''
            self.topo_1 = f.variables['SRF_bam13'][:]
            if verbose:
                for d in f.dimensions.items():
                    print(d)
            
            ice_frac = np.ma.masked_where(f.variables['MSK_MAR'][:]<=0,
                                          f.variables['MSK_MAR'][:] - 1.)
            
            if plot_anomalies:
                ice_area = ice_frac * f.variables['AREA'][:] * 10**6 # m^2
                total_ice_area = np.ma.sum(ice_area)
            
                monthly_gis_smb.append(np.ma.sum(ice_area *
                                                 f.variables['SMBcorr'][:],
                                                 axis=(1,2)) / total_ice_area)
                monthly_gis_precip.append(np.ma.sum(ice_area *
                                                    (f.variables['SF'][:] +
                                                     f.variables['RF'][:]),
                                                     axis=(1,2)) / total_ice_area)
                monthly_gis_melt.append(np.ma.sum(ice_area *
                                                 f.variables['MEcorr'][:],
                                                 axis=(1,2)) / total_ice_area)
            
            if plot_ssmb_maps:
                monthly_gis_ssmb_m_yr_T = (365.25 * f.variables['SMBcorr'][:].T) / (1000. 
                                                                         * DAYS_PER_MONTH)
                for j, month_idx in enumerate(f.variables['time'][:]):
                    ssmb_ax = setup_map(central_longitude=-39)
                    self.map2gris(ssmb_ax,
                                  np.ma.masked_where(ice_frac<0.95,
                                                     monthly_gis_ssmb_m_yr_T.T[j]),
                                  cbar_orientation='horizontal', tick_spacing=0.5,
                                  vmin=-1.5, vmax=1.5)
                    plt.title('Reconstructed monthly climatic mass balance\n'
                              '%d CE' % (args.start_year + i))
                    plt.savefig(os.path.join('results', 'graphics', 'gris_ssmb_mar3_%s_%05d.png'
                                             % (args.mar_forcing, gis_ssmb_time_idx)),
                                             dpi=300)
                    plt.close()
                    gis_ssmb_time_idx += 1
                    
            else: # calculate sublimation
                monthly_gis_sub_m_yr_T = (365.25 * f.variables['SUB'][:].T) / (1000.
                                                            * DAYS_PER_MONTH)
                if i==0:
                    mean_annual_gis_sub_m_yr = np.ma.mean(monthly_gis_sub_m_yr_T.T,
                                                          axis=0)
                elif i>0:
                    mean_annual_gis_sub_m_yr += np.ma.mean(monthly_gis_sub_m_yr_T.T,
                                                           axis=0)
            
            f.close()
            year_idx+=1.
        self.mean_gis_sub_m_yr = mean_annual_gis_sub_m_yr / year_idx
            
        if plot_anomalies:
        
            self.monthly_gis_smb_m_yr = (365.25 *
                                         np.ma.asarray(monthly_gis_smb)) / (1000. *
                                                                     DAYS_PER_MONTH)
            self.monthly_gis_smb_anomaly = self.monthly_gis_smb_m_yr - np.ma.mean(
                                            self.monthly_gis_smb_m_yr[:args.ref_years],
                                            axis=0)
            self.monthly_gis_precip_m_yr = (365.25 *
                                         np.ma.asarray(monthly_gis_precip)) / (1000. *
                                                                     DAYS_PER_MONTH)
            self.monthly_gis_precip_anomaly = self.monthly_gis_precip_m_yr - np.ma.mean(
                                            self.monthly_gis_precip_m_yr[:args.ref_years],
                                            axis=0)
            self.monthly_gis_melt_m_yr = (365.25 *
                                         np.ma.asarray(monthly_gis_melt)) / (1000. *
                                                                     DAYS_PER_MONTH)
            self.monthly_gis_melt_anomaly = self.monthly_gis_melt_m_yr - np.ma.mean(
                                            self.monthly_gis_melt_m_yr[:args.ref_years],
                                            axis=0)
        
            for time_idx, smb_anomaly in enumerate(self.monthly_gis_smb_anomaly.flatten()):
                months_2_plot = np.arange(time_idx+1)
                years_ce = self.anl_years[0] + months_2_plot / 12.
                fig, ax_arr = plt.subplots(nrows=3,ncols=1, sharex=True, sharey=True)
                ax_arr[-1].set_xticks(np.arange(self.anl_years[0], self.anl_years[-1]),
                                      minor=True)
                ax_arr[-1].set_yticks(np.arange(-10, 10, 0.5))
                ax_arr[-1].set_yticks(np.arange(-10,10,0.1),minor=True)
                ax_arr[-1].autoscale(axis='y')
                ax_arr[0].autoscale(axis='y')
                ax_arr[1].autoscale(axis='y')
                if years_ce[-1] <= self.anl_years[0] + args.ref_years:
                    ax_arr[-1].set_xlim(self.anl_years[0], self.anl_years[0] + args.ref_years)
                    #ax_arr[0].set_xlim(self.anl_years[0], self.anl_years[0] + args.ref_years)
                    #ax_arr[1].set_xlim(self.anl_years[0], self.anl_years[0] + args.ref_years)
                else:
                    ax_arr[-1].set_xlim(years_ce[-1] - args.ref_years, years_ce[-1])
                    #ax_arr[1].set_xlim(years_ce[-1] - args.ref_years, years_ce[-1])
                    #ax_arr[0].set_xlim(years_ce[-1] - args.ref_years, years_ce[-1])
                
                ax_arr[-1].set_xlabel('year (CE)')
                ax_arr[-1].set_ylabel('climatic mass balance\n(m w.eq. yr$^{-1}$)')
                ax_arr[0].set_ylabel('precipitation rate\n(m w.eq. yr$^{-1}$)')
                ax_arr[1].set_ylabel('melt rate\n(m w.eq. yr$^{-1}$)')
                fig.suptitle('Greenland ice sheet monthly anomalies (vs. %d-%d)' %
                             (self.anl_years[0], self.anl_years[0] + args.ref_years))
                ax_arr[0].hlines(0, self.anl_years[0], self.anl_years[-1], color='black',
                                 linewidths=0.5)
                ax_arr[1].hlines(0, self.anl_years[0], self.anl_years[-1], color='black',
                                 linewidths=0.5)
                ax_arr[2].hlines(0, self.anl_years[0], self.anl_years[-1], color='black',
                                 linewidths=0.5)
            
                ax_arr[0].tick_params(axis='both', which='both', top=True, right=True)
                ax_arr[1].tick_params(axis='both', which='both', top=True, right=True)
                ax_arr[2].tick_params(axis='both', which='both', top=True, right=True)
            
                ax_arr[-1].bar(years_ce,
                               self.monthly_gis_smb_anomaly.flatten()[:time_idx+1],
                               width=1./12., align='edge',
                               color=monthly_cmap((months_2_plot%12) / 12.))
                ax_arr[0].bar(years_ce,
                               self.monthly_gis_precip_anomaly.flatten()[:time_idx+1],
                               width=1./12., align='edge',
                               color=monthly_cmap((months_2_plot%12) / 12.))
                           
                ax_arr[1].bar(years_ce,
                               self.monthly_gis_melt_anomaly.flatten()[:time_idx+1],
                               width=1./12., align='edge',
                               color=monthly_cmap((months_2_plot%12) / 12.))
                '''
                ax.plot(time_monthly[:12*(year_idx+1)], gis_smb_m_yr_filt[:12*(year_idx+1)],
                        color='#ffd200')
                '''
            
                plt.savefig(os.path.join('results', 'graphics', 'gris_anom_mar3_%s_%05d.png'
                                         % (args.mar_forcing, time_idx)),
                                         dpi=300)
                plt.close()
        '''
        self.mean_annual_gis_smb_m_yr = (365.25 / 1000.) * (mean_annual_gis_smb /
                                                            (year_idx + 1.))
        '''
        #ipdb.set_trace()
        
    def get_data(self, input_dir, start_year, stop_year, mar_forcing='NCEP1'):
        data_stream = list()
        self.anl_years = np.arange(start_year, stop_year)
        for i, file_name in enumerate(sorted(os.listdir(input_dir))):
            if fnmatch.fnmatch(file_name, 'MARv3.5.2-20km-monthly-%s-????.nc'
                               % mar_forcing):
                # Data stream file
                year = np.array([int(file_name.split('.')[-2].split('-')[-1])])
                if np.isin(year, self.anl_years)[0]:
                    data_stream.append(os.path.join(input_dir, file_name))

        return data_stream

    def map_smb(self):
        self.ax = setup_map()
        self.map2gris(self.ax, self.mean_annual_gis_smb_m_yr)

    def map2gris(self, ax, var_arr, vmin=-1.4, vmax=1.4, tick_spacing=0.2,
                 cbar_max = None,
                 cmap=colors.ListedColormap(cc.CET_D10),
                 #cmap=colors.ListedColormap(cc.CET_CBTD1),
                 label="specific balance rate (m w.eq. yr$^{-1}$)",
                 cbar_orientation='vertical',
                 ais=False):
        """ Map var_arr onto the given Greenland map (ax)
        """
        if vmin < 0:
            cbar_min = vmin
        else:
            cbar_min = 0
            
        if cbar_max is None:
            cbar_max = vmax
        
        if ais:
            level_max = 4070
            zmax = 4000
        else:
            level_max = 3207
            zmax = 3207
        
        print('Plotting RCM GrIS data...')
        if cbar_orientation=='horizontal' or cbar_orientation=='vertical':
            cbar = plt.colorbar(ScalarMappable(norm=colors.Normalize(vmin=vmin,
                                                                     vmax=vmax,
                                                                     clip=True),
                                           cmap=cmap),
                                           ticks=np.arange(cbar_min, cbar_max+0.001, tick_spacing),
                                           boundaries=np.arange(cbar_min, cbar_max+0.001, tick_spacing/5.),
                                           values=np.arange(cbar_min + 0.5 * (tick_spacing/5.),
                                                            cbar_max, tick_spacing/5.),
                                           label=label,
                                           orientation=cbar_orientation)
        
        
        ax.contourf(self.lons, self.lats,np.ones(var_arr.shape),
                    #np.ma.clip(var_arr,vmin,vmax),
                    #cmap=cmap,
                    #levels=np.arange(vmin,vmax+0.01, 0.01),
                    levels=1,
                    #vmin=vmin, vmax=vmax, #edgecolors='None',
                    colors='#c6beb5',
                    transform=ccrs.PlateCarree())
        '''
        ax.pcolormesh(self.lons, self.lats,
                  var_arr, cmap=cmap,
                  shading='nearest',
                  vmin=vmin, vmax=vmax, edgecolors='None',
                  transform=ccrs.PlateCarree())
        '''
        ax.contourf(self.lons, self.lats, var_arr,
                    cmap=cmap,
                    levels=np.arange(vmin, vmax+0.01, 0.1),
                    extend='both',
                    transform=ccrs.PlateCarree())
        '''              
        ax.contour(self.lons, self.lats,
                   self.topo_mar,
                       levels=np.arange(0, level_max, 500)[1:],
                       cmap=colors.ListedColormap(cc.linear_grey_0_100_c0),
                       linewidths=0.5,
                       linestyles='solid',
                       alpha=1.0,
                       vmin=-zmax,vmax=zmax,
                       transform=ccrs.PlateCarree())
        '''
        ax.contour(self.lons, self.lats,
                   self.topo_1,
                       levels=np.arange(0, level_max, 500),
                       #cmap=colors.ListedColormap(cc.linear_grey_0_100_c0),
                       colors='black',
                       linewidths=0.5,
                       linestyles='solid',
                       alpha=0.5,
                       vmin=-zmax,vmax=zmax,
                       transform=ccrs.PlateCarree())

def run():
    pass

def main():
    run()

if __name__=='__main__':
    main()