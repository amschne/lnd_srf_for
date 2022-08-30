#!/usr/bin/env python

""" Taylor diagram of mean 1980-1990 precipitation rates over ice sheets
    wrt SUMup dataset
"""

import os

from matplotlib import pyplot as plt

import numpy as np
import mpl_toolkits.axisartist.floating_axes as fa
import mpl_toolkits.axisartist.grid_finder as gf
from matplotlib.projections import PolarAxes

AREA_FACTOR = 15 #/ 25.

class TaylorDiagram(object):
    """
    Taylor diagram.
    Plot model standard deviation and correlation to reference (data)
    sample in a single-quadrant polar plot, with r=stddev and
    theta=arccos(correlation).
    """

    def __init__(self, refstd=1,
                 fig=None, rect=111, label='_', srange=(0, 1.55), extend=False):
        """
        Set up Taylor diagram axes, i.e. single quadrant polar
        plot, using `mpl_toolkits.axisartist.floating_axes`.
        Parameters:
        * refstd: reference standard deviation to be compared to
        * fig: input Figure or None
        * rect: subplot definition
        * label: reference label
        * srange: stddev axis extension, in units of *refstd*
        * extend: extend diagram to negative correlations
        """

        self.refstd = refstd            # Reference standard deviation

        tr = PolarAxes.PolarTransform()

        # Correlation labels
        rlocs = np.array([0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99, 1.])
        if extend:
            # Diagram extended to negative correlations
            self.tmax = np.pi
            rlocs = np.concatenate((-rlocs[:0:-1], rlocs))
        else:
            # Diagram limited to positive correlations
            self.tmax = np.pi/2
        tlocs = np.arccos(rlocs)        # Conversion to polar angles
        gl1 = gf.FixedLocator(tlocs)    # Positions
        mydict = dict(zip(tlocs, map(str, rlocs)))
        for key,val in mydict.items():
            mydict[key] = val[:4]
        
        tf1 = gf.DictFormatter(mydict)

        # Standard deviation axis extent (in units of reference stddev)
        self.smin = srange[0] #* self.refstd
        self.smax = srange[1] #* self.refstd

        ghelper = fa.GridHelperCurveLinear(
            tr,
            extremes=(0, self.tmax, self.smin, self.smax),
            grid_locator1=gl1, tick_formatter1=tf1)

        if fig is None:
            fig = plt.figure()

        ax = fa.FloatingSubplot(fig, rect, grid_helper=ghelper)
        fig.add_subplot(ax)

        # Adjust axes
        ax.axis["top"].set_axis_direction("bottom")   # "Angle axis"
        ax.axis["top"].toggle(ticklabels=True, label=True)
        ax.axis["top"].major_ticklabels.set_axis_direction("top")
        ax.axis["top"].label.set_axis_direction("top")
        ax.axis["top"].label.set_text("correlation")

        ax.axis["left"].set_axis_direction("bottom")  # "X axis"
        ax.axis["left"].label.set_text("standard deviation (cm w.eq. yr$^{-1}$)")
        

        ax.axis["right"].set_axis_direction("top")    # "Y-axis"
        ax.axis["right"].toggle(ticklabels=True)
        ax.axis["right"].major_ticklabels.set_axis_direction(
            "bottom" if extend else "left")

        if self.smin:
            ax.axis["bottom"].toggle(ticklabels=False, label=False)
        else:
            ax.axis["bottom"].set_visible(False)          # Unused

        self._ax = ax                   # Graphical axes
        self.ax = ax.get_aux_axes(tr)   # Polar coordinates

        # Add reference point and stddev contour
        l, = self.ax.plot([0], self.refstd, '*',
                          ls='', ms=10, label=label,
                          color='#555759')
        t = np.linspace(0, self.tmax)
        r = np.zeros_like(t) + self.refstd
        self.ax.plot(t, r, ls='dotted', color=(0,0,0,0.75),label='_')

        # Collect sample points for latter use (e.g. legend)
        self.samplePoints = [l]
        
        self.t_linspace = t

def get_gris_results():
    """ Results derived separately from analysis_era5 and analysis_gswp3 programs
    """
    
    sumup_radius = [24.049718323938052, #23.719029486103437, # ERA5 cm /yr
                    23.996971676499307, # 23.371683440537733, # MERRA-2
                    24.051362434541037, # 23.719029486103437, # CRUNCEP
                    24.051362434541037, # 23.719029486103437, # GSWP3
                    24.051362434541037] # 23.719029486103437] # WFDE5

    
    era5_theta = 0.2045482908902662 #0.2048826648895172 # radians
    era5_radius = 26.893792683943417 # cm /yr
    
    merra2_theta = 0.21896414957472363 # 0.28582511855431747 # radians
    merra2_radius = 22.574283660341457 # 22.6284310247389 # cm/yr
    
    wfde5_theta = 0.9503020169170465 # 0.9377025075172972 # radians
    wfde5_radius = 35.271312887883006 # cm per year
    
    gswp3_theta = 1.0392600516288513# 1.0234720512343265 # radians
    gswp3_radius = 43.06229791795151 # cm /yr
    
    cruncep_theta = 1.0389109501180485 # 1.0171079721060172 # radians
    cruncep_radius = 22.5894250371343 # cm /yr
    
    return({'era5': (era5_theta, era5_radius),
            'merra2':(merra2_theta, merra2_radius),
            'wfde5': (wfde5_theta, wfde5_radius),
            'gswp3': (gswp3_theta, gswp3_radius),
            'cruncep': (cruncep_theta, cruncep_radius),
            'sumup': (0, np.mean(sumup_radius))})
    
def get_ais_results():

    sumup_radius = [11.214891060861696,#11.269603662650837, # ERA5 cm/yr
                    11.337603133420226, # 11.496933568622762, # MERRA2 cm/yr
                    11.347537483099194, # 11.505834868534699, # CRUNCEP cm/yr
                    11.347537483099194, # 11.505834868534699, # GWSP3
                    11.347537483099194] # 11.505834868534699] # WFDE5
                    
    era5_theta = 0.5420903918917168 # 0.5816203992277824 # radians
    era5_radius = 10.987223390578112 # cm /yr
    
    merra2_theta = 0.5176767766930164 # 0.5591699951113855 # radians
    merra2_radius = 10.54466748441476 # cm/yr
    
    wfde5_theta =  0.5547841649037506 # 0.6035660695589737 # radians
    wfde5_radius = 10.63814578363175 # cm/yr
    
    gswp3_theta = 0.784882815993962 # 0.8255956397327044 # radians
    gswp3_radius = 12.86432345537742
    
    cruncep_theta = 1.19774845566386 # 1.2607223588974623# radians
    cruncep_radius = 13.558739442821484 # cm/yr
    
    return({'era5': (era5_theta, era5_radius),
            'merra2':(merra2_theta, merra2_radius),
            'wfde5': (wfde5_theta, wfde5_radius),
            'gswp3': (gswp3_theta, gswp3_radius),
            'cruncep': (cruncep_theta, cruncep_radius),
            'sumup': (0, np.mean(sumup_radius))})

def run():
    gris_results = get_gris_results()
    ais_results = get_ais_results()
    
    plt.style.use('agu_quarter')
    #plt.style.use('grl')
    taylor_diagram = TaylorDiagram(refstd=gris_results['sumup'][1],
                                   srange=(0,45))
    ax = taylor_diagram.ax
    ax.plot([0], ais_results['sumup'][1], 'X',
            ls='', ms=10, color='#c6beb5')
    ax.plot(taylor_diagram.t_linspace,
            np.zeros_like(taylor_diagram.t_linspace) + ais_results['sumup'][1],
            ls='dashed', color=(0,0,0,0.25),label='_')
    
    # WFDE5
    ax.scatter(gris_results['wfde5'][0], gris_results['wfde5'][1],
               s=AREA_FACTOR*33,
               c='black', marker='$\mathrm{WE5_G^*}$')
    ax.scatter(gris_results['wfde5'][0], gris_results['wfde5'][1],
               c='#1b3d6d', s= 0.5*AREA_FACTOR*33, alpha=0.5,
               marker='o')
    ax.scatter(ais_results['wfde5'][0], ais_results['wfde5'][1],
               s=AREA_FACTOR *60,
               c='black', marker='$\mathrm{WE5_A^{x}}$')
    ax.scatter(ais_results['wfde5'][0], ais_results['wfde5'][1],
               c='#f7eb5f', s=0.5*AREA_FACTOR * 60, alpha=0.5,
               marker='o')

    # CRUNCEP7
    ax.scatter(gris_results['cruncep'][0], gris_results['cruncep'][1],
               s=AREA_FACTOR*33,
               c='black', marker='$\mathrm{CN7_G^*}$')
    ax.scatter(gris_results['cruncep'][0], gris_results['cruncep'][1],
               c='black', s= 0.5*AREA_FACTOR*33, alpha=0.5,
               marker='o')
    ax.scatter(ais_results['cruncep'][0], ais_results['cruncep'][1],
               s=AREA_FACTOR *60,
               c='black', marker='$\mathrm{CN7_A^{x}}$')
    ax.scatter(ais_results['cruncep'][0], ais_results['cruncep'][1],
               c='#f78d2d', s=0.5*AREA_FACTOR * 60, alpha=0.5,
               marker='o')
    
    # GSWP3
    ax.scatter(gris_results['gswp3'][0], gris_results['gswp3'][1],
               s=AREA_FACTOR * 33,
               c='black', marker='$\mathrm{GP3_G^*}$')
    ax.scatter(gris_results['gswp3'][0], gris_results['gswp3'][1],
               c='#555759', s=0.5*AREA_FACTOR * 33, alpha=0.5,
               marker='o')
    ax.scatter(ais_results['gswp3'][0], ais_results['gswp3'][1],
               s=AREA_FACTOR * 60,
               c='black', marker='$\mathrm{GP3_A^{x}}$')
    ax.scatter(ais_results['gswp3'][0], ais_results['gswp3'][1],
               c='#6aa2b8', s=0.5*AREA_FACTOR * 60, alpha=0.5,
               marker='o')
               
    # MERRA-2
    ax.scatter(gris_results['merra2'][0], gris_results['merra2'][1],
               c='black', s=AREA_FACTOR*34,
               marker='$\mathrm{MR2_G^*}$')
    ax.scatter(gris_results['merra2'][0], gris_results['merra2'][1],
               c='#7ab800', s=0.5*AREA_FACTOR*34, alpha=0.5,
               marker='o')
    ax.scatter(ais_results['merra2'][0], ais_results['merra2'][1],
               s=AREA_FACTOR*61,
               c='black', marker='$\mathrm{MR2_A^{x}}$')
    ax.scatter(ais_results['merra2'][0], ais_results['merra2'][1],
               c='#c6beb5', s=0.5*AREA_FACTOR*61, alpha=0.5,
               marker='o')
               
    # ERA5
    ax.scatter(gris_results['era5'][0], gris_results['era5'][1],
               c='black', s=AREA_FACTOR*33,
               marker='$\mathrm{ER5_G^*}$')
    ax.scatter(gris_results['era5'][0], gris_results['era5'][1],
               c='#0064a4', s=0.5*AREA_FACTOR*33, alpha=0.5,
               marker='o')
    ax.scatter(ais_results['era5'][0], ais_results['era5'][1],
               s=AREA_FACTOR*65,
               c='black', marker='$\mathrm{ER5_A^{x}}$')
    ax.scatter(ais_results['era5'][0], ais_results['era5'][1],
               c='#ffd200', s=0.5*AREA_FACTOR*65, alpha=0.5,
               marker='o')
               

    plt.grid(alpha=0.25)
    #ax.set_thetamin(0)
    #ax.set_thetamax(90)

    plt.savefig(os.path.join('results', 'taylor_e5_we_g3_cn_m2_sumup.pdf'))

def main():
    run()

if __name__=='__main__':
    main()