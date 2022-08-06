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

AREA_FACTOR = 1. / 25.

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
        self.smin = srange[0] * self.refstd
        self.smax = srange[1] * self.refstd

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
        ax.axis["left"].label.set_text("standard deviation (normalized)")
        

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
        l, = self.ax.plot([0], self.refstd, 'kX',
                          ls='', ms=10, label=label)
        t = np.linspace(0, self.tmax)
        r = np.zeros_like(t) + self.refstd
        self.ax.plot(t, r, ls='dotted', color=(0,0,0,0.5),label='_')

        # Collect sample points for latter use (e.g. legend)
        self.samplePoints = [l]

def get_gris_results():
    """ Results derived separately from analysis_era5 and analysis_gswp3 programs
    """
    
    era5_theta = 0.19073151285457224 # radians
    era5_radius = 1.2459804340913663
    
    wfde5_theta = 0.9283789072489057 # radians
    wfde5_radius = 1.6341081475375194
    
    gswp3_theta = 1.0128260796391035 # radians
    gswp3_radius = 1.9950618822467103
    
    cruncep_theta = 1.0046885784830488 # radians
    cruncep_radius = 1.0465605184220506
    
    return({'era5': (era5_theta, era5_radius),
            'wfde5': (wfde5_theta, wfde5_radius),
            'gswp3': (gswp3_theta, gswp3_radius),
            'cruncep': (cruncep_theta, cruncep_radius)})
    
def get_ais_results():
    era5_theta = 0.5795554925778409 # radians
    era5_radius = 0.9684034966207717
    
    wfde5_theta = 0.5793364633318597 # radians
    wfde5_radius = 0.9237726825519608
    
    gswp3_theta = 0.8255956397327044 # radians
    gswp3_radius = 1.1170847654555394
    
    cruncep_theta = 1.2612001662696082 # radians
    cruncep_radius = 1.1773849843635356
    
    return({'era5': (era5_theta, era5_radius),
            'wfde5': (wfde5_theta, wfde5_radius),
            'gswp3': (gswp3_theta, gswp3_radius),
            'cruncep': (cruncep_theta, cruncep_radius)})

def run():
    plt.style.use('agu_quarter')
    taylor_diagram = TaylorDiagram(srange=(0,2.2))
    ax = taylor_diagram.ax

    gris_results = get_gris_results()
    ais_results = get_ais_results()
    
    ax.scatter(gris_results['era5'][0], gris_results['era5'][1],
               c='black', s=2*AREA_FACTOR * np.pi * 33**2,
               marker='$\mathrm{E5}$')
    ax.scatter(gris_results['era5'][0], gris_results['era5'][1],
               c='#1b3d6d', s=AREA_FACTOR * np.pi * 33**2, alpha=0.5,
               marker='o')
    ax.scatter(ais_results['era5'][0], ais_results['era5'][1],
               s=2*AREA_FACTOR * np.pi * 65**2,
               c='black', marker='$\mathrm{E5}$')
    ax.scatter(ais_results['era5'][0], ais_results['era5'][1],
               c='#f78d2d', s=AREA_FACTOR * np.pi * 65**2, alpha=0.5,
               marker='o')
               
    ax.scatter(gris_results['wfde5'][0], gris_results['wfde5'][1],
               s=2*AREA_FACTOR * np.pi * 33**2,
               c='black', marker='$\mathrm{WE}$')
    ax.scatter(gris_results['wfde5'][0], gris_results['wfde5'][1],
               c='#1b3d6d', s=AREA_FACTOR * np.pi * 33**2, alpha=0.5,
               marker='o')
    ax.scatter(ais_results['wfde5'][0], ais_results['wfde5'][1],
               s=2*AREA_FACTOR * np.pi * 60**2,
               c='black', marker='$\mathrm{WE}$')
    ax.scatter(ais_results['wfde5'][0], ais_results['wfde5'][1],
               c='#f78d2d', s=AREA_FACTOR * np.pi * 60**2, alpha=0.5,
               marker='o')
               
    ax.scatter(gris_results['gswp3'][0], gris_results['gswp3'][1],
               s=2*AREA_FACTOR * np.pi * 33**2,
               c='black', marker='$\mathrm{G3}$')
    ax.scatter(gris_results['gswp3'][0], gris_results['gswp3'][1],
               c='#1b3d6d', s=AREA_FACTOR * np.pi * 33**2, alpha=0.5,
               marker='o')
    ax.scatter(ais_results['gswp3'][0], ais_results['gswp3'][1],
               s=2*AREA_FACTOR * np.pi * 60**2,
               c='black', marker='$\mathrm{G3}$')
    ax.scatter(ais_results['gswp3'][0], ais_results['gswp3'][1],
               c='#f78d2d', s=AREA_FACTOR * np.pi * 60**2, alpha=0.5,
               marker='o')

    #plt.grid()
    #ax.set_thetamin(0)
    #ax.set_thetamax(90)

    plt.savefig(os.path.join('results', 'taylor_era5_wfde5_gswp3_sumup.pdf'))

def main():
    run()

if __name__=='__main__':
    main()