#!/usr/bin/env python

""" Defining a coordinate space, from "Image Processing and GIS for Remote
    Sensing: Techniques and Applications" Second Edition. Jian Guo Liu and
    Philippa J. Mason.

    Coded by Adam M. Schneider (amschne@umich.edu)
"""

from os import path

from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import colorcet as cc

class CoordinateSpace(object):
    def __init__(self, spheroid='WGS84'):
        if spheroid == 'WGS84':
            self.wgs84()
        elif spheroid == 'MERIT':
            self.merit()
        
        self.get_transform()
    
    def wgs84(self):
        self.proj_crs = CRS.from_epsg(4326)
        self.spheroid()
    
    def merit(self):
        self.proj_crs = CRS.from_proj4('+proj=stere +lat_0=90 +ellps=MERIT')
        self.spheroid()
    
    def spheroid(self):
        self.globe = ccrs.Globe(ellipse=None,
                        semimajor_axis=self.proj_crs.ellipsoid.semi_major_metre,
                        semiminor_axis=self.proj_crs.ellipsoid.semi_minor_metre,
                  inverse_flattening=self.proj_crs.ellipsoid.inverse_flattening)
                           
        self.proj_dict = self.proj_crs.to_dict()
        self.proj_dict["pm"] = self.proj_crs.prime_meridian.longitude
        
    def get_transform(self):
        proj = Transformer.from_crs(self.proj_crs.geodetic_crs, self.proj_crs)
        self.transform = proj.transform

def nh_horizontal_comparison(map_lat_min=15, map_lon_0=-80,
                             lon_lines = np.arange(-180, 180, 30),
                             lat_lines = np.arange(-90, 90, 30)):
    """
    """
    map_proj = ccrs.LambertAzimuthalEqualArea(central_longitude=map_lon_0,
                                              central_latitude=90.0,
                                              false_easting=0.0,
                                              false_northing=0.0)
                             
    plt.style.use(path.join('schneida_tools', 'gmd_movie_frame.mplstyle'))
    
    nrows=1
    ncols=4
    geo_axes1 = plt.subplot(nrows, ncols, 1, projection=map_proj)
    geo_axes2 = plt.subplot(nrows, ncols, 2, projection=map_proj)
    geo_axes3 = plt.subplot(nrows, ncols, 3, projection=map_proj)
    geo_axes4 = plt.subplot(nrows, ncols, 4, projection=map_proj)
    axes = [geo_axes1, geo_axes2, geo_axes3, geo_axes4]
    for i, ax in enumerate(axes):
        ax.set_extent((-180, 180, map_lat_min, 90), crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.LAND, color='#555759')
        ax.add_feature(cfeature.OCEAN, color='#555759')
        
    return axes
        
def draw_elevation_contours(ax, elevation_levels=np.arange(0,8849,1500),
                            vmin=-7500, vmax=7500):
    """ Add evelvation contours using WFDE5 data
    """
    rootgrp = Dataset(path.join('data_raw', 'ASurf_WFDE5_CRU_v1.1.nc'))
    ax.contour(rootgrp.variables['lon'][:], rootgrp.variables['lat'][:],
               np.ma.filled(rootgrp.variables['ASurf'][:], fill_value=0),
               levels=elevation_levels, vmin=vmin, vmax=vmax,
               linewidths=0.5,
               cmap='cet_CET_L1', 
               transform=ccrs.PlateCarree())

def test():
    axes = nh_horizontal_comparison()
    for i, ax in enumerate(axes):
        draw_elevation_contours(ax)
        
    plt.show()

def run():
    test()
    
def main():
    run()

if __name__=='__main__':
    main()