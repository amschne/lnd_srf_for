#!/usr/bin/env python

""" Read in Greenland Ice Sheet Mapping Project digital elevation model and plot
    contours.

    Data provided by the Greenland Ice Sheet Mapping Project (Howat, I.M., 
    Negrete, A., Smith, B.E., 2014).
"""

import warnings
from os import path

from osgeo import gdal
from pyproj.crs import CRS
from pyproj import Transformer
import cartopy.crs as ccrs
from netCDF4 import Dataset
import numpy as np

from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
import colorcet as cc

WSG_84 = 4326
EPSG_NSIDC = 3413

class GrisDEM(object):
    def __init__(self, filename=path.join('data_raw','gimpdem_90m_v01.1.tif'),
                 wide=False):
        self.dataset = gdal.Open(filename, gdal.GA_ReadOnly)
        if not self.dataset:
            warnings.warn('GDAL failed to open %s' % filename)
        self.wide = wide
        
        org_crs = CRS.from_epsg(EPSG_NSIDC)
        self.spatial2coordinate(self.dataset.GetSpatialRef())
        self.transformer = Transformer.from_crs(org_crs, self.proj_crs)
        self.prepare_cartopy()
                
    def draw_contours(self, ax, filename=path.join('data_clean',
                                                   "gimpdem_90m_v01.1.nc"),
                      cmap=ListedColormap(cc.linear_grey_0_100_c0)):
        """
        """
        downsample=np.gcd(self.dataset.RasterXSize, self.dataset.RasterYSize)
        gimpdem = Dataset(filename)
        (xx, yy) = np.meshgrid(gimpdem.variables['x'][::downsample],
                               gimpdem.variables['y'][::downsample])
        (lats, lons) = self.transformer.transform(xx, yy)
        print('Drawing GIMP DEM contours...')
        ax.contour(lons, lats,
                   gimpdem.variables['Band1'][::downsample,::downsample],
                   levels=np.arange(0, 3207, 500),
                   cmap=cmap,
                   linewidths=0.5,
                   linestyles='solid',
                   alpha=1.0,
                   vmin=-4000,vmax=4000,
                   label='GIMP DEM',
                   transform=ccrs.PlateCarree())
                   
        gimpdem.close()
                   
        return ax
    
    def spatial2coordinate(self, osr_crs):
        """ Convert from osgeo.osr.SpatialReference to pyproj.crs.CRS
    
            http://pyproj4.github.io/pyproj/stable/crs_compatibility.html#
        """
        osr_crs.ImportFromEPSG(WSG_84)
        self.proj_crs = CRS.from_wkt(osr_crs.ExportToWkt())
    
    def prepare_cartopy(self):
        """ Prepare pyproj.crs.CRS for cartopy.crs.CRS
    
            http://pyproj4.github.io/pyproj/stable/crs_compatibility.html#cartopy
        """
    
        globe = ccrs.Globe(ellipse=None,
                           semimajor_axis=self.proj_crs.ellipsoid.semi_major_metre,
                           semiminor_axis=self.proj_crs.ellipsoid.semi_minor_metre,
                           inverse_flattening=self.proj_crs.ellipsoid.inverse_flattening)
        proj_dict = self.proj_crs.to_dict()
        proj_dict["pm"] = self.proj_crs.prime_meridian.longitude
        try:
            self.cart_crs = ccrs.CRS(proj_dict, globe=globe)
            print("Returning cartopy CRS.")
        except AttributeError:
            warnings.warn('Can not initialize cartopy CRS. Returning "globe" '
            'instead.')
            self.globe = globe
        
    def setup_map(self, central_longitude=-42.1, central_latitude=71.4,
                  wide=False):
        """ Setup the map projection
        """
        if True:
            ax = plt.axes(projection=ccrs.LambertAzimuthalEqualArea(
                            central_longitude=central_longitude,
                            central_latitude=central_latitude))
        else:    
            ax = plt.axes(projection=ccrs.AlbersEqualArea(central_longitude=
                                                      central_longitude, 
                                                      central_latitude=71.4,
                                                      standard_parallels=(60, 83), 
                                                      globe=None))
        if self.wide:
            ax.set_extent((-55, -29, 54, 84))
        elif False:    
            ax.set_extent((-61, -29, 59, 84))
        else:
            ax.set_extent((-55, -29, 59, 84))
            
        return ax
    
    def print_dataset_info(self):
        print("Driver: {}/{}".format(self.dataset.GetDriver().ShortName,
                                     self.dataset.GetDriver().LongName))
        print("Size is {} x {} x {}".format(self.dataset.RasterXSize,
                                            self.dataset.RasterYSize,
                                            self.dataset.RasterCount))
        print("Projection is {}".format(self.dataset.GetProjection()))
        geotransform = self.dataset.GetGeoTransform()
        if geotransform:
            print("Origin = ({}, {})".format(geotransform[0], geotransform[3]))
            print("Pixel Size = ({}, {})".format(geotransform[1], geotransform[5]))
    
        band = self.dataset.GetRasterBand(1)
        print("Band Type={}".format(gdal.GetDataTypeName(band.DataType)))

        min = band.GetMinimum()
        max = band.GetMaximum()
        if not min or not max:
            (min,max) = band.ComputeRasterMinMax(True)
        print("Min={:.3f}, Max={:.3f}".format(min,max))

        if band.GetOverviewCount() > 0:
            print("Band has {} overviews".format(band.GetOverviewCount()))

        if band.GetRasterColorTable():
            print("Band has a color table with {} "
                  "entries".format(band.GetRasterColorTable().GetCount()))
    
def test():
    #plt.style.use('agu_half_vertical')
    plt.style.use('gmd_movie_frame')
    test = GrisDEM()
    test.print_dataset_info()
    ax = test.setup_map()
    test.draw_contours(ax)
    
    plt.show()

def run():
    test()

def main():
    run()

if __name__=='__main__':
    main()
