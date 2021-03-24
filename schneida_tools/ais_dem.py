#!/usr/bin/env python

from osgeo import gdal

import os
import warnings

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
EPSG_NSIDC = 3976

class AisDem(object):
    def __init__(self, filename=os.path.join('data_raw', 'krigged_dem_nsidc.bin')):
        self.dataset = gdal.Open(filename, gdal.GA_ReadOnly)
        if not self.dataset:
            warnings.warn('GDAL failed to open %s' % filename)
        
        #org_crs = CRS.from_epsg(EPSG_NSIDC)
        org_crs=CRS.from_proj4("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 "
                               "+k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m "
                               "+no_defs")
        self.spatial2coordinate(self.dataset.GetSpatialRef())
        self.transformer = Transformer.from_crs(org_crs, self.proj_crs)
        self.prepare_cartopy()
            
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

    def draw_contours(self, ax, filename=os.path.join('data_clean',
                                                      "krigged_dem_nsidc.nc"),
                      errormap=os.path.join('data_clean', "krigged_dem_errormap_nsidc.nc"),
                      downsample=10,
                      cmap=ListedColormap(cc.linear_grey_0_100_c0)):
        aisdem = Dataset(filename)
        aisdem_errormap = Dataset(errormap)
        (xx, yy) = np.meshgrid(aisdem.variables['x'][::downsample],
                               aisdem.variables['y'][::downsample])
        (lats, lons) = self.transformer.transform(xx, yy)
        print('Drawing AIS DEM contours...')
        ax.contour(lons, lats,
                   np.ma.masked_where(lats <= -80,
                   aisdem.variables['Band1'][::downsample,::downsample]) - 
                   np.ma.masked_where(lats <= -80,
                   aisdem_errormap.variables['Band1'][::downsample,::downsample]),
                   levels=np.arange(0, 4070, 500),
                   cmap=cmap,
                   linewidths=0.5,
                   linestyles='solid',
                   alpha=1.0,
                   vmin=-4000,vmax=4000,
                   label='AIS DEM$^1$',
                   transform=ccrs.PlateCarree())
                   
        aisdem.close()
        aisdem_errormap.close()
                   
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
    
    def setup_map(self):
        """ Draw elevation contours on the map background
        """
        
        #fig = plt.figure(figsize=[10, 5])
        #ax1 = fig.add_subplot(1, 1, 1, projection=ccrs.SouthPolarStereo())
        #ax1 = plt.axes(projection=ccrs.SouthPolarStereo(true_scale_latitude=-71))
        ax1 = plt.axes(projection=ccrs.LambertAzimuthalEqualArea(
                       central_longitude=0.0, central_latitude=-90))
        # Limit the map to -60 degrees latitude and below.
        ax1.set_extent([-180, 180, -90, -65], ccrs.PlateCarree())

        return ax1
    
def write_netcdf():
    ds = gdal.Translate('krigged_dem_errormap_nsidc.nc', 'krigged_dem_errormap_nsidc.bin', format='NetCDF')

def test():
    test = AisDem()
    test.print_dataset_info()
    #write_netcdf()
    ax = test.setup_map()
    test.draw_contours(ax)
    
    plt.show()

def run():
    test()

def main():
    run()

if __name__=='__main__':
    main()
