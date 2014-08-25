# Start of function file 

import numpy as np
import sys
import os
from netCDF4 import Dataset
from osgeo import ogr
import datetime as dt
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

def get_lakes():
	fn = 'Lakes/comsat_fetch.geojson'
	driver = ogr.GetDriverByName('GeoJSON')
	ds = driver.Open(fn)
	lyr = ds.GetLayer('OGRGeoJSON')
	dfn = lyr.GetLayerDefn()
	#fieldnames = [dfn.GetFieldDefn(i).GetName() for i in range(dfn.GetFieldCount())]
	return

## use comsat.Altitude for the lake elevation for now... it will be different 
## for the many other lakes. This is only a pracice set. 





def ncdump(nc_fid):
    '''
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.

    Parameters
    ----------
    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object

    Returns
    -------
    nc_attrs : list
        A Python list of the NetCDF file global attributes
    nc_dims : list
        A Python list of the NetCDF file dimensions
    nc_vars : list
        A Python list of the NetCDF file variables
    '''
    def print_ncattr(key):
        """
        Prints the NetCDF file attributes for a given key

        Parameters
        ----------
        key : unicode
            a valid netCDF4.Dataset.variables key
        """
        try:
            print "\t\ttype:", repr(nc_fid.variables[key].dtype)
            for ncattr in nc_fid.variables[key].ncattrs():
                print '\t\t%s:' % ncattr,\
                      repr(nc_fid.variables[key].getncattr(ncattr))
        except KeyError:
            print "\t\tWARNING: %s does not contain variable attributes" % key

    # NetCDF global attributes
    print "NetCDF Global Attributes:"
    nc_attrs = nc_fid.ncattrs()
    for nc_attr in nc_attrs:
        print '\t%s:' % nc_attr, repr(nc_fid.getncattr(nc_attr))
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    nc_vars = [var for var in nc_fid.variables]  # list of nc variables
    # Dimension shape information.
    print "NetCDF dimension information:"
    for dim in nc_dims:
        print "\tName:", dim 
        print "\t\tsize:", len(nc_fid.dimensions[dim])
        print_ncattr(dim)
    # Variable information.
    print "NetCDF variable information:"
    for var in nc_vars:
        if var not in nc_dims:
            print '\tName:', var
            print "\t\tdimensions:", nc_fid.variables[var].dimensions
            print "\t\tsize:", nc_fid.variables[var].size
            print_ncattr(var)
    return nc_attrs, nc_dims, nc_vars






# Read in some Land-surface data to use for test plotting
#def ls_msk():
#print 'Reading 1x1 degree ls-mask'
tmp = Dataset('Extra/landsea.nc','r')

nc_attrs, nc_dims, nc_vars = ncdump(tmp)


lsm = tmp.variables.get('LSMASK')
lon = tmp.variables.get('lon')
lat = tmp.variables.get('lat')

#lats = tmp.variables['lat'][:]
#lons = tmp.variables['lon'][:]	
#lsm = tmp.variables['LSMASK'][:]

plt.show()


	#def tst():
	#nc_dims = [dim for dim in tmp.dimensions]  # list of nc dimensions
    #nc_vars = [var for var in tmp.variables]  # list of nc variables
    # To print info about the unknown..
    #for dim in nc_dims:
    ##    print "\tName:", dim 
     #   print "\t\tsize:", len(tmp.dimensions[dim])
        #print ncattr(dim)
    # Variable information.
    #print "NetCDF variable information:"
    #for var in nc_vars:
    #    if var not in nc_dims:
    #        print '\tName:', var
    #        print "\t\tdimensions:", tmp.variables[var].dimensions
    #        print "\t\tsize:", tmp.variables[var].size
    #        #print_ncattr(var)
    #        '''
	#return lats