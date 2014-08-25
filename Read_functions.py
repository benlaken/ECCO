# Start of function file 

import numpy as np
from osgeo import ogr

fn = 'comsat.geojson'
driver = ogr.GetDriverByName('GeoJSON')
ds = driver.Open(fn)
lyr = ds.GetLayer('OGRGeoJSON')
dfn = lyr.GetLayerDefn()
fieldnames = [dfn.GetFieldDefn(i).GetName() 
              for i in range(dfn.GetFieldCount())]

## use comsat.Altitude for the lake elevation for now... it will be different 
## for the many other lakes. This is only a pracice set. 
