# Main python program to create time-series of Lake data
import numpy as np
import ECCO_functions as ecco

print
print 'Start of ECCO data creation software'
print
lake_data = 'Lakes/comsat_fetch.geojson'         # <-- Specify the path and file to read
lake_geometry, lake_id, lake_name,lake_path = ecco.Read_Lakes(lake_data)

# As an example of how the software works, pick a lake and run the code.
num = 6                                                           # Pick an integer number to represent a lake
ecco.Show_Lake_Info(num,lake_id,lake_name,lake_path)    # Show the info for a specified lake
ecco.Show_Lake_Plot(num,lake_path,lake_name,lake_id)    # Preview a plot of the selected lake
lake_obj,x1,x2,y1,y2 = ecco.Pth_Create(lake_path[num][0][0][:])     # Create the lake object
print 'Lake Area (units km^2): ', ecco.Poly_Area2D(ecco.EqArea(lake_obj.vertices))

setlat = 5; setlon = 5
testdata,lon_atts,lat_atts=ecco.Gen_Test_Dat(setlon,setlat,lake_obj)  # Simulated data of 5x5 pix over a lake
pix_wts = ecco.Pixel_Weights(lake_obj, testdata,lon_atts,lat_atts)    # Calculate weight of pixels
ecco.Preview_Weights(lake_obj,pix_wts,lon_atts,lat_atts)              # A test plot to show weights vs. lake

# THe above can easily be Parallelized

# It is at this point in the code where a weighted mean could be calculated for a specified lake
# A 2D dataset (lon,lat) of any variable, can be fed into the below function, and produce the weighted
# mean of all pixels the lake covers. I.e. For a 3D (lon,lat,time) dataset. The below must be run
# in a loop, for a subset of each days values:

eg_data = np.random.rand(setlon,setlat) 

print 'Example mean: ',ecco.Weighted_Mean(pix_wts,eg_data)

print 'End of program'