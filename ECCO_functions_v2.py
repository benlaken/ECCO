import numpy as np
import pandas as pd
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import mpl_toolkits.basemap.pyproj as pyproj

import glob
import os
import time
import h5py
import time as clock
import osgeo.ogr
import subprocess
import sys
import json
import simplejson
import csv

from netCDF4 import Dataset, num2date, date2num
from matplotlib import cm, path
from matplotlib.path import Path
from matplotlib.transforms import Bbox
from math import pi, cos, sin, radians, atan, asin


def Area_Lake_and_Islands(lake_poly):
    '''Purpose  -  Calculate the area of a given Lake polygon, taking into account islands.
    Input - Lake polygon object (which may or may not include islands). Polygon object 
            should be generated with the function Path_Lake_and_Islands().
    Output- If islands are present the output is the area (in km^2) of the lake boundary
            minus the sum of the island(s) area. If no islands are present then the output
            is simply the area of the lake boundary.
    '''
    aaa = np.where(lake_poly.codes == 1)
    bbb = np.where(lake_poly.codes == 79)
    #print 'number of 1s:',len(aaa[0])
    #print 'number of 79s',len(bbb[0])
    island_num = (len(aaa[0]))-1
    #print 'Found %i Islands'%island_num
    area_start = [i for n,i in enumerate(aaa[0])]
    area_end = [i for n,i in enumerate(bbb[0])]
    #Seperated paths, pos.0 is the boundary and rest are islands
    sep_paths = [lake_poly.vertices[area_start[n]:area_end[n]] for n,i in enumerate(area_start)]
    sep_areas = [float(Poly_Area2D(EqArea(sep_paths[n]))) for n in xrange(len(sep_paths))]
    if island_num > 0:
        total = sep_areas[0] - sum(sep_areas[1:])
        #print 'Area minus Islands:',total
    else:
        total = sep_areas[0]
        #print 'No islands found, total area is simply lake boundary'
    #print 'Area without considering islands:',sep_areas[0]
    return total

def BBOX_PathCode_Fix(lake_poly):
    '''
    Need to solve a problem where BBOX trimmed paths have no end 79 code. This
    breaks my super-area calculation code for islands...
    '''
    #print 'before hack',lake_poly.codes 
    if 79 in lake_poly.codes:  # Gives false if no 79 code exists (no end polypath)
        #print 'True 79 code exists, taking no action'
        x = 0  # doesnt do anything, just lets the code pass to the else
    else:
        #print 'False 79 code is not present, adding it to correct places'
        aaa = np.where(lake_poly.codes == 1)
        #print len(aaa[0]),'start codes'
        # Need to insert 79 codes at: 1) before all 1 codes after position 0, and 
        # 2) in the last element of the array
        lake_poly.codes[-1] = 79                 # Start with the last element...
        if(len(aaa[0]) > 1):   # for sitatuion with at least one island...
            for n,i in enumerate(aaa[0]):
                if n > 0:
                    lake_poly.codes[i-1]=79
    return lake_poly



def Calc_Coordinates(lon_2transform,lat_2transform):
    ''' Returns lat lon coordinates on a polar rotated sphere, from
    the input of the North Pole longitude, and North Pole latitude
    (that is the rotated position of the pole), and the lon and lat
    which you wish to transform (as two speperate floating point
    values).
    Note   - Currently this has the CORDEX EUR-11 pole shift
    hardcoded into the routine, as lo_polo = 198. (Cartesian Lontidue
    of N. Pole Shift), and la_polo = 39.25 (Latitude of N.Pole shift)
    '''
    lo = lon_2transform
    la = lat_2transform
    lo_polo=198.                  # Lon and lat of the new position of the north pole
    la_polo=39.25
    lon=lo*(pi/180)               # Transform into radians
    lat=la*(pi/180)
    lon_polo=lo_polo*pi/180       # Transformation into radians
    lat_polo=la_polo*pi/180
    phi=pi-lon_polo               # Calcuus of the angles we are rotating to move
    teta=-(pi/2-lat_polo)         # from the real north pole to the new one
    x=cos(lon)*cos(lat)           # Change in coordinates from lon, lat to cardinates
    y=sin(lon)*cos(lat)
    z=sin(lat)
    xr=cos(teta)*cos(phi)*x-cos(teta)*sin(phi)*y-sin(teta)*z    # Calculus of the new rotated cordinates in cartesians
    yr=sin(phi)*x+cos(phi)*y
    zr=sin(teta)*cos(phi)*x-sin(teta)*sin(phi)*y+cos(teta)*z
    lonr=atan(yr/xr)               # Transformation from cartesians into lon and lat again 
    latr=asin(zr)
    return lonr*180/pi, latr*180/pi


def catchment_timeseries(nc_path, outputprefix,plots = False,rprt=False,sbar=False):
    '''
    Function to produce time series over catchment areas for a given input file of
    CORDEX data. Metadata and weights have been pre-calculated with the function
    Catchment_Weights_Meta().
    Output prefix is the file path where you want the output data to go e.g.
    outputprefix='Folder1/Subfolder/'
    Note, this function requires metadata (csv) and weights (hdf5) located in a
    relative path of Catchments/Metadata/ and Catchments/Weights respectivley.
    These files can be obtained from http://www.files.benlaken.com/documents/
    and are Catchment_meta.csv and catchment_weights.h5.
    '''
    # Section for loading of data and opening of an output file
    if rprt:
        atime = time.time()
    
    clim_dat,rlat,rlon,timeCDX,metadata,txtfname = ECCO.Read_CORDEX_V2(nc_path) # CORDEX NetCDF Read file
    vname, m1, m2, dexp, m3, m4, m5, m6, drange_orignial = metadata        # Metadata of fname string
    var_type = clim_dat.standard_name                                      # What kind of CORDEX data?
    dat_loaded = clim_dat[:,:,:]                                           # Load CORDEX data into RAM
    rlat_loaded = rlat[:]
    rlon_loaded = rlon[:]
    
    thefilename = 'Catchment_'+str.split(nc_path,'/')[-1][:-3] 
    FILE= outputprefix + thefilename +'.h5'                # Set up HDF5 file output
    if os.path.isfile(FILE) == True:
        print 'File already exists: Overwriting...'
        os.remove(FILE)
    else:
        print 'Creating file: ',FILE
    f = h5py.File(FILE,'w')
    
    # Metadata contains EB_id(index), area, npix, xpix, ypix
    # Where the number of pixels is 1, we can use the xpix, ypix directly as the time series
    # Where the npix is larger than 1, we can call the precalculated array of weights.
    catch_meta = pd.read_csv('Catchments/Metadata/Catchment_meta.csv') # Metadata
    catch_meta = catch_meta.set_index('EB_id')              # Hex code is used as id
    
    lake_file='Catchments/ecco_biwa_catchments_part_1.shp'  # just for testing...
    ShapeData = osgeo.ogr.Open(lake_file)                   # Connection to catchment shapes
    TheLayer = ShapeData.GetLayer(iLayer=0)
        
    if rprt:
        btime = time.time()
    if sbar:
        icnt = 0
    tot_shapes = 274791  # total number of shapes in the three shape files
    
    # For the three catchment files - loop over each file and do every feature element
    catchflist = []                               # Gather precalculated surface weights 
    for fnm in glob.glob("Catchments/*.shp"):     # N.b. You can precalculate as many as you
        catchflist.append(fnm)
        #print fnm
        ShapeData = osgeo.ogr.Open(fnm)           # Make a link to Lake Shape Files
        TheLayer = ShapeData.GetLayer(iLayer=0)
        dolakes=range(TheLayer.GetFeatureCount()) # Create a range to loop over lake features  
   
        for n in dolakes:
            tlist = []                        # Will hold the time series for each n of dolake
            feature1 = TheLayer.GetFeature(n)                        # Get catchtment data
            lake_feature = feature1.ExportToJson(as_object=True)     # Convert to JSON
            EB_id = hex(int(lake_feature['properties']['ebint']))[2:]# Extract id number
            wgs85_xy = Reproj_Catchment(lake_feature=lake_feature,
                                       chatty=False)                  # Convert to WGS system
            lake_cart = ECCO.Path_LkIsl_ShpFile([wgs85_xy])          # Create shape object
            lake_rprj = ECCO.Path_Reproj(lake_cart,False)            # Reproj 2 Plr rotated
                    
            if plots:     
                ECCO.Preview_Lake(lake_cart)        
                print 'Area in km^2 (not inc. islands):',ECCO.Area_Lake_and_Islands(lake_poly=lake_cart)         
                print ', No. xy bound. points:',len(lake_cart.vertices)
            
            if catch_meta.npix[EB_id] == 1:
                #print 'One pixel check:',EB_id,catch_meta.npix[EB_id]
                ypix = catch_meta.ypix[EB_id]       # Get the pre-calc. pixel indexes...
                xpix = catch_meta.xpix[EB_id]       # ...calc in MT_Gen_SWeights() earlier
                tlist = dat_loaded[:, ypix, xpix]
                #plt.plot(tlist-273.15,alpha=0.4)
                #plt.show()
            elif catch_meta.npix[EB_id] > 1:
                weight_mask = get_catchweight(tmp_ebid=EB_id)  # Get the catchment weights from precalc file
                # If you were to calculate the weighted mask again you can do it from this...
                # weight_mask = ECCO.Pixel_Weights(lake_in=lake_rprj,datin=sub_clim,
                #                                 lat_atts=sub_rlat,lon_atts=sub_rlon)
                
                #plt.imshow(weight_mask,interpolation='none',cmap=plt.cm.gray,origin='lower')
            
                sub_clim,sub_rlat,sub_rlon = ECCO.TrimToLake3D(lake_rprj,dat_loaded,rlat_loaded,rlon_loaded,
                                                      off = 3, show = False)
                
                #print 'More pixels check:',EB_id,catch_meta.npix[EB_id]
                tlist = ECCO.Weighted_Mean_3D(weight_mask=weight_mask,all_time_clim=sub_clim,chatty=False)
                #tlist2 = ECCO.Weighted_Mean_3D_old(weight_mask=weight_mask,all_time_clim=sub_clim,chatty=False)       
                #plt.plot(tlist-273.15,tlist2-273.15alpha=0.4)  # shows old vs new weighting func to test if 
                #plt.show() # that was a potential error source. Seems like it is fine.          
                
            if plots:
                ECCO.Show_LakeAndData(lake_rprj,dat_loaded[0,:,:],rlat,rlon,zoom=6.)
                ECCO.Preview_Weights(lake_rprj,weight_mask,sub_rlat,sub_rlon) 
            
            if sbar:
                icnt=icnt+1
                if (float(icnt) % 10.) == 0.0:
                    ECCO.Update_Progress(float(icnt)/(tot_shapes-1.0))
                    
            write_hdf_catchment(fw=f,EB_id=EB_id,tseries=tlist)  #Write timeseries data to HDF file
            # Loop of n (all features)
        # Loop of 3 shape files (containing features)
    # Outside of all loops
    f.close()   # Close the HDF5 file after looped through all 3 shape files
    subprocess.call(["gzip", FILE])     # Compress the file and remove original with gzip             
    #-------------------------------------------------------------------------------------------
    if rprt:      # Finish and report time if requested
        print '\nTime to read data: %4.2f sec'%(btime - atime)
        print 'Time to Process %i shapes: %4.2f sec'%(tot_shapes,time.time() - btime)
    return


def Catchment_Weights_Meta(nc_path,sbar=False):
    '''
    From catchment data, generate surface weights if requested, and meta
    data files also.
    '''
    # 1. LOAD Climate DATA
    clim_dat,rlat,rlon,timeCDX,metadata,txtfname = Read_CORDEX_V2(nc_path) # CORDEX NetCDF Read file
    vname, m1, m2, dexp, m3, m4, m5, m6, drange_orignial = metadata        # Metadata of fname string
    var_type = clim_dat.standard_name                                      # What kind of CORDEX data?
    dat_loaded = clim_dat[:,:,:]                                           # Load CORDEX data into RAM
    rlat_loaded = rlat[:]
    rlon_loaded = rlon[:]
    # Create a hdf5 file of catchment weights
    thefilename = 'catchment_weights'
    FILE= 'Catchments/Weights/' + thefilename +'.h5'                # Set up HDF5 file output
    if os.path.isfile(FILE):
        print 'HDF5 File already exists. Leaving loop so you dont clobber it by accident.'
        print 'To run this function, decide manually if you want to remove it or not.'
        return
        #print 'hdf weights file exists, removing it...'
        #os.remove(FILE)
    else:
        print 'Creating file: ',FILE
        fweights = h5py.File(FILE,'w')
    
    # set and write header info for the metadata file
    metacsv = 'Catchments/Metadata/Catchment_meta.csv'
    if os.path.isfile(metacsv) == True:
        print 'Earlier metadata exists. Erasing it...'
        os.remove(metacsv)
    tmplist = ['EB_id','area','npix','ypix','xpix'] 
    write_metadata_csv(mfname=metacsv,meta_list=tmplist)
    
    if sbar:
        icnt = 0
    
    # For the three catchment files - loop over each file and do every feature element...
    catchflist = []                                     # Gather precalculated surface weights 
    for fnm in glob.glob("Catchments/*.shp"):           # N.b. You can precalculate as many as you
        catchflist.append(fnm)
        ShapeData = osgeo.ogr.Open(fnm)                  # Make a link to Lake Shape Files
        TheLayer = ShapeData.GetLayer(iLayer=0)
        dolakes=range(TheLayer.GetFeatureCount())   # Create a range to loop over lake features   
        for n in dolakes:
            tlist = []
            feature1 = TheLayer.GetFeature(n)                        # Get catchtment data
            lake_feature = feature1.ExportToJson(as_object=True)     # Convert to JSON
            EB_id = hex(int(lake_feature['properties']['ebint']))[2:]# Extract id number
            wgs85_xy = Reproj_Catchment(lake_feature=lake_feature,
                                       chatty=False)                 # Convert to WGS system
            lake_cart = Path_LkIsl_ShpFile([wgs85_xy])          # Create shape object
            lake_rprj = Path_Reproj(lake_cart,False)            # Reproj 2 Plr rotated
            sub_clim,sub_rlat,sub_rlon = TrimToLake(lake_in=lake_rprj,Cdat=dat_loaded[0,:,:],
                                                         rlat=rlat_loaded,rlon=rlon_loaded,
                                                         off = 3, show = False) 
            weight_mask = Pixel_Weights(lake_in=lake_rprj,datin=sub_clim,
                                             lat_atts=sub_rlat,lon_atts=sub_rlon)
            # For making the surface-weights saved .npy arrays and metadata (pixel counts and xy)
            pix_truth = (weight_mask > 0.0)    # Count how many times the weight mask is
            pxnum = len(weight_mask[pix_truth])  #  above 0.0 (i.e. how many pixels of data are needed)
        
            if pxnum > 1:
                #np.save('Catchments/Weights/'+EB_id,weight_mask)
                Write_HDF_weights(fw=fweights,EB_id=EB_id,weights=weight_mask) 
                ypix = -99
                xpix = -99

            if pxnum < 1:
                pxnum = 1                                               # Small bug fix, no biggy...
            if pxnum == 1:
                xxx,yyy = Get_LatLonLim(xypath=lake_rprj.vertices)  # Find upp./low.lake lims.
                ypix = (Closest(array=rlat,value=yyy[0]))                # For lakes of one pixel  
                xpix = (Closest(array=rlon,value=xxx[0]))
            tmplist=[EB_id, Area_Lake_and_Islands(lake_poly=lake_cart),pxnum,ypix,xpix]
            write_metadata_csv(mfname=metacsv,meta_list=tmplist)
            if sbar:
                icnt=icnt+1
                if (float(icnt) % 10.) == 0.0:
                    Update_Progress(float(icnt)/274791.)
    fweights.close()
    return


def Closest(array, value):
    '''Purpose    -  Functions like the IDL routine CLOSEST.
    Essentially just returns the value of an array closest
    to a specified value.
    Input    -  array: a np.array
             -  value: the value to which you are looking
                for the closest match to
    Output   - Returns an integer value index to the input array, of
               the data point
               most closely corresponding to the input value.
    '''
    y = [0] * len(array)            # Declare a list to hold (array - value) numbers.
    #out = [0] * 2                  # A list to hold the output.
    for i in xrange(len(array)):
        y[i] = np.abs(array[i] - value)
    mval = np.min(y)
    #out[1] = np.where(y == mval)   # Identify where smallest diffrence occurs in y.
    #out[0] = array[out[1]]
    out = np.where(y == mval)
    return int(out[0])


def EqArea(verts):
    '''Purpose - Take gridded data and, using the assumption of a 
    spherical earth, re-project it to a spherical coordinate system.
    Input   - The Matplolib.Path.Path object vertexes
    Output  - x,y coordinates projected onto a sphere
    Notes   - Solution from stackoverflow.com/questions/4681737/how
              -to-calculate-the-area-of-a-polygon-
              on-the-earths-surface-using-python
    By multiplying the latitude by the length of one degree of 
    latitude, and the longitude by the length of
    a degree of latitude and the cosine of the latitude. Then, 
    calculate the area of an arbitrary polygon in a plane.
'''
    earth_radius = 6367.4447 # Earth avg. radius (km) Wolfram Alpha
    lat_dist = pi * earth_radius / 180.0  
    eqout =[]
    for n,i in enumerate(verts):
        longitude = i[0] 
        latitude = i[1]
        x = longitude * lat_dist * cos(radians(latitude))
        y = latitude * lat_dist
        eqout.append([x,y])
    return eqout


def Fast_v3(nc_path, lake_file, outputprefix,lstart=0,lstop=275265,
                       hexlist=None,tt=None,plots = False,rprt=False,sbar=False,
                       rprt_loop=False):
    '''
    Input:
    
           hexlist      If not None, this must be a list of hexcodes which match lake codes.
                        These codes will be the list of lakes to be processed, with
                        (regardless of lstart / lstop settings). The list should be ascii values
                        in any order. E.g. hexlist = ['a2204','155980','d23e4a','7aa917']
    
     k.w.agrs:   
    
           plots        True or False(default). If True, preview plots are created. I reccomend
                        using this feature carefully. Unless you have set a very small number of
                        lakes, this will produce a huge number of plots!
    
           rprt         Returns information on how long the program took to load the data,
                        complete, and how many lakes were processed.
        
           rprt_loop    Returns info on how long each specific lake took (can be a lot...)
    
           sbar         Create a status bar, to show the progression through a large loop. Not
                        shown by default. As, when this pro is on MPI it is not a good feature.
    
     Notes:       Out of 275265 total lakes, 264532 of them are within one pixel of EUR-11 data.
                  Metadata called from lake hexcode as index: not shapefile feature number.
    '''
    # 1. LOADING DATA SECTION
    if rprt:
        atime = clock.time()

    lk_processed_inf = pd.read_csv('Metadata/Meta_Lakes.csv')  # Pre-processed lake metadata (CSV)
    lk_processed_inf.index = lk_processed_inf.hex           # Use the hex-code column as the index 
    
    ShapeData = osgeo.ogr.Open(lake_file)                  # Make a link to Lake Shape Files
    TheLayer = ShapeData.GetLayer(iLayer=0)
    
    clim_dat,rlat,rlon,timeCDX,metadata,txtfname = Read_CORDEX_V2(nc_path) # CORDEX NetCDF Read file
    vname, m1, m2, dexp, m3, m4, m5, m6, drange_orignial = metadata     # Metadata of fname string
    var_type = clim_dat.standard_name                                   # What kind of CORDEX data?
    dat_loaded = clim_dat[:,:,:]                                        # Load CORDEX data into RAM
    rlat_loaded = rlat[:]
    rlon_loaded = rlon[:]
    
    orog = Height_CORDEX()                                 # NetCDF EUR-11 surface height data 
    
    precalculated = []                                     # Gather precalculated surface weights 
    for fnm in glob.glob("Metadata/Weights/*.npy"):           # N.b. You can precalculate as many as you
        precalculated.append(fnm[14:-4])                   # like: place in folder to run (for speed)
    precalculated = np.array(precalculated)                # Make it a np.array (needed for functions)
    
    if hexlist == None:                                    # Set up the list of lakes to process:
        dolakes=np.arange(lstart,lstop,1)          #If no Hexcodes, use lstart/lstop to form a list
    else:
        dolakes= lk_processed_inf.num[test]     #If hexcodes, then gen. list of nums from PD object
    
    thefilename = 'Lakes_'+str.split(nc_path,'/')[-1][:-3]
    # Set up HDF5 file output
    FILE= outputprefix + thefilename +'.h5'
    #print FILE
    if os.path.isfile(FILE) == True:
        print 'Earlier file already exists: Overwriting...'
        os.remove(FILE)
    else:
        print 'No file found. Creating: ',FILE
    f = h5py.File(FILE,'w')

    # 2. LOOP OVER ALL LAKES (or specified lakes from lstart to lstop)
    if rprt:
        btime = clock.time()
    if sbar:
        icnt = 0
    for n in dolakes:
        tlist = []
        #if rprt_loop == True:
        #    ltime = clock.time()
        feature1 = TheLayer.GetFeature(n)           # Get individ. lake in shapefile
        lake_feature = feature1.ExportToJson(as_object=True)
        lake_cart = Path_LkIsl_ShpFile(lake_feature['geometry']['coordinates'])
        lake_altitude=lake_feature['properties']['stf_mean']
        EB_id = lake_feature['properties']['EBhex']
        EB_id = EB_id[2:]                           # Strip off the hexcode label 0x
        lake_rprj = Path_Reproj(lake_cart,False)    # Reproj. lake to CORDEX plr. rotated
        #if rprt_loop == True:
        #    print '\rCheck:',n,EB_id,lk_processed_inf.hex[EB_id],lk_processed_inf.npix[EB_id],
        if EB_id != lk_processed_inf.index[n]:      # Some handy error check
            print 'Warning! Lake feature and metadata miss-match for some reason. Check it out:'
            print 'Problem at:',num,lk_processed_inf.num[n],EB_id[2:],lk_processed_inf.index[n]
        if plots:     
            Preview_Lake(lake_cart)        
            print 'Area in km^2 (not inc. islands):', Area_Lake_and_Islands(lake_cart),         
            print ', No. xy bound. points:',len(lake_cart.vertices)
        if lk_processed_inf.npix[EB_id] == 1:         # ONE PIXEL LAKES <<<
            ypix = lk_processed_inf.ypix[EB_id]       # Get the pre-calc. pixel indexes...
            xpix = lk_processed_inf.xpix[EB_id]       # ...calc in MT_Gen_SWeights() earlier
            if lake_altitude == None:                 # Some lakes don't have alitude values
                offset = -999.
            else:
                offset = OnePix_HOffset(lake_altitude,orog[ypix, xpix],var_type)
            tlist = dat_loaded[:, ypix, xpix]
            #if rprt_loop == True:
            #    print '1pix, only slicing. Time:',clock.time() - ltime
        else:                                         # LAKES OF MORE THAN ONE PIXEL <<<
            pre_test = (lk_processed_inf.hex[EB_id] == precalculated)
            if(any(pre_test) == True):                # Scipy's any() evalautes list truth
                weightfile = 'Metadata/Weights/'+precalculated[pre_test][0]+'.npy'
                weight_mask = np.load(weightfile)
            else:  # If no pre-calculated weight mask file then calculate it now
                sub_clim,sub_rlat,sub_rlon = TrimToLake(lake_rprj,dat_loaded[0,:,:],rlat_loaded,
                                                        rlon_loaded,off = 3, show = False) 
                weight_mask = Pixel_Weights(lake_rprj,sub_clim,sub_rlat,sub_rlon)
            if ((var_type == 'air_temperature')| (var_type == 'surface_air_pressure')): 
                sub_orog,sub_rlat,sub_rlon = TrimToLake(lake_rprj,orog,rlat_loaded,
                                                            rlon_loaded,off = 3, show = False)
                #print 'Stats 2:',n,EB_id,lake_altitude
                if lake_altitude == None:                 # Some lakes don't have alitude values
                    offset = -999.
                else:
                    hght,offset = Orographic_Adjustment(weight_mask,sub_orog,
                                                        lake_altitude,clim_dat,chatty=False)
            else:
                hght = -999.                         # If no offset calculated then
                offset = -999.                       # just set them to missing data
            
            sub_clim,sub_rlat,sub_rlon = TrimToLake3D(lake_rprj,dat_loaded,rlat_loaded,rlon_loaded,
                                                      off = 3, show = False)
            tlist = Weighted_Mean_3D(weight_mask, sub_clim, chatty=False)  # Here's the t-series
            tlist = np.squeeze(tlist)                                      # Remove empty dimension
            
            if plots:
                Show_LakeAndData(lake_rprj,dat_loaded[0,:,:],rlat,rlon,zoom=6.)
                Preview_Weights(lake_rprj,weight_mask,sub_rlat,sub_rlon) 
            
            #if rprt_loop == True:
            #    print '\r>2pix, weighting needed. Time:',clock.time() - ltime,
                
        if rprt_loop:
            print '\rStats:',(float(n)/float(lstop))*100.,'% ',n,EB_id,offset,lake_altitude,
        if sbar:
            icnt=icnt+1
            if (float(icnt) % 10.) == 0.0:
                Update_Progress(float(icnt)/float(len(dolakes)-1))

        Write_HDF(f,EB_id,tlist,offset,lk_processed_inf.area[n])  # Write inside function
        feature1=0
        lake_feature = 0

    f.close()                                # Close the HDF5 file after the lake loop finishes
    subprocess.call(["gzip", FILE])             # Compress the file and remove original with gzip             
    if rprt:
        ctime = clock.time()

    if rprt:
            print '\nTime to read data: %4.2f sec'%(btime - atime)
            print 'Time to Process %i lakes: %4.2f sec'%(len(dolakes),ctime - btime)
    return



def Folder_Create(outputprefix,fnm_head,EBid):
    '''Purpose   -  A function to generate the file paths for a 
    given lake (from the EB-ID and user-specified output path).
    '''
    out = os.path.join(outputprefix, fnm_head)       
    ideb = EBid[2:]
    #print(ideb)
    l = len(ideb)
    if l == 6:
        idnew = ideb
    elif l == 5:
        idnew = '0%s' % ideb
    elif l == 4:
        idnew = '00%s' % ideb
    elif l == 3:
        idnew = '000%s' % ideb
    elif l == 2:
        idnew = '0000%s' % ideb
    elif l == 1:
        idnew = '00000%s' % ideb
    return


def Get_LatLonLim(xypath):
    '''Purpose - Find the bounding (i.e. max/min) lat/lon of a lake
    Input   - Path.verticies of Lake
    Output  - x[max,min],y[max,min]
    '''
    xset = []
    yset = []
    for i in xypath:
        xset.append(i[0])
        yset.append(i[1])        
    return [max(xset),min(xset)],[max(yset),min(yset)]


def Gen_FileName(pth,lname,drange,ext):
    '''Purpose   - This is needed to use the unicode lake names to 
    generate a string used as a file name and folder path to save 
    output files. Idea is to modify this as needed.
    
    Inputs    - lname: Lake_name[n] unicode name
              - pth: path of where the file should go
              - ext: extention type of the desired file
    Outputs   - filename: A string with the path, filename, and type 
                all together. 
    ''' 
    #pth = '/uio/kant/geo-metos-u1/blaken/Work/Python/ECCO/Outputs/'
    #ext = '.pdf'
    filename = u''.join((pth,lname,'_',drange,ext)).encode('utf-8').strip()  
    return filename


def Gen_FName_Head(long_name):
    '''Purpose   - Set the string to start the file names with 
    (these strings should be the variable)
    Input    - clim_dat.long_name (from the NetCDF data)
    Output   - fnm_head: a short string to indicate the variable 
    '''
    fnm_head=''
    if long_name == 'Near-Surface Air Temperature':
        fnm_head = 'TAS_'
    return fnm_head

def get_catchweight(tmp_ebid):
    '''
    From the hex EB_id (KID) value of a catchment covering more than one
    pixel, fetch the np array of the pre-calculated weight from the hdf5 
    file. N.b. this file has been pre-calculated, using the function
    Catchment_Weights_Meta(), and is a 65mb file required for the code to
    run correctly.
    Requies h5py module.
    Assumes relative path.
    '''
    with h5py.File('Catchments/Weights/catchment_weights.h5','r') as fp:
        tmparray = fp[tmp_ebid]
        tmparray = tmparray[:,:]
    return tmparray

def Height_CORDEX():
    '''Obtain the orographic height of the CORDEX data from a
    metadata file.
    '''
    NCfpth = 'Metadata/orog_EUR-11_ECMWF-ERAINT_evaluation_r1i1p1_KNMI-RACMO22E_v1_fx.nc'
    #print 'Reading file:',NCfpth
    tmp = Dataset(NCfpth,'r')          # Use NetCDF4 to read the file
    tmp.close                          # Close connection to the file
    orog = tmp.variables['orog']
    return orog


def MT_Means_Over_Lake(nc_path, lake_file, lake_num, outputprefix,
                       tt=None,plots = False,rprt_tme=False):
    '''Purpose             
    This program is the main wrapper to execute the ECCO 
    project functions. It is designed to be executed within the MultiProcessing
    module. Where one Lake is processed by one Process. All other functions in
    this file are called within this function. It essnetially loads data, gets
    a weight mask, calculates the orographic offset, error correction value for 
    the header, weights the climate data, and creates a folder/file system with
    the appropriate name (generated from metadata), and the writes out a copre-
    ssed text file.
    Inputs
    Nc_path is the file path for the CORDEX NetCDF file
    lake_data is file path and filename to GeoJSON file
    lake_num is the lake number to be processed (the int feature number from 
    within the lake_data file) outputprefix is the directory for the outputs
    (it will be then nested according to netcdf file name)
    Outputs 
    Compressed text files with a time-series of values weighted averages 
    over each lake.
    Notes 
    Multiprocessing is not shared memory, so need to load data for each process.
    A much more heavily commented version of this code exists (Speed_MeUp.ipnyb)
    which was where this function was developed and tested.

    This code is in version 2 as it has now been alterd to use the shapefile
    lake data (rather than the GeoJSON data) which contains 127,000 lakes.
    '''
    if rprt_tme:
        a = clock.time()
    num = lake_num
    ShapeData = osgeo.ogr.Open(lake_file)
    TheLayer = ShapeData.GetLayer(iLayer=0)
    feature1 = TheLayer.GetFeature(num) 
    lake_feature = feature1.ExportToJson(as_object=True)
    lake_cart = Path_LkIsl_ShpFile(lake_feature['geometry']['coordinates']) 
    EB_id = lake_feature['properties']['EBhex']
    clim_dat,rlat,rlon,time,metadata,txtfname = Read_CORDEX_V2(nc_path)
    din = clim_dat[:,:,:]
    vname, m1, m2, dexp, m3, m4, m5, m6, drange_orignial = metadata 
    if rprt_tme:
        b = clock.time()    
    if plots:
        Preview_Lake(lake_cart)        
        print 'Island aware Area(km^2)=', Area_Lake_and_Islands(lake_cart),         
        print ', No. xy lake boundary points=',len(lake_cart.vertices)

    lake_rprj = Path_Reproj(lake_cart,False)             
    sub_clim,sub_rlat,sub_rlon = TrimToLake(lake_rprj,clim_dat[0,:,:],rlat,
                                            rlon,off = 3, show = False) 
    tcheck1= clock.time()
    if num == 266232:
        weight_mask = np.load('Lakes/Weights/4a8a86.npy')
    else:
        weight_mask = Pixel_Weights(lake_rprj,sub_clim,sub_rlat,sub_rlon)
    print 'Weighted mask:',clock.time() - tcheck1
    # If the data can easily be adjusted by height then send it for this calculation
    type_of_data = clim_dat.standard_name
    #print 'working on:',type_of_data
    if ((type_of_data == 'air_temperature')| (type_of_data == 'surface_air_pressure')):
        orog = Height_CORDEX()                               # Access the height data 
        lake_altitude=lake_feature['properties']['vfp_mean']
        sub_orog,sub_rlat,sub_rlon = TrimToLake(lake_rprj,orog,rlat,
                                            rlon,off = 3, show = False)
        hght,offset = Orographic_Adjustment(weight_mask,sub_orog,
                                        lake_altitude,clim_dat,chatty=False)
    else:
        hght = -999.   # If no offset is to calculated then
        offset = -999. # just set them to missing data.
    if plots:
        Show_LakeAndData(lake_rprj,clim_dat[0,:,:],rlat,rlon,zoom=8.)
        Preview_Weights(lake_rprj,weight_mask,sub_rlat,sub_rlon)
    if rprt_tme:
        c = clock.time()
    # Here is where the weighting occurs
    pix_truth = (weight_mask > 0.0)      # Count how many times the weight mask is
    pxnum = len(weight_mask[pix_truth])  #  above 0.0 (i.e. how many pixels of data are needed)

    if pxnum == 1: 
        sub_clim,sub_rlat,sub_rlon = TrimToLake3D(lake_rprj,clim_dat,rlat,rlon,
                                                    off = 3, show = False)
        keypix = weight_mask == 1.
        tlist = sub_clim[:,keypix]
    if pxnum > 1:  
        # Now, If there is more than one pixel, weighting needs to be applied...
        tcheck2 = clock.time()
        sub_clim,sub_rlat,sub_rlon = TrimToLake3D(lake_rprj,din,rlat,rlon,
                                                    off = 3, show = False)
        print 'Trim clim dat to Lake:',clock.time() - tcheck2
        tcheck3 = clock.time()
        tlist = Weighted_Mean_3D(weight_mask, sub_clim, chatty=False)
        print 'Calc weighted mean:',clock.time() - tcheck3
    if rprt_tme:
        d = clock.time()
    idnew = EB_id[2:]
    fnm_head = vname+'_'
    hcreate = 'Height offset = %f  Data = %s, Time range = %s  Scenario = %s  Lake = %s'%(\
                offset, clim_dat.long_name, drange_orignial, dexp, EB_id[2:])
    Folder_Create(outputprefix,fnm_head,EB_id)
    pathname = os.path.join(outputprefix, '_'.join([m1, m2, m3, m4, m5, m6]),
                            idnew[:2], idnew[:4], idnew)
    if not os.path.exists(pathname): os.makedirs(pathname)
    np.savetxt(os.path.join(pathname, txtfname+'txt.gz'),
               tlist,fmt='%7.3f',newline='\n', header=hcreate)
    if rprt_tme:
        e = clock.time()
        #print '\n'
        print num,' ',EB_id[2:],' Area(km^2)=', Area_Lake_and_Islands(lake_cart),' no pix:',pxnum
        print ('%4.2f sec : Read Data:'%(b-a))
        print ('%4.2f sec : Calculated height offset and Weighting Mask'%(c-b))
        print ('%4.2f sec : Weighted time-series'%(d-c))
        print ('%4.2f sec : Folder path and file creation'%(e-d))
        print ('%4.2f sec : Total'%(e-a))
    return



def MT_Gen_SWeights(nc_path, lake_file, lake_num, outputprefix, threeD=True,
                       tt=None,rprt_tme=False):
    '''
    Purpose:          
    This program Generates metadata files used to speed up the final runs.
    This was forthe lakes (not catchments), and created a metadata text
    file, to be used as a mini-database in pandas.
    '''
    estart = clock.time()
    if rprt_tme == True:
        a = clock.time()
    num = lake_num
    orog = Height_CORDEX()
    ShapeData = osgeo.ogr.Open(lake_file)
    TheLayer = ShapeData.GetLayer(iLayer=0)
    feature1 = TheLayer.GetFeature(num) 
    lake_feature = feature1.ExportToJson(as_object=True)
    lake_cart = Path_LkIsl_ShpFile(lake_feature['geometry']['coordinates']) 
    EB_id = lake_feature['properties']['EBhex']
    

    lake_altitude=lake_feature['properties']['vfp_mean']
    clim_dat,rlat,rlon,time,metadata,txtfname = Read_CORDEX_V2(nc_path)
    vname, m1, m2, dexp, m3, m4, m5, m6, drange_orignial = metadata 
    if rprt_tme == True:
        b = clock.time()    

    lake_rprj = Path_Reproj(lake_cart,False)

    sub_clim,sub_rlat,sub_rlon = TrimToLake(lake_rprj,clim_dat[0,:,:],rlat,
                                            rlon,off = 3, show = False) 
    weight_mask = Pixel_Weights(lake_rprj,sub_clim,sub_rlat,sub_rlon)
    if rprt_tme == True:
        c = clock.time()   

    pix_truth = (weight_mask > 0.0)    # Count how many times the weight mask is
    pxnum = len(weight_mask[pix_truth])  #  above 0.0 (i.e. how many pixels of data are needed)
    
    etime = clock.time() - estart # (Calculate elapsed time in sec)
    # Calculate the CORDEX pixels (useful only for 1-pixel lakes)
    ypix = -99
    xpix = -99
    if pxnum == 1:
        xxx,yyy = Get_LatLonLim(lake_rprj.vertices)  # Find upp./low.lake lims.
        ypix = (Closest(rlat,yyy[0]))                # For lakes of one pixel  
        xpix = (Closest(rlon,xxx[0]))
    if pxnum < 1:
        pxnum = 1  # Small bug where it thinks lakes dont exist, no biggy...

    return num,EB_id[2:], Area_Lake_and_Islands(lake_cart),pxnum,etime,ypix,xpix


def nc_write_func(fileout,EB_id,standard_name,long_name,units,area,tlist):
    ''' Write a variable to an open NC file. This is done inside
    a function because, when the function is left, all entries in
    the namespace dissapear. Hence, the non-bug in NetCDF4 of
    bleeding memory in loops is no issue.
    This was just for testing. It turns out the HDF5 is way more efficent here.
    '''
    lake = fileout.createGroup(EB_id)                       
    values = lake.createVariable('values','f4',('time',))   #,zlib=True,least_significant_digit=3)
    values.standard_name=standard_name
    values.long_name=long_name
    values.units=units
    values.lake_area=str(area)+'km^2'
    values.height_adjustment=str(-99.)
    values[:] = tlist
    return

def OnePix_HOffset(lake_altitude,orog,var_type):
    ''' Purpose    -   To replace the height offset calculation for cases where only one pixel
                       is needed (and has already been identified in preprocessed metadata).
                       This speeds up program execution.I.e. Replaces Orographic_Adjustment()
                       in most cases.
    '''
    if (var_type == 'air_temperature'): 
            ELR = 6.49/1000.                        # Environmental lapse rate of 6.49K/1,000m
            offset = (orog - lake_altitude) * ELR
    elif(var_type == 'surface_air_pressure'):
            Pchng = 1200./100.                      # Pressue change of 1.2kPa/100m
            offset = (orog - lake_altitude) * Pchng
    else:
        offset = -999.                              # Missing val. if not temp. or press. data
    return offset


def Orographic_Adjustment(weight_mask,orog,lake_altitude,clim_dat,chatty):
    '''Purpose    - Reads in the 2D weight mask (from Pixel_Weights 
    function) and the trimmed data from TrimToLake and returns a
    weighted mean. This should be iterated for each time-step.
    Input    - weight_mask: Pixel weights
             - sub_clim: subset of the climate data (one time slice,
               lat lon subset)
             - chatty: a true or false statement, if true, it will
               print some info on the weighting process.
    Output   - val_out: the weighted mean value of height diffrence
             - offset : the calculated paramater offset (e.g. surface
               air temp) calculated for each pixel
                        individually, weighted into a single value.
    '''
    type_of_data = clim_dat.standard_name
    ELR = 6.49/1000    # Environmental lapse rate of 6.49K per 1,000m
    Pchng = 1200/100    # Pressue change of 1.2kPa per 100m
    offset_made = False   # Initialise a boolean test as False
    aaa = np.where(weight_mask > 0.000)   # Index where weights exist
    if (len(aaa[0]) == 0):
        print 'Error: no lake cover identified! :('
    if (chatty == True):
        print 'Lake covers',len(weight_mask[aaa]),' pixels'
        print 'Actual weight values are :',weight_mask[aaa]
        print 'GeoJSON height of lake (m):',lake_altitude
        print 'Height of individual pixels in data (m) is:',orog[aaa]
        print 'Cum. sum of pixel weights (should end as 1.0):',cumsum(weight_mask[aaa])
    val_out = 0
    tmp_diffs =[]
    offset = 0
    for n in xrange(len(aaa[0])):
        val_out = val_out + weight_mask[aaa[0][n],aaa[1][n]] * orog[aaa[0][n],aaa[1][n]]
        
    for n in xrange(len(aaa[0])):  # For every lake pixel, calc height adjustment
        if type_of_data == 'air_temperature':
            #print 'found air temp data'
            #print 'blah2',(orog[aaa[0][n],aaa[1][n]]- lake_altitude) * ELR
            tmp_diffs.append((orog[aaa[0][n],aaa[1][n]]- lake_altitude) * ELR)
            offset_made =True

    for n in xrange(len(aaa[0])):  # For every lake pixel, calc height adjustment
        if type_of_data == 'surface_air_pressure':
            tmp_diffs.append((orog[aaa[0][n],aaa[1][n]]- lake_altitude) * Pchng)
            offset_made =True

    if offset_made == True:
        for i,n in enumerate(xrange(len(aaa[0]))):
            offset = offset + weight_mask[aaa[0][n],aaa[1][n]] * tmp_diffs[i]
            #print i, tmp_diffs[i],weight_mask[aaa[0][n],aaa[1][n]],offset
    else:
        offset = 0.0   # Set the offset to nothing if no offset has been calculated
        
    return val_out, offset


def Path_Make(coord):
    '''Purpose  - Create a Polygon from the Lake vectors, as a Matplotlib.Path.Path object, and nothing else.\n",
    Input    - Lake coordinate data as an x,y list (from the GeoJSON file).\n",
    Output   - Path_out (A Matplotlib.Path object) \n",
    Requires - Matplotlib.Path import Path\n",
    Notes   - The way the lake addressing goes is Lake[0][0][0] where the first 0 is the individual lake.\n",
    Lake[0][0][:] will give you every x,y element. Lake[0][0][0][0] will give the first x-element of the first lake.\n",
    Lake[0][0][0][1] will give the first y-element etc.
    nb. This is version 2 created Nov 10th 2014
    '''
    verts =[] ; codes=[]
    pth_ln = int(len(coord) - 1)
    for i,n in enumerate(coord[0:pth_ln]):
        verts.append(n)
        if i == 0:
            codes.append(Path.MOVETO)
        if ((i >= 1)&(i < (pth_ln - 1))):
            codes.append(Path.LINETO)
        if i == (pth_ln - 1):
            codes.append(Path.CLOSEPOLY)
    path_out = Path(verts, codes)
    return path_out


def Path_LkIsl_ShpFile(lake_path):
    '''Purpose  - This function replaces a direct call to Path_Make() with a function
    able to add islands. This is reading data from the shapefile, not GeoJSON.
    Input - Feature['geometry']['coordinates'] (taken from a shp file layer)
    Output- Path (a Matplotlib.Path object)
    nb. edited Nov 10th 2014
    This is the function which should be used to create island-smart lakes from shp file.
    '''
    num_rings = len(lake_path)
    pathouts=[Path_Make(lake_path[ringi]) for ringi in range(num_rings)]
    psep =[]
    csep =[]
    for n in xrange(num_rings):
        ptmp = pathouts[n]
        psep.append(ptmp.vertices)
        csep.append(ptmp.codes)
    lk_stack = np.concatenate(psep, axis=0)
    lk_codes = np.concatenate(csep, axis=0)
    path_wisl = Path(lk_stack, lk_codes)
    return path_wisl


def PMake(coord):
    '''Purpose  - Create a Polygon from the Lake vectors, as a Matplotlib.Path.Path object, and nothing else.
    Input    - Lake coordinate data as an x,y list (from the GeoJSON file).
    Output   - x,y vectors of path, and codes (used to create a Path Matplotlib object)
    Requires - Matplotlib.Path import Path
    Notes   - The way the lake addressing goes is Lake[0][0][0] where the first 0 is the individual lake.
    Lake[0][0][:] will give you every x,y element. Lake[0][0][0][0] will give the first x-element of the first lake.
    Lake[0][0][0][1] will give the first y-element etc.
    This works with the GeoJSON format file.
    '''
    verts =[] ; codes=[]
    pth_ln = int(len(coord) - 1) 
    for i,n in enumerate(coord[0:pth_ln]):
        #print i,n
        verts.append(n)
        if i == 0:
            codes.append(Path.MOVETO)
        if ((i >= 1)&(i < (pth_ln - 1))):
            codes.append(Path.LINETO)  
        if i == (pth_ln - 1):
            codes.append(Path.CLOSEPOLY)
    path_out = Path(verts, codes)  
    return verts,codes


def Path_Lake_and_Islands_old(num,lake_path):
    '''Purpose  - This function replaces Path_Make() with a function able to add islands.
    Input - Number of Lake
          - List of Lake coordinates read in from the datafile
    Output- Path (a Matplotlib.Path object)
    '''
    num_rings = len(lake_path[num][0]) 
    pathouts = [PMake(lake_path[num][0][ringi]) for ringi in range(num_rings)]
    coordinates = [e1 for e1, e2 in pathouts]
    codeouts = [e2 for e1, e2 in pathouts]
    lk_stack = np.concatenate(coordinates, axis=0)
    lk_codes = np.concatenate(codeouts, axis=0)
    path_wisl = Path(lk_stack, lk_codes)
    return path_wisl

def Path_Lake_and_Islands(num,lake_path):
    '''Purpose  - This function replaces Path_Make() with a function able to add islands.
    Input - Number of Lake
          - List of Lake coordinates read in from the datafile
    Output- Path (a Matplotlib.Path object)
    '''
    num_rings = len(lake_path[num]) 
    pathouts = [PMake(lake_path[num][ringi]) for ringi in range(num_rings)]
    coordinates = [e1 for e1, e2 in pathouts]
    codeouts = [e2 for e1, e2 in pathouts]
    lk_stack = np.concatenate(coordinates, axis=0)
    lk_codes = np.concatenate(codeouts, axis=0)
    path_wisl = Path(lk_stack, lk_codes)
    return path_wisl

def Path_Reproj(path_in,INV):
    '''This program reprojects a path object from cartesian lat lon 
    to Spherical reprojected coordinates.
    Note, this could easily be changed in the future to use diffrent
    coordinates as requirements of model input changes.
    Currently works on Lat Lon to rotated Lat Lon of EUR-11 model.
    Input    -   Lake Path object
             -   change_proj, a projection object from Get_Proj() to
                 transform coordinates.
             -   INV keyword as a boolean. If false the 
                 transformation goes forwards
                 meaning cartesian lon/lat are converted by the 
                 reprojection. If it
                 is true, however, the tranformation goes in reverse.
                 Changing reprojected data back to cartesian data.
    Output   -   Transformed Path object. Note the codes are intact.        
    '''
    verts =[] ; codes=[] 
    pth_ln = int(len(path_in.vertices) - 1) 
    for i,n in enumerate(path_in.vertices):
        if (INV == False):
            rx,ry = Calc_Coordinates(n[0],n[1])
        if (INV == True):
            print 'Still must add inverse projection functionaility!'
        tmp_rverts = [rx,ry]
        verts.append(tmp_rverts)         # Use reprojected verticies
        codes.append(path_in.codes[i])   # Use path codes already set
    path_out = Path(verts, codes)
    return path_out




def Pixel_Weights(lake_in, datin,lat_atts,lon_atts):
    '''
    Purpose - Provides a 2D array matching the input array, which has
    the weights of pixels to calc. lake mean
    
    Input   - lake_in: A Path object
            - datin: A 2D array of values: either real data, or fake 
              data generated by Gen_Test_Dat()
            - Lat subscripts: list of y values in rotated coordinates
            - Lon subscripts: list of x values in rotated coordinates
            
    Output  - pix_weights: A 2D array matching the input data of the 
              fractional area of lake per pixel
    
    Notes   - For each grid cell of the provided 2D data the 
              individual bounding box of each pixel is calculated
                as a bounding box object (Matplotlib.Transform.Bbox).
              Using the intersects_bbox() function, a logic test 
              determines if any vectors of the Lake's path are within
              the bounding box. If the condition is false (no lake 
              within pixel), the value of the mask is set to 0.0. 
              If the condition is true (lake within pixel), the
              fractional area of the lake polygon within the polygon
              within the pixel is calculated. This requires the area 
              of the lake polygon within the pixel (in km^2) is 
              calculated, and divided by the area of the total lake
              (also in km^2). The areas are calculated by the custom
              functions EqArea() and PolyArea2D().

              Nb. As of 23rd Oct 2014, this was updated to use Area_
              Lake_and_Islands() function to calculate area. This
              requires the BBOX_PathCode_Fix() function also, to 
              correct an error, where after trimming, the 79 (close)
              path codes may be removed.
    '''
    fillone = 0  # A small fix to ensure a bug-fix wont break (fix was for cases of puddles)
    lout = lake_in
    pix_weights = np.zeros(np.shape(datin))
    cnt = 0
    latstep = (float(lat_atts[-1]) - float(lat_atts[0])) / float(len(lat_atts))
    lonstep = (float(lon_atts[-1]) - float(lon_atts[0])) / float(len(lon_atts))
    xvar = lon_atts[0] ; yvar= lat_atts[0]       # Initialize the looping indexes to first positions
    for x in xrange(len(lon_atts)):
        if(x == 0):  
            x_ll = lon_atts[0]            # ..then set the lower left corner to be the first lon val.
            x_ur = lonstep + x_ll         # and the upper right, to be the x_ll + increment
        if(x > 0):                        # If it is already in the loop, use the older x_ll, and x_ur
            x_ll = x_ur                   # to continue calculating the bounding box area. 
            x_ur = x_ll + lonstep    
        for y in xrange(len(lat_atts)):
            if(y == 0):
                y_ll = lat_atts[0]
                y_ur = latstep + y_ll
            if(y > 0):
                y_ll = y_ur
                y_ur = y_ll + latstep  
            lims =np.array([[x_ll,y_ll],[x_ur,y_ur]])  # Construct a 2D np.array for the BB object
            tmpbb = Bbox(lims)
            
            #acheck = 0  # My own logical test to see if pixels overlap with lake
            #for n in lout.vertices:
            #    if ((n[0] > lims[0][0])&(n[0]<lims[1][0])&(n[1] >lims[0][1])&\
            #        (n[1]<lims[1][1])):
            #        acheck = 1

            ben_test = 0
            acheck = [((n[0] > lims[0][0])&(n[0]<lims[1][0])&(n[1] >lims[0][1])&\
                (n[1]<lims[1][1])) for n in lout.vertices]
            acheck = np.array(acheck)
            tst = np.where(acheck == True)  # Use this snippet to test if the coordinates
            if len(tst[0]) > 0:             # Fall within a bounding box range 
                ben_test = 1
            #print 'Values in BB?',ben_test
            #print 'Condition ',acheck
            #test = lout.intersects_bbox(tmpbb)     # This isn't working right with the islands
            if (ben_test == 1):                     # If Lake within the Bbox test will be 1 (no lake, test = 0)
                cnt = cnt + 1                       # Counter is just for testing.
                #pix_weights[y,x] = 0.0             # For pixels where Lake is, calc. and write the frac area (%)
                sub_lout =[] ; area_sub = []
                #print 'Here!',lims
                sub_lout = lout.clip_to_bbox(tmpbb,inside ='True')
                #area_sub = Poly_Area2D(EqArea(sub_lout.vertices))
                #lkarea = Poly_Area2D(EqArea(lout.vertices))
                sub_lout = BBOX_PathCode_Fix(sub_lout) # Corrects the 79 code error!
                area_sub = Area_Lake_and_Islands(sub_lout)
                lkarea = Area_Lake_and_Islands(lout)
                if lkarea > 0.1:  # Added 4 Dec 2014 to combat a bug from small lakes
                    pix_weights[y,x]= float(area_sub)/float(lkarea)
                if lkarea < 0.1:
                    if fillone == 0:
                        #print 'CHECK: INSIDE YES!'
                        pix_weights[y,x] = 1.
                    fillone = fillone + 1   
                    # Fractional area (0-1.0) of lake within a given pixel
    return pix_weights


def Poly_Area2D(poly):
    '''Purpose - This function implements Green's Theorem to
    calculate area of a polygon in a 2D co-ordinate system.
    Input   - A polygon, as the verticies of a Matplotlib.Path 
    Notes   - More info at http://code.activestate.com/recipes/
                           578275-2d-polygon-area
    Join this function with EqArea function to find the area of a
    lake in km^2.
    Example - area = Poly_area2D(EqArea(A_Lake_Path.vertices))
    '''
    total = 0.0
    N = len(poly)
    for i in range(N):
        v1 = poly[i]
        v2 = poly[(i+1) % N]
        total += v1[0]*v2[1] - v1[1]*v2[0]
    return '%6.2f'%abs(total/2.)  


def Plot_LakeAndData_Save(lake_in,cdat,rlat,rlon,zoom):
    ''' Like Show_LakeAndData(), but this one outputs the plot object
    without displaying it to the screen. Creates a nice plot for the
    user to check the lake is in the right place on the data. Specify
    a zoom to change perspective as required.
    Input   - lake_in : a Path object containing a lake
            - cdat: a 2D array of lat/lon holding climate variables
            - rlon and rlat : rotated lon/lat vectors of cdat
            - zoom : a positive floating point value speciying the
              area in decimal degrees to plot around the lake.You can
              use this to zoom on the lake!
    Output  - A plot object of the lake overlaid on the data, ready
              to save to file.
    '''
    if ((zoom < 0) | (zoom != zoom)):
        zoom = 0.
    xmaxmin,ymaxmin = Get_LatLonLim(lake_in.vertices)

    fig2 = plt.figure()
    ax1 = fig2.add_subplot(111)
    patch = patches.PathPatch(lake_in, facecolor='#06ebf6', lw=1)
    ax1.add_patch(patch)      # ADD LAKE
    ax1.set_xlim(xmaxmin[1]-zoom,xmaxmin[0]+zoom)
    ax1.set_ylim(ymaxmin[1]-zoom,ymaxmin[0]+zoom)
    ax1.set_ylabel('Lat. (Deg. N)')
    ax1.set_xlabel('Lon. (Deg. E)')
    ax1.set_title('Lake and Climate data Overlaid'+'\n'+'Zoom of '+str(zoom) +' degrees around lake')
    ax1.imshow(cdat,interpolation='none', cmap=cm.RdBu,extent=[rlon[0],rlon[-1],rlat[0],rlat[-1]],origin='lower')
    #plt.show(fig2)
    return fig2

def Preview_Lake(lake_in):
    '''Purpose   - Plot to screen a specified lake from path object.
    Only meant for preview purposes.
    Input     - Matplotlib.Path object
    '''
    xtmp=[] ; ytmp=[]
    for i in xrange(len(lake_in.vertices)):
        xtmp.append(lake_in.vertices[i][0])
        ytmp.append(lake_in.vertices[i][1]) 
    fig2 = plt.figure()
    ax1 = fig2.add_subplot(111)
    patch = patches.PathPatch(lake_in, facecolor='#06ebf6', lw=1)
    ax1.set_axis_bgcolor('#87c540')   # Hexcolor is a green
    ax1.set_xlim(min(xtmp)-0.05,max(xtmp)+0.05)
    ax1.set_ylim(min(ytmp)-0.05,max(ytmp)+0.05)
    ax1.add_patch(patch)              # Add the lake object here
    ax1.set_ylabel('Lat. (Deg. N)')
    ax1.set_xlabel('Lon. (Deg. E)')
    ax1.set_title('Lake preview')
    plt.show()
    return


def Preview_Weights(lake_in,pix_weights,lat_atts,lon_atts):
    '''
    Purpose  - To plot a preview of the pixel weight mask with the lake overlaid.
                 Weights calculated from the custom Pixel_Weights() function.
    Input    - Lake Path
             - 2D array of pixel weights
             - subscripts of lon and lats from SubsetClimDat() function.
    '''
    lout = lake_in
    fig2 = plt.figure()
    ax1 = fig2.add_subplot(111)
    patch = patches.PathPatch(lout, facecolor='#06ebf6', lw=1)
    ax1.add_patch(patch)     # ADD LAKE
    ax1.set_xlim(lon_atts[0],lon_atts[-1])
    ax1.set_ylim(lat_atts[0],lat_atts[-1])
    ax1.set_ylabel('Lat. (Deg. N)')
    ax1.set_xlabel('Lon. (Deg. E)') 
    ax1.imshow(pix_weights,interpolation='none', cmap=cm.Greys, extent=[lon_atts[0],lon_atts[-1],lat_atts[0],
                lat_atts[-1]],origin='lower')
    plt.show()
    return

def Read_CORDEX_V2(nc_path):
    ''' PURPOSE
    After NetCDF4 data has been downloaded from a WGet script (obtained
    from ESGF website), this subroutine takes the path/fname as an input,
    and returns the data and metadata. In this case, metadata re the details
    of the file are taken from the filename rather than the netcdf contents.

    INPUT
    The path to the netcdf file only

    OUTPUT
    clim_dat : the actual climate data from the file
    rlat / rlon : rotated lat/lon coordinate dimensions
    time : time dimension
    metadata : a list of metadata extracted from the filename
    drange : the daterange of the file (extracted from filename)
    txtfname : the filename of the output text file (made of variable name
                and date range).
    '''
    nc_fname = os.path.basename(nc_path)
    metadata = os.path.splitext(nc_fname)[0].split('_')
    vname, m1, m2, dexp, m3, m4, m5, m6, drange_orignial = metadata
    drange = '_'.join([dexp, drange_orignial])
    txtfname = '_'.join([vname, drange])
    nc = Dataset(nc_path)                 # Read the NetCDF Data 
    clim_dat = nc.variables[vname]
    rlat = nc.variables['rlat']
    rlon = nc.variables['rlon']
    time = nc.variables['time']
    return clim_dat,rlat,rlon,time,metadata,txtfname

def Read_Lakes(file_in):
    '''Purpose - Use Json module to read GeoJSON Lake data
    Input   - File name including path (string)
    Output  - Various arrays containing Lake data
    '''
    with open(file_in) as f:
        data = json.load(f)
    lake_id =[] ; lake_path=[] ; lake_geometry=[] ; lake_name=[]
    lake_altitude=[] ; feno_lake_id=[]
    for feature in data['features']:
        lake_geometry.append(feature['geometry']['type'])
        lake_id.append(feature['id'])
        lake_name.append(feature['properties']['comsat.Lake'])
        lake_path.append(feature['geometry']['coordinates'])
        lake_altitude.append(feature['properties']['comsat.Altitude'])
        feno_lake_id.append(feature['properties']['fennoscan.lake_id'])
    #print 'Read data for', len(lake_id), ' lakes.'
    return lake_geometry, lake_id, lake_name,lake_path,lake_altitude,feno_lake_id


def Read_LakesV2(file_in):
    '''Purpose - Use Json module to read GeoJSON Lake data
    Input   - File name including path (string)
    Output  - Various arrays containing Lake data
    '''
    with open(file_in) as f:
        data = json.load(f)
    EB_id = [feature['properties']['EBhex'] for feature in data['features']]
    lake_path = [feature['geometry']['coordinates'] for feature in data['features']]
    lake_altitude = [feature['properties']['vfp_mean'] for feature in data['features']]
    return EB_id, lake_path, lake_altitude 

def Read_LakesV3(file_in):
    '''Purpose - Use Json module to read GeoJSON Lake data (ran into trouble with v2)
    Input   - File name including path (string)
    Output  - Various arrays containing Lake data
    '''
    with open(file_in) as f:
        data = simplejson.loads(f)   # diffrence here on the loads
    EB_id = [feature['properties']['EBhex'] for feature in data['features']]
    lake_path = [feature['geometry']['coordinates'] for feature in data['features']]
    lake_altitude = [feature['properties']['vfp_mean'] for feature in data['features']]
    return EB_id, lake_path, lake_altitude 

def Reproj_Catchment(lake_feature,chatty=False):
    '''
    Takes the lowest level data from the Catchment shape file feature, and
    reprojects the coordinate system using PROJ (PyProj) to a WGS85 system.
    '''
    import pyproj

    p1 = pyproj.Proj("+init=EPSG:25833")   # Projected system for lake catchments
    p2 = pyproj.Proj("+init=EPSG:4326")    # World Geodetic system
    tmp1=[]
    
    if lake_feature['geometry']['type'] == 'MultiPolygon':
        if chatty == True: print 'Multi-polygon type'
        for n in lake_feature['geometry']['coordinates'][0][0]:  #For every vertice pair...
            x,y=(pyproj.transform(p1,p2,n[0],n[1],z=None,radians=False)) # trans. coords. to WGS85
            tmp1.append([x,y])
        return tmp1
    else:
        if chatty == True: print 'Simple Polygon type'
        for n in lake_feature['geometry']['coordinates'][0]:  #For every vertice pair...
            x,y=(pyproj.transform(p1,p2,n[0],n[1],z=None,radians=False)) # trans. coords. to WGS85
            tmp1.append([x,y])
        return tmp1

def Show_LakeAndData(lake_in,cdat,rlat,rlon,zoom):
    ''' This creates a nice plot for the user to check the lake is in
    the right place on the data. Specify a zoom to change perspective
    as required.
    Input   - lake_in : a Path object containing a lake
            - cdat: a 2D array of lat/lon holding climate variables
            - rlon and rlat : rotated lon/lat vectors of cdat
            - zoom : a positive floating point value speciying the 
            area in decimal degrees to plot around the lake. You can
            use this to zoom on the lake!
    Output  - A nice plot of the lake overlaid on the data
    '''
    if ((zoom < 0) | (zoom != zoom)):
        zoom = 0.
    xmaxmin,ymaxmin = Get_LatLonLim(lake_in.vertices)

    fig2 = plt.figure()
    ax1 = fig2.add_subplot(111)
    patch = patches.PathPatch(lake_in, facecolor='#06ebf6', lw=1)
    ax1.add_patch(patch)      # ADD LAKE
    ax1.set_xlim(xmaxmin[1]-zoom,xmaxmin[0]+zoom)
    ax1.set_ylim(ymaxmin[1]-zoom,ymaxmin[0]+zoom)
    ax1.set_ylabel('Lat. (Deg. N)')
    ax1.set_xlabel('Lon. (Deg. E)')
    ax1.set_title('Lake and Climate data Overlaid'+'\n'+'Zoom of '+str(zoom) +' degrees around lake')
    ax1.imshow(cdat,interpolation='none', cmap=cm.RdBu,extent=[rlon[0],rlon[-1],rlat[0],rlat[-1]],origin='lower')
    plt.show(fig2)
    return

def TrimToLake(lake_in,Cdat,rlat,rlon,off,show):
    ''' Purpose   - To go from the full CORDEX lat lon array (one
        time slice), to a small np.array subset over the immediate
        pixels surrounding the lake (with a few eitherside). This
        will speed up the execution of the pixel weight calculation
        code.
        
    Input     -  lake_in : A Path object, holding the lake data
              -  Cdat : Climate Data, as a 2D array (lat,lon) sliced
                 from the cordex array
              -  rlat : Rotated Latitude attributes of the Cdat array
              -  rlon : Rotated Lontiude attributes of the Cdat array
              -  off  : an offset value (in pixels) to expand the 
                 area around the lake. By default, if the offset
                 enterd is less than 3, it will be set to 3 pixels
                 (as from testing this prevents errors).
              -  show : A keyword set to either True or False 
                 depending on if you want to see the result plotted
                 to the screen or not.
                        
    Output    -  data_sub : Np.Array subset of Cdat. This will be
                 used to create a weighted mask.
              -  rlat_subs: subset of the rotated latitude
              -  rlon_subs: subset of the rotated lontiude
              
    Example   - grid_subset,subset_rlats,subset_rlons = 
                Subset_ClimDat(lake_in,climdata[0,:,:],rlat,rlon,
                show=True)
    
    Notes     - The zoom2 factor, is a bit arbitrary, is is designed
                to include several pixels comfortably either side
                of the lake. While this may be a waste, and perhaps
                could be cut, I perfer to leave it in, as the lakes
                have some unusual shapes, and leaving pixels to work
                with seems like a good idea for now.
    '''
    if ((off < 3) | (off != off)):
        off = 3
    xxx,yyy = Get_LatLonLim(lake_in.vertices)  # Get bounds of a lake
    ymx = (Closest(rlat,yyy[0])) + off         # Gather the max and minimum range
    ymn = (Closest(rlat,yyy[1])) - off         # also add an offset (measured in pixels)
    xmx = (Closest(rlon,xxx[0])) + off         # Gather the max and minimum range
    xmn = (Closest(rlon,xxx[1])) - off         # also add an offset (measured in pixels)
    sub_rlat = rlat[ymn:ymx]
    sub_rlon = rlon[xmn:xmx]
    data_sub = Cdat[ymn:ymx,xmn:xmx]
    if show == True:          # (If show is set to True, then make a plot to show what's what)
        fig3 = plt.figure()
        ax1 = fig3.add_subplot(111)
        patch = patches.PathPatch(lake_in, facecolor='#06ebf6', lw=1)
        ax1.add_patch(patch)      # ADD LAKE
        ax1.set_ylabel('Lat. (Deg. N)')
        ax1.set_xlabel('Lon. (Deg. E)')
        ax1.set_title('Preview of trimmed climate data with lake overlaid')
        ax1.imshow(data_sub,interpolation='none', cmap=cm.RdBu,
                   extent=[sub_rlon[0],sub_rlon[-1],sub_rlat[0],sub_rlat[-1]],origin='lower')
        plt.show(fig3)    
        print np.shape(data_sub),type(data_sub)
    return data_sub,sub_rlat,sub_rlon

def Tmp_CORDEX_Read():
    '''Temporay way to read the CORDEX data, while I am testing the 
    software and getting the remote access working still.
    '''
    NCfpth = '/uio/kant/geo-metos-u1/blaken/Downloads/tas_EUR-11_ICHEC-EC-EARTH_rcp85_r3i1p1_DMI-HIRHAM5_v1_day_20960101-21001231.nc'
    #print 'Reading file:',NCfpth
    tmp = Dataset(NCfpth,'r')              # Use NetCDF4 to read the file
    tmp.close                              # Close the connection to the file 
    tas = tmp.variables['tas']
    rlat = tmp.variables['rlat']
    rlon = tmp.variables['rlon']
    time = tmp.variables['time']
    drange = NCfpth  # This is a unicode date range of the data stripped from the filename, pass it to the text
    drange = drange[-20:-3]   # files as part of    their naming convention
    dexp = tmp.driving_experiment  # variable to hold experiment info
    return tas,rlat,rlon,time,drange,dexp


def TrimToLake3D(lake_in,Cdat,rlat,rlon,off,show):
    ''' Purpose   - To go from the full CORDEX lat lon array (one
        time slice), to a small np.array subset over the immediate
        pixels surrounding the lake (with a few eitherside). This
        will speed up the execution of the code in theory...
        
    Input     -  lake_in : A Path object, holding the lake data
              -  Cdat : Climate Data, as a 3D array (time, lat, lon) sliced
                 from the cordex array
              -  rlat : Rotated Latitude attributes of the Cdat array
              -  rlon : Rotated Lontiude attributes of the Cdat array
              -  off  : an offset value (in pixels) to expand the 
                 area around the lake. By default, if the offset
                 enterd is less than 3, it will be set to 3 pixels
                 (as from testing this prevents errors).
              -  show : A keyword set to either True or False 
                 depending on if you want to see the result plotted
                 to the screen or not.
                        
    Output    -  data_sub : Np.Array subset of Cdat. This will be
                 used to create a weighted mask.
              -  rlat_subs: subset of the rotated latitude
              -  rlon_subs: subset of the rotated lontiude
              
    Example   - grid_subset,subset_rlats,subset_rlons = 
                Subset_ClimDat(lake_in,climdata[:,:,:],rlat,rlon,
                show=True)
    
    Notes     - The zoom2 factor, is a bit arbitrary, is is designed
                to include several pixels comfortably either side
                of the lake. While this may be a waste, and perhaps
                could be cut, I perfer to leave it in, as the lakes
                have some unusual shapes, and leaving pixels to work
                with seems like a good idea for now.
    '''
    if ((off < 3) | (off != off)):
        off = 3
    xxx,yyy = Get_LatLonLim(lake_in.vertices)  # Get bounds of a lake
    ymx = (Closest(rlat,yyy[0])) + off         # Gather the max and minimum range
    ymn = (Closest(rlat,yyy[1])) - off         # also add an offset (measured in pixels)
    xmx = (Closest(rlon,xxx[0])) + off         # Gather the max and minimum range
    xmn = (Closest(rlon,xxx[1])) - off         # also add an offset (measured in pixels)
    sub_rlat = rlat[ymn:ymx]
    sub_rlon = rlon[xmn:xmx]
    data_sub = Cdat[:, ymn:ymx, xmn:xmx]
    return data_sub,sub_rlat,sub_rlon



def TrimToLake3D_plus(lake_in,Cdat,rlat,rlon,off,show):
    ''' Purpose   - To go from the full CORDEX lat lon array (one
        time slice), to a small np.array subset over the immediate
        pixels surrounding the lake (with a few eitherside). This
        will speed up the execution of the code in theory.
        This is a sped-up version of the TrimToLake3D() function.
        As it appears that that was where some big delay was occuring.
        
    Input     -  lake_in : A Path object, holding the lake data
              -  Cdat : Climate Data, as a 3D array (time, lat, lon) sliced
                 from the cordex array
              -  rlat : Rotated Latitude attributes of the Cdat array
              -  rlon : Rotated Lontiude attributes of the Cdat array
              -  off  : an offset value (in pixels) to expand the 
                 area around the lake. By default, if the offset
                 enterd is less than 3, it will be set to 3 pixels
                 (as from testing this prevents errors).
              -  show : A keyword set to either True or False 
                 depending on if you want to see the result plotted
                 to the screen or not.
                        
    Output    -  data_sub : Np.Array subset of Cdat. This will be
                 used to create a weighted mask.
              -  rlat_subs: subset of the rotated latitude
              -  rlon_subs: subset of the rotated lontiude
              
    Example   - grid_subset,subset_rlats,subset_rlons = 
                Subset_ClimDat(lake_in,climdata[:,:,:],rlat,rlon,
                show=True)
    
    Notes     - The zoom2 factor, is a bit arbitrary, is is designed
                to include several pixels comfortably either side
                of the lake. While this may be a waste, and perhaps
                could be cut, I perfer to leave it in, as the lakes
                have some unusual shapes, and leaving pixels to work
                with seems like a good idea for now.
    '''
    if ((off < 3) | (off != off)):
        off = 3
    xxx,yyy = Get_LatLonLim(lake_in.vertices)  # Get bounds of a lake
    ymx = (Closest(rlat,yyy[0])) + off         # Gather the max and minimum range
    ymn = (Closest(rlat,yyy[1])) - off         # also add an offset (measured in pixels)
    xmx = (Closest(rlon,xxx[0])) + off         # Gather the max and minimum range
    xmn = (Closest(rlon,xxx[1])) - off         # also add an offset (measured in pixels)
    sub_rlat = rlat[ymn:ymx]
    sub_rlon = rlon[xmn:xmx]
    data_sub = Cdat[:, ymn:ymx, xmn:xmx]
    return data_sub,sub_rlat,sub_rlon





def Update_Progress(progress):
    '''A nice solution to progress bars, all contained here (no need
    to load packages). Update_Progress() : Displays or updates a
    console progress bar Accepts a float between 0 and 1. Any int
    will be converted to a float. A value under 0 represents a 'halt'
    A value at 1 or bigger represents 100%
    Took this code from http://stackoverflow.com/questions/3160699/
                        python-progress-bar
    Testing shows the p-bar doesn't slow programs (tests showed +0.2
    sec in loops of n=77)
    '''
    barLength = 40 # Modify this to change the length of the p-bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "\r\n"
    block = int(round(barLength*progress))
    text = "\rProgress: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), int(progress*100), status)
    sys.stdout.write(text)
    sys.stdout.flush()
    return

def Weighted_Mean(weight_mask,sub_clim,chatty):
    '''Purpose    - Reads in the 2D weight mask (from Pixel_Weights
    function) and the trimmed data from TrimToLake and returns a
    weighted mean. This should be iterated for each time-step.
    Input    - weight_mask: Pixel weights
             - sub_clim: subset of the climate data (one time slice,
               lat lon subset)
             - chatty: a true or false statement, if true, it will 
               print some info on the weighting process.
    Output   - val_out: the weighted mean value
    '''
    aaa = np.where(weight_mask > 0.000)   # Indx where weights exist
    if (len(aaa[0]) == 0):
        print 'Error: no lake cover identified! :('
    if (chatty == True):
        print 'Lake covers',len(weight_mask[aaa]),' pixels'
        print 'Actual weight values are :',weight_mask[aaa]
        print 'Cum. sum of pixel weights (should end as 1.0):',np.cumsum(weight_mask[aaa])
    val_out = 0
    for n in xrange(len(aaa[0])):
        val_out = val_out + weight_mask[aaa[0][n],aaa[1][n]] * sub_clim[aaa[0][n],aaa[1][n]]
    return val_out


def Weight_Speedup(nc_path, lake_file, lake_num, outputprefix,
                       tt=None,plots = False,rprt_tme=False):
    '''
    Purpose             
    This program generates the weighted mean files for specified lakes
    (the slow ones), to be read back in when processing later for a 
    big speed boost! This only needs to be run one time. 
    '''
    a = clock.time()
    num = lake_num
    ShapeData = osgeo.ogr.Open(lake_file)
    TheLayer = ShapeData.GetLayer(iLayer=0)
    feature1 = TheLayer.GetFeature(num) 
    lake_feature = feature1.ExportToJson(as_object=True)
    lake_cart = Path_LkIsl_ShpFile(lake_feature['geometry']['coordinates']) 
    EB_id = lake_feature['properties']['EBhex']
    clim_dat,rlat,rlon,time,metadata,txtfname = Read_CORDEX_V2(nc_path)

    vname, m1, m2, dexp, m3, m4, m5, m6, drange_orignial = metadata 
    if rprt_tme == True:
        b = clock.time()    
    if plots == True:
        Preview_Lake(lake_cart)        
        print 'Island aware Area(km^2)=', Area_Lake_and_Islands(lake_cart),         
        print ', No. xy lake boundary points=',len(lake_cart.vertices)

    lake_rprj = Path_Reproj(lake_cart,False)             
    sub_clim,sub_rlat,sub_rlon = TrimToLake(lake_rprj,clim_dat[0,:,:],rlat,
                                            rlon,off = 3, show = False) 
    tcheck1= clock.time()
    weight_mask = Pixel_Weights(lake_rprj,sub_clim,sub_rlat,sub_rlon)

    np.save('Lakes/Weights/'+EB_id[2:],weight_mask)
    print 'Generated ',EB_id[2:],' in ',clock.time()-a
    return



def Weighted_Mean_3D(weight_mask,all_time_clim,chatty):
    '''Purpose    - Reads in the 2D weight mask (from Pixel_Weights
    function) and the trimmed data from TrimToLake and returns a
    weighted mean. This should be iterated for each time-step.
    This was fixed (27/04/2015) by BAL. There was a bug in the .dot
    approach that was giving bad values. This requires re-calculating
    as multi-pixel data is probably wrong...
    Input    - weight_mask: Pixel weights
             - all_time_clim: the climate data (lat lon subset)
             from TrimToLake3D().
             - chatty: a true or false statement, if true, it will 
               print some info on the weighting process.
    Output   - val_out: the weighted mean value
    '''
    tvals = weight_mask != 0.   # Indx where weights exist
    if len(weight_mask[tvals]) < 1:
        print 'Error, no shape cover identified in weighted mask.'
        return
    if chatty:
        print 'Shape covers',len(weight_mask[tvals]),' pixels'
        print 'Actual weight values are :',weight_mask[tvals]
        print 'Pix weight cum. sum (should end as 1.0):',np.cumsum(weight_mask[tvals])
    wm = weight_mask[tvals].flatten()  # Get a small array of the weight pixels only
    tmp_series = []
    
    # Roll the axis so that from a loop I can pull back all time dimension for each pixl
    ctemp = np.rollaxis(all_time_clim,2) 
    ctemp = np.rollaxis(ctemp,2)
    
    for n,pix in enumerate(ctemp[tvals,:]):  # For each row (timeseries) of pixel data
        tmp_series.append(wm[n]*pix)         # append it
    tmp_series = np.array(tmp_series)
    time_series = np.sum(tmp_series,axis=0)
    return time_series


def Weighted_Mean_3D_old(weight_mask,all_time_clim,chatty):
    '''
    I THINK THIS IS BUGGY - so I created this old version to track 
    how it used to look. New changes in the non-old version above.
    Purpose    - Reads in the 2D weight mask (from Pixel_Weights
    function) and the trimmed data from TrimToLake and returns a
    weighted mean. This should be iterated for each time-step.
    Input    - weight_mask: Pixel weights
             - all_time_clim: the climate data (lat lon subset)
             - chatty: a true or false statement, if true, it will 
               print some info on the weighting process.
    Output   - val_out: the weighted mean value
    '''
    aaa = (weight_mask > 0.000).flatten()   # Indx where weights exist
    if (len(aaa) == 0):
        print 'Error: no lake cover identified! :('
    if (chatty == True):
        print 'Lake covers',len(weight_mask[aaa]),' pixels'
        print 'Actual weight values are :',weight_mask[aaa]
        print 'Cum. sum of pixel weights (should end as 1.0):',cumsum(weight_mask[aaa])
    dim1, dim2, dim3 = all_time_clim.shape
    d = all_time_clim.reshape((dim1, dim2 * dim3))[:, aaa] ## n_time x n_pixels
    wm = weight_mask.flatten()[aaa]
    val_out = np.dot(d, wm.reshape(len(wm), 1)) ## n_time x n_pixels dot n_pixels x 1
    return val_out

def Write_HDF(f,EB_id,tlist,offset,area):
    '''
    Simple writing using H5py module in HDF5 format.
    (This should be used for the final program)
    '''
    dset = f.create_dataset(EB_id,(len(tlist),),dtype='f')
    dset.attrs['Area']=area
    dset.attrs['Offset']=offset
    dset[...]=tlist
    return

def write_hdf_catchment(fw,EB_id,tseries):
    '''
    Simple writing using H5py module in HDF5 format.
    (This should be used for the final program) to
    write the catchment time series data to hdf file.
    '''
    dim1 = len(tseries)
    dset = fw.create_dataset(EB_id,(dim1,),dtype='f')
    dset[...]=tseries
    return

def write_metadata_csv(mfname,meta_list):
    '''
    Writes a list of metadata (meta_list) to a csv file one
    line at a time, appending to the file. This means I dont
    need to keep all the data in memory at once.
    requires > import csv
    '''
    with open(mfname, 'a') as csvfile:
        metawriter = csv.writer(csvfile,delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        metawriter.writerow(meta_list)
    return
    
def Write_HDF_weights(fw,EB_id,weights):
    '''
    Simple writing using H5py module in HDF5 format.
    (This should be used for the final program)
    '''
    dim1,dim2 = np.shape(weights)
    dset = fw.create_dataset(EB_id,(dim1,dim2,),dtype='f')
    dset[...]=weights
    return




# VERSION 1: DEPRECIATED!
#
#    PRIMARY FUNCTIONS ALL COMBINED INTO A SINGLE FUNCTION BELOW
#             READY TO BE RUN USING MUTLITHREADDING MODULE
#
def OLD_MT_Means_Over_Lake(lake_num):
    ''' This function combines all the routines of the function file,
    to create time series (or images) of the individual lakes, with 
    data overlain. It outputs to a path created within the function,
    unuqie to each lake. This is ready to be run on parallalel 
    processors (non-shared memory): each instance loads its own
    version of the data. Takes about 7min per lake (for one EUR-11
    CORDEX file), and produces a 14kb text file (5kb in gz format):
    speed is on an approx 3GHz machine.

    WARNING - THIS VERSION OF THE CODE IS DEPRECIATED. PLEASE USE 
    MT_Means_Over_Lake  (which appears at the top of this file).

    '''
    num = lake_num
    # Load data - Note Multiprocessing does not use shared memory
    lake_data = 'Lakes/comsat_fetch.geojson'
    lake_geometry, lake_id, lake_name,lake_path,lake_altitude,\
        feno_lake_id = Read_Lakes(lake_data)
    clim_dat,rlat,rlon,time,drange,dexp = Tmp_CORDEX_Read()  
    orog = Height_CORDEX()
    
    # Create Lake shapes and subset the CORDEX data
    lake_cart= Path_Make(lake_path[num][0][0][:])
    lake_rprj = Path_Reproj(lake_cart,False)             
    sub_clim,sub_rlat,sub_rlon= TrimToLake(lake_rprj,\
        clim_dat[0,:,:],rlat,rlon,off = 3, show = False)
    weight_mask = Pixel_Weights(lake_rprj,sub_clim,sub_rlat,sub_rlon)
   # print ' area=',Poly_Area2D(EqArea(lake_rprj.vertices)),'km^2'

    # The below two lines are where the height adjustment comes from
    sub_orog,sub_rlat,sub_rlon = TrimToLake(lake_rprj,orog[:,:],\
        rlat,rlon,off = 3, show = False)
    hght,offset = Orographic_Adjustment(weight_mask,sub_orog,\
        lake_altitude[num],clim_dat,chatty=False)
   # print 'Weighted mean height of lake pixels and lake (m):',
   # print hght,' vs ',lake_altitude[num]
   # print 'Height diffrence between lake and EUR-11 pixels (m):',
   # print lake_altitude[num]-hght
   # print 'Offset to be applied (in Kelvin):',offset

    # Plot out the lake shape for evaluation
    fig_lake = Plot_LakeAndData_Save(lake_rprj,clim_dat[0,:,:],rlat,\
                rlon,zoom=0.10)
    out = 'Outputs/'
    figname = Gen_FileName(out,lake_name[num],'fig','.pdf')
    fig_lake.savefig(figname,dpi=72)

    # Process the time-series
    tlist =[]
    #for t in xrange(30):   # Shorter, test time
    for t in xrange(len(time)):
        sub_clim,sub_rlat,sub_rlon = TrimToLake(lake_rprj,\
            clim_dat[t,:,:],rlat,rlon,off = 3, show = False)
        final_val = Weighted_Mean(weight_mask,sub_clim,chatty=False)
        tlist.append(final_val)
        #print 'Timestep:',t, '  Weighted temperature =',
        #print '%6.2f'%(final_val) #-272.15),'Deg C'
    
    # Finally, generate a file name and output the time-series data
    fnm_head=Gen_FName_Head(clim_dat.long_name)
    out = 'Outputs/'+fnm_head
    hcreate = 'Height offset = '+ '%6.2f'%offset+'  Data = '+\
        clim_dat.long_name+', time range of file = '+drange+\
        ' Model = '+dexp
    txtfname = Gen_FileName(out,lake_name[num],drange,'.txt.gz')       
    np.savetxt(txtfname,tlist,fmt='%7.3f',newline='\n', header=hcreate)
    return