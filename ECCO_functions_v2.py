import numpy as np
import sys, time, os, json
from netCDF4 import Dataset  
from matplotlib.path import Path
from matplotlib import cm
import matplotlib.patches as patches
from matplotlib import path
from matplotlib.transforms import Bbox
from math import pi, cos, sin, radians, atan, asin
import mpl_toolkits.basemap.pyproj as pyproj
import matplotlib.pyplot as plt
import os
import time
import time as clock
import osgeo.ogr

# NB the only diffrence between v2 (this file) and the original are that these
# functions are set up to work with the input file for the lakes as a shape file
# whereas the earlier code relied on a GeoJSON format.

def MT_Means_Over_Lake(nc_path, lake_file, lake_num, outputprefix, threeD=True,
                       tt=None,plots = False,rprt_tme=False):
    '''Purpose             This program is the main wrapper to execute the ECCO 
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
    if plots == True:
        Preview_Lake(lake_cart)        
        print 'Island aware Area(km^2)=', Area_Lake_and_Islands(lake_cart),         
        print ', No. xy lake boundary points=',len(lake_cart.vertices)
    lake_rprj = Path_Reproj(lake_cart,False)             
    sub_clim,sub_rlat,sub_rlon = TrimToLake(lake_rprj,clim_dat[0,:,:],rlat,
                                            rlon,off = 3, show = False) 
    weight_mask = Pixel_Weights(lake_rprj,sub_clim,sub_rlat,sub_rlon)
    sub_orog,sub_rlat,sub_rlon = TrimToLake(lake_rprj,orog[:,:],rlat,
                                            rlon,off = 3, show = False)
    hght,offset = Orographic_Adjustment(weight_mask,sub_orog,
                                        lake_altitude,clim_dat,chatty=False)
    if plots == True:
        Show_LakeAndData(lake_rprj,clim_dat[0,:,:],rlat,rlon,zoom=8.)
        Preview_Weights(lake_rprj,weight_mask,sub_rlat,sub_rlon)
    if rprt_tme == True:
        c = clock.time()
    if threeD:
        sub_clim,sub_rlat,sub_rlon = TrimToLake3D(lake_rprj,clim_dat[:,:,:],rlat,rlon,
                                                    off = 3, show = False)
        tlist = Weighted_Mean_3D(weight_mask, sub_clim, chatty=False)
    else:
        tlist =[]
        if tt is None:
            tt = clim_dat.shape[0]
        for t in xrange(tt):
            sub_clim,sub_rlat,sub_rlon = TrimToLake(lake_rprj,clim_dat[t,:,:],
                                                    rlat,rlon,off = 3, show = False)
            if t == 0 :
                final_val = Weighted_Mean(weight_mask,sub_clim,chatty=True)
            else:
                final_val = Weighted_Mean(weight_mask,sub_clim,chatty=False)
            tlist.append(final_val)
            print 'Timestep:',t, ' Weighted tmp. =','%6.2f'%((final_val)-272.15),'C'
        tlist = np.array(tlist)
    if rprt_tme == True:
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
    if rprt_tme == True:
        e = clock.time()
        print '\n'
        print ('%4.2f sec : Read Data:'%(b-a))
        print ('%4.2f sec : Calculated height offset and Weighting Mask'%(c-b))
        print ('%4.2f sec : Weighted time-series'%(d-c))
        print ('%4.2f sec : Folder path and file creation'%(d-c))
        print ('%4.2f sec : Total'%(e-a))
    return





#
#     ROUTINES FOR READING DATA
#
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
    #print(metadata)
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



def Tmp_CORDEX_Read():
    '''Temporay way to read the CORDEX data, while I am testing the 
    software and getting the remote access working still.
    '''
    #cordex_dlist = 'file_list.txt'    # A list of CORDEX data held locally (path and file name)
    #with open(cordex_dlist) as f:
    #    cordex_files = f.readlines()
    #[l.strip('\n\r') for l in cordex_files]
    
    #print cordex_files  # NB having some problem, CORDEX file names are too big for the module to read!
    # I have temporarily renamed tas_EUR-11_ICHEC-EC-EARTH_rcp45_r1i1p1_KNMI-RACMO22E_v1_day_20960101-21001231.nc
    # to tas_20960101-21001231.nc while i am testing...
    #NCfpth = 'Data/CORDEX/tas_20960101-21001231.nc'
    NCfpth = '/uio/kant/geo-metos-u1/blaken/Downloads/tas_EUR-11_ICHEC-EC-EARTH_rcp85_r3i1p1_DMI-HIRHAM5_v1_day_20960101-21001231.nc'
    #print 'Reading file:',NCfpth
    tmp = Dataset(NCfpth,'r')              # Use NetCDF4 to read the file
    tmp.close                              # Close the connection to the file 
    # Gather info about the NetCDF file
    #print 'Experiment ID ',tmp.experiment_id
    #print 'Domain',tmp.CORDEX_domain
    #print 'Driving Experiment:',tmp.driving_experiment
    #print 'Experiment name:',tmp.driving_experiment_name
    #print tmp.filepath()
    #print 'From Institute:',tmp.institute_id

    #for dimobj in tmp.dimensions.values():        # Examine the dimensions
    #    print dimobj
    #for varobj in tmp.variables.values():       # Examine the variable
    #    print varobj
    tas = tmp.variables['tas']
    rlat = tmp.variables['rlat']
    rlon = tmp.variables['rlon']
    time = tmp.variables['time']
    
    drange = NCfpth  # This is a unicode date range of the data stripped from the filename, pass it to the text
    drange = drange[-20:-3]   # files as part of    their naming convention
    dexp = tmp.driving_experiment  # variable to hold experiment info
    return tas,rlat,rlon,time,drange,dexp


#
#  ROUTINES FOR PROCESSING DATA
#

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
    #print 'after hack',lake_poly.codes
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
    #if (lonr < 0.):
    #    lonr=2* pi+lonr            # If the longitude is negative
    return lonr*180/pi, latr*180/pi


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
            #test = lout.intersects_bbox(tmpbb)  # THIS isnt working right with the islands
            if (ben_test == 1):                         # If Lake within the Bbox test will be 1 (no lake, test = 0)
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
                pix_weights[y,x]= float(area_sub)/float(lkarea)   # Fractional area (0-1.0) of lake within a given pixel
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
        patch = patches.PathPatch(lake_rprj, facecolor='#06ebf6', lw=1)
        ax1.add_patch(patch)      # ADD LAKE
        ax1.set_ylabel('Lat. (Deg. N)')
        ax1.set_xlabel('Lon. (Deg. E)')
        ax1.set_title('Preview of trimmed climate data with lake overlaid')
        ax1.imshow(data_sub,interpolation='none', cmap=cm.RdBu,
                   extent=[sub_rlon[0],sub_rlon[-1],sub_rlat[0],sub_rlat[-1]],origin='lower')
        plt.show(fig3)    
        print shape(data_sub),type(data_sub)
    return data_sub,sub_rlat,sub_rlon


def TrimToLake3D(lake_in,Cdat,rlat,rlon,off,show):
    ''' Purpose   - To go from the full CORDEX lat lon array (one
        time slice), to a small np.array subset over the immediate
        pixels surrounding the lake (with a few eitherside). This
        will speed up the execution of the pixel weight calculation
        code.
        
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
        status = "Done...\r\n"
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

def Weighted_Mean_3D(weight_mask,all_time_clim,chatty):
    '''Purpose    - Reads in the 2D weight mask (from Pixel_Weights
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