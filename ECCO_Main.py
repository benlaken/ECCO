# This script is a multiprocessor wrapper for the ECCO work
# it breaks down all lakes in a given file to n groups of 
# processors, splitting the task, and ensures that the memory does not
# over-fill

# NB. this should be run as: nice ipython ECCO_Main.py

from ECCO_functions_v2 import *
from multiprocessing import Process,Pool,cpu_count
import time
import math
import glob

lake_data = '/uio/kant/geo-metos-u1/blaken/datadisk/ECCO/Lakes/ecco-biwa_lakes_v.0.1.shp'
#lake_data = 'Lakes/largest100.geojson'
out_path = '/uio/kant/geo-metos-u1/blaken/datadisk/ECCO/Outputs'

# Search for all files in the below location with *.nc and generate a list.
nc_folder = '/uio/kant/geo-metos-u1/blaken/datadisk/ECCO/CORDEX/Data_CORDEX/'
nc_list = glob.glob(nc_folder+'*.nc')  



if __name__ == "__main__":
    start = time.time()

    for nfile in nc_list:   # For each netcdf file in a list to processes
        print 'Running on:', nfile[50:]
        processes= [Process(target=MT_Means_Over_Lake, args=(nfile,lake_data,\
                    i,out_path,), kwargs={'plots':False,'rprt_tme':False,})\
                                    for i in range(10)] # Run x lakes 

        cpus = cpu_count() -2  # CPU's to work on (leave some free)
        p_groups= int(math.ceil(float(len(processes))/float(cpus))) # Num groups
        srng = 0
        frng = cpus
        for n in xrange(p_groups):
            [p.start() for p in processes[srng:frng]]
            [p.join() for p in processes[srng:frng]]
            #print 'loop',n,' of ', p_groups,'  group:',srng,':',frng
            srng+= cpus
            frng+= cpus  # increment start and finish by the number of CPUs
            Update_Progress(float(n)/(float(p_groups)-1.))  # A small progress bar...
    fin_time = time.time()
    print '\n','Finished processing in %6.2f min.'%((fin_time-start)/60.)

