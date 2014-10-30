from ECCO_functions import *
from multiprocessing import Process,Pool,cpu_count
import time
import math


nc_path = '/uio/kant/geo-metos-u1/blaken/Downloads/'+\
          'tas_EUR-11_ICHEC-EC-EARTH_rcp85_r3i1p1_DMI-'+\
          'HIRHAM5_v1_day_20060101-20101231.nc'
lake_data = 'Lakes/largest100.geojson'
out_path = 'Outputs'



if __name__ == "__main__":
    start = time.time()
    processes= [Process(target=MT_Means_Over_Lake, args=(nc_path,lake_data,\
                    i,out_path,), kwargs={'plots':False,'rprt_tme':False,})\
                                    for i in range(24)] # Run x lakes 

    cpus = cpu_count() -2  # CPU's to work on (leave some free...)?
    p_groups= int(math.ceil(float(len(processes))/float(cpus))) # How many groups?
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