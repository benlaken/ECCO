from ECCO_functions import *
from multiprocessing import Process,Pool,cpu_count
import time
import math

start = time.time()

nc_path = '/uio/kant/geo-metos-u1/blaken/Downloads/'+\
          'tas_EUR-11_ICHEC-EC-EARTH_rcp85_r3i1p1_DMI-'+\
          'HIRHAM5_v1_day_20060101-20101231.nc'

lake_data = 'Lakes/largest100.geojson'
lake_num = 5
out_path = 'Outputs'


#MT_Means_Over_Lake(nc_path,lake_data,\
#                4,out_path,plots=True,rprt_tme=True)


#for n in range(5):
#
#    MT_Means_Over_Lake(nc_path,lake_data,\
#                    n,out_path,plots=False,rprt_tme=False)




# Nb. Below, i, an incrementing integer, represents the lake number
processes= [Process(target=MT_Means_Over_Lake, args=(nc_path,lake_data,\
                    i,out_path,), kwargs={'plots':False,'rprt_tme':False,})\
                                    for i in range(16)]


cpus = cpu_count() -1  # number of CPU's and leave some free...

cpus = 4

p_groups= int(math.ceil(float(len(processes))/float(cpus))) # how many groups 


srng = 0
frng = cpus
for n in xrange(p_groups):
    [p.start() for p in processes[srng:frng]]
    [p.join() for p in processes[srng:frng]]
    print 'loop',n,' of ', p_groups,'  group:',srng,':',frng
    srng+= cpus
    frng+= cpus  # increment start and finish by the number of CPUs


#for p in processes:
#    p.start()
#for p in processes:
#    p.join()           #Kill zombies


end = time.time()
print 'Finished processing in %6.2f min.'%((end-start)/60.)

