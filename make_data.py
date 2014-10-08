# Main python program to create time-series of Lake data
import numpy as np
import time
from multiprocessing import Process,Pool,cpu_count
import ECCO_functions as ecco

start_time = time.time()         # Note the start time of the program
print '\n','Start of ECCO data creation software','\n'
print 'Number of cores:',cpu_count()

print 'Note, check system status to see if processing is occuring', 
print 'there will be no progress written to terminal for a while...'
# Set range below to be equal to the number of lakes in the lake file
processes= [Process(target=ecco.MT_Means_Over_Lake, args=(i,)) for i in range(10)]
#processes= [Process(target=ecco.MT_Means_Over_Lake, args=(i,)) for i in range(5)]
for p in processes:
    p.start()                                 # Start the processes
for p in processes:
    p.join()                                  # Kill zombies

print '\n'+'Lake processing complete'
elapsed = (time.time() - start_time)/60.
print 'Completed in','%6.f'%elapsed,'miniutes'
print 'End of program'