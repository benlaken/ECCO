import pandas as pd
import pg8000
from sqlalchemy import create_engine

engine = create_engine('postgresql+pg8000://user:password@vm-srv-finstad.vm.ntnu.no/ecco_biwa_db')

simheadsql = 'select * from sim limit 20;'
julymixdepthsql = '''
select 
       mm.ya0 as y0, 
       mm.yb1 as y1, 
       mm.rcma as rcm, 
       t1.ebint, 
       t1.JulyMixDepth
from 
     (select floor(sim_id / 2e7) as mm, 
     sim_id % 2e7 as ebint, avg(pcd07) as JulyMixDepth 
     from sim 
     group by sim_id) t1, 
     mm   
where mm.mm = t1.mm 
order by ebint, mm.mm
limit 30;
'''

simhead = pd.read_sql(simheadsql, engine)
julymixdepth = pd.read_sql(julymixdepthsql, engine)



