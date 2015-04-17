from __future__ import print_function,division,generators
import numpy as np
import pandas as pd
import sys
import h5py
import datetime as dt
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

def calc_density(yvals,histxvals,lake_eg,DynamicScale=False):
    ''' Create a normalized time-series density plot. The data is binned to 
    monthly time-steps. yvals = list of y-values, histxvals = list of integers
    for x-axis (not dt values)
    '''
    histout = np.histogram2d(yvals,histxvals,bins=[50,120])
    return histout

def center_bins(original):
    ''' Returns a list of the bin centers
    '''
    bin_center =[]
    bin_size = original[0] - original[1]
    for n in range(len(original)-1):
        bin_center.append(original[n] + bin_size/2.)
    return bin_center

def fromID_Gen2DInput(lakeids,kind='abs'):
    ''' Takes a list of lake ID codes (meant for single_lakes). 
    Also a keyword specifying the type of data to return (kind = 'absoloute'
     or 'current_anomaly' or 'historical_anomaly').
    The diffrent kinds reffer to the abosoloute values of the CORDEX data, a
    de-seasonalised anomaly, or a deasezonalised anomaly downsampled to 
    monthly resolution. Returned values are a list of the y-values, an 
    integer list used to subscript the y-values and a real list of the 
    datatime objects for the x-values. The primary purpouse of this function
    is to prepare the input for np.histogram2d()
    '''
    yvals=[]
    histxvals=[]
    xvals=[]
    print("")
    for n,lakeid in enumerate(lakeids):
        #status(n,len(lakeids)-1, bar_length=40)
        print("Progress:{0:6.2f}%, Lake:{1}, Index:{2}".format(n/len(lakeids)*100,lakeid,n))
        #(n=None,tmp_dts=dts_comb,f1=f1,f2=f2)
        tmp_lake = read_fast(n=lakeid,tmp_dts=dts_comb,f1=f1,f2=f2)
        if kind == 'abs':     
            yvals += list(tmp_lake.values.flatten())
            histxvals += range(len(tmp_lake))
            xvals += tmp_lake.index      
        elif kind == 'anom':
            yvals += list(Seasonality(tmp_lake,monthly=False).values.flatten())
            histxvals +=range(len(tmp_lake.values)) 
            xvals +=tmp_lake.index
        elif kind == 'hist_anom':
            tmp = Seasonality_Historical(tmp_lake,monthly=False,historical=True)
            yvals += list(tmp.values.flatten())
            histxvals += range(len(tmp))
            xvals += tmp.index
        else:
            print('\n Error: {0} kind not reccognised.'.format(kind))
            print('Enter: abs(absoloute), anom(current_anomaly),or hist_anom (historical_anomaly)')
            return None,None,None
    return yvals,histxvals,xvals

def get_CORDEXmeta_from_fnm(tmpfnm):
    '''From the input path/file of CORDEX output processed in the ECCO project, a dic
    of metadata about the specific file are produced.
    '''
    output={} 
    tmpfnm = tmpfnm.split('/')[-1]        # strip the file name from the directories
    fnm_meta = tmpfnm.split('_')         # break apart the file name into metadata
    #print(fnm_meta)
    output['variable'] = fnm_meta[1]
    output['RCM'] = fnm_meta[3]
    output['scenario'] = fnm_meta[4]
    output['ensemble'] = fnm_meta[5]
    output['GCM'] = fnm_meta[6]
    output['dates'] = fnm_meta[9].split('.h5')[0]  # trim off the hdf file ending too
    return output

def Gen_df_list(file1,file2,idlist):
    ''' Generate a list of dictionaries (pd df objects)
    from an input list of lake id values and two input files.
    '''
    df_list=[]
    print('Generating list of lake dictionaries:')
    for n,idval in enumerate(idlist):
        df_list.append(read_fast(n=idval,tmp_dts=dts_comb,f1=f1,f2=f2))
        #(n=None,tmp_dts=dts_comb,f1=f1,f2=f2)
        status(n,len(idlist)-1, bar_length=40)
    return df_list

def gen_timelist(timestring):
    '''This is expecting a string of dates (as is obtained directly from the dictionary when
    running the function get_CORDEXmeta_from_fnm(), variable Output['dates']. Once input, it
    will return a list of datetime objects at a daily resolution covering the time range.
    '''
    drange=[]
    start_date = timestring.split('-')[0]
    end_date = timestring.split('-')[1]
    #print(start_date,end_date)
    start = dt.datetime.strptime(start_date,'%Y%m%d')
    end = dt.datetime.strptime(end_date,'%Y%m%d')
    for n in range(((end - start).days)+1):
        drange.append((start + dt.timedelta(days=n)).date())
    return drange

def status(current,end_val, bar_length=20):
    ''' Spawns a status bar
    '''
    percent = float(current) / end_val
    hashes = '#' * int(round(percent * bar_length))
    spaces = ' ' * (bar_length - len(hashes))
    sys.stdout.write("\rPercent: [{0}] {1}% ".format(hashes + spaces, int(round(percent * 100))))
    sys.stdout.flush()

def read_fast(n,tmp_dts,f1,f2):
#def read_fast(n=None,tmp_dts,f1=f1,f2=f2):
    '''
    Return a dataframe of 10-years of data for a given lake ID.
    Requires the file links, dt index, etc to be pre-produced.
    this is an optimized version of read_data().
    '''
    n = str(n)
    tmp_data = list(f1[n]) + list(f2[n])  
    df = pd.DataFrame(tmp_data[:],index=tmp_dts,columns=[n])
    return df

def read_historical(n,tmp_dts,hf1,hf2):
#def read_historical(n=None,tmp_dts=histdts_comb,hf1=hf1,hf2=hf2):
    '''
    Return a dataframe of 10-years of data for a given lake ID.
    Requires the file links, dt index, etc to be pre-produced.
    this is an optimized version of read_data()
    This is specifically to obtain Historical data (which should
    be fixed to the decade of 1970-1980). Need to read both data
    in question (with read_fast) and historical data (with this 
    function) as calculation of delta should be relative to a 
    historical period.
    '''
    n = str(n)
    tmp_data = list(hf1[n]) + list(hf2[n])  
    df = pd.DataFrame(tmp_data[:],index=tmp_dts,columns=[n])
    return df

def Seasonality_Historical(datin, seazon=False, monthly=False,historical=True):
    ''' From a pd df object, with dt as the index, this calculates seasonality.
    The anomalies (deseazonalised data) are then returned, also as a pd df.
    The returned values are at a Monthly sampling if monthly=True. Else, they
    are at daily timescales.
    If seazon = True, the seazonal cycle is returned, not the anomalies.
    Does this for individual pd dictionaryies (ie. one lake at a time).
    '''
    desz = []
    seazon_mean = []
    tmpdic = pd.DataFrame(data=np.zeros(len(datin.values),dtype=float),
                            index = datin.index)
    if historical == True:
        for key in datin:
            hist_data = read_historical(n=key,tmp_dts=histdts_comb,hf1=hf1,hf2=hf2)
        DOY_index_hist = [dnum.timetuple().tm_yday for dnum in hist_data.index]
        DOY_index = [dnum.timetuple().tm_yday for dnum in datin.index] 
    if historical == False:
        DOY_index = [dnum.timetuple().tm_yday for dnum in datin.index]

    for yday in range(min(DOY_index),max(DOY_index)+1):
        if historical == False:
            mask = [doy == yday for doy in DOY_index]
            seazon_mean.append(np.mean(datin[mask])) 
            tmpdic[mask]= datin[mask] - np.mean(datin[mask])
        if historical == True:
            mask = [doy == yday for doy in DOY_index]
            mask_hist = [doy == yday for doy in DOY_index_hist]
            seazon_mean.append(np.mean(hist_data[mask_hist]))
            tmpdic[mask]= datin[mask] - np.mean(hist_data[mask_hist])
    if seazon == True:
        return seazon_mean
    if monthly == False:
        return tmpdic
    elif monthly == True:
        return tmpdic.asfreq('M')


if __name__ == "__main__":
    #fpth = '/uio/kant/geo-metos-u1/blaken/datadisk/ECCO/Outputs/'
    fpth = 'Outputs/tmp_data/'
    file1 = 'Abel_outputs/Lakes_tas_EUR-11_ICHEC-EC-EARTH_historical_r3i1p1_DMI-HIRHAM5_v1_day_19710101-19751231.h5'
    file2 = 'Abel_outputs/Lakes_tas_EUR-11_ICHEC-EC-EARTH_historical_r3i1p1_DMI-HIRHAM5_v1_day_19760101-19801231.h5'

    #file1  = 'Abel_outputs/Lakes_tas_EUR-11_ICHEC-EC-EARTH_historical_r3i1p1_DMI-HIRHAM5_v1_day_20010101-20051231.h5'
    #file2  = 'Abel_outputs/Lakes_tas_EUR-11_ICHEC-EC-EARTH_rcp45_r3i1p1_DMI-HIRHAM5_v1_day_20060101-20101231.h5'

    #file1  = 'Abel_outputs/Lakes_tas_EUR-11_ICHEC-EC-EARTH_rcp45_r3i1p1_DMI-HIRHAM5_v1_day_20310101-20351231.h5'
    #file2  = 'Abel_outputs/Lakes_tas_EUR-11_ICHEC-EC-EARTH_rcp45_r3i1p1_DMI-HIRHAM5_v1_day_20360101-20401231.h5'

    #file1  = 'Abel_outputs/Lakes_tas_EUR-11_ICHEC-EC-EARTH_rcp45_r3i1p1_DMI-HIRHAM5_v1_day_20610101-20651231.h5'
    #file2  = 'Abel_outputs/Lakes_tas_EUR-11_ICHEC-EC-EARTH_rcp45_r3i1p1_DMI-HIRHAM5_v1_day_20660101-20701231.h5'

    #file1  = 'Abel_outputs/Lakes_tas_EUR-11_ICHEC-EC-EARTH_rcp45_r3i1p1_DMI-HIRHAM5_v1_day_20910101-20951231.h5'
    #file2  = 'Abel_outputs/Lakes_tas_EUR-11_ICHEC-EC-EARTH_rcp45_r3i1p1_DMI-HIRHAM5_v1_day_20960101-21001231.h5'

    meta = pd.read_csv('Metadata/Lake_Stats.csv')              # Load Lake Metadata as a pandasd dataframe object 
    meta.index = meta.Lake_ID                                  # Index the metadata by lake ID's
    single_lakes = np.load('Metadata/single_pixel_lakeid.npy') # Open the lake id metdatada of the single pixel lake id's
    f1 = h5py.File(fpth+file1, "r") 
    f2 = h5py.File(fpth+file2, "r")
    fmeta1 = get_CORDEXmeta_from_fnm(file1)
    fmeta2 = get_CORDEXmeta_from_fnm(file2)
    dts1 = list(gen_timelist(fmeta1['dates']))
    dts2 = list(gen_timelist(fmeta2['dates'])) 
    dts_comb = dts1+ dts2
    # Historical decade (needed to calculate anomalies)
    hfile1 = 'Abel_outputs/Lakes_tas_EUR-11_ICHEC-EC-EARTH_historical_r3i1p1_DMI-HIRHAM5_v1_day_19710101-19751231.h5'
    hfile2 ='Abel_outputs/Lakes_tas_EUR-11_ICHEC-EC-EARTH_historical_r3i1p1_DMI-HIRHAM5_v1_day_19760101-19801231.h5'
    hf1 = h5py.File(fpth+hfile1, "r") 
    hf2 = h5py.File(fpth+hfile2, "r")
    hfmeta1 = get_CORDEXmeta_from_fnm(hfile1)
    hfmeta2 = get_CORDEXmeta_from_fnm(hfile2)  
    histdts1 = list(gen_timelist(hfmeta1['dates']))
    histdts2 = list(gen_timelist(hfmeta2['dates']))
    histdts_comb = histdts1+ histdts2
    # Create list of df objects and anomaly against historical (4hr for 8k)
    lake_dics = Gen_df_list(file1,file2,single_lakes[0:2])
    yvals,histxvals,xvals = fromID_Gen2DInput(single_lakes[0:],kind='hist_anom')  
    zout = calc_density(yvals=yvals,histxvals=histxvals,lake_eg=lake_dics[0])
    outfnm='Hist2D_'+fmeta1['variable']+'_'+fmeta1['dates']+'_'+fmeta1['GCM']
    np.save(fpth+'Sverdrup_outputs/'+outfnm,zout)
    #plt.imshow(zout[0])
    print('Finished')
    #plt.show()