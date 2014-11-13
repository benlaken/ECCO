#!/bin/bash
##############################################################################
# Main script to execute the ECCO workflow
#
# Version: 1
# Will call a wget script to download CORDEX files, held on ESGF system account
# (this step requires OpenID credentials on ESGF), in netcdf format.
# Then calls ECCO_Main.py executed with a nice command, to run through the list
# of downloaded netcdf files.
# Then erases the NetCDF files.
# 
# Benjamin A. Laken (blaken@geo.uio.no)
###############################################################################

pth_vann="/mn/vann/metos-d1/blaken/ECCO" 

wget_folder=$pth_vann/"CORDEX"
nc_path=$pth_vann/"CORDEX/Data_CORDEX"

#cd $pth_vann
# cd -
# pwd
# read


### DOWNLOAD THE DATA
pushd $wget_folder
pwd
bash $wget_folder/"wget_EUR-11.KNMI_all.sh"    # RUN A WGET SCRIPT
popd


### PROCESS THE DATA
nice python ECCO_Main.py  # run main python script from working dir

### REMOVE THE NETCDF FILES
pushd $nc_path
pwd
#read
echo "Erasing the following files:"
for i in $(ls *.nc); do echo $i; done
#for i in $(ls *.nc); do rm $i; done      # REMOVE THE .NC FILES
#popd