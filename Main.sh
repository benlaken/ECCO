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
echo "Commencing ECCO workflow..."

wget_file="/uio/kant/geo-metos-u1/blaken/datadisk/ECCO/CORDEX/wget_EUR-11.KNMI_all.sh"
nc_path="/uio/kant/geo-metos-u1/blaken/datadisk/ECCO/CORDEX/Data_ECCO/"



echo $wget_file

#bash $action_folder$wget_script

ls $nc_path"*.nc" -G

pushd $nc_path
pwd
for f in *; do
  echo "File -> $f"
done

popd
pwd