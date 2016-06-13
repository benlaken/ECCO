import ECCO_functions_v2 as ECCO

ncfile = 'CORDEX/tas_EUR-11_ICHEC-EC-EARTH_rcp45_r1i1p1_KNMI-RACMO22E_v1_day_20960101-21001231.nc'

ECCO.Catchment_Weights_Meta(nc_path=ncfile, sbar=True)