#!/bin/bash

files=`ls /datascope/pradal-esm/${1}/history/*ocean_month.nc`
#files=`ls /datascope/pradal-esm/${1}/history/*ocean_bling.nc`

# Ocean Month:
#var="temp_xflux_adv temp_yflux_adv temp_zflux_adv temp temp_tendency neutral_temp temp_submeso temp_nonlocal_KPP sw_heat temp_vdiffuse_impl rho_dzt"
var='mld'

# Ocean BLING: 
#var="dic"


for v in $var
  do
  echo $v
  i=0
  for f in $files
    do
    num=`printf %03d $i`
    ncra -v $v $f ~/control_aredi800/data/var${num}_tmp.nc
    let i=i+1
  done
  files_ann=`ls ~/control_aredi800/data/var*_tmp.nc`
  ncrcat $files_ann ~/control_aredi800/data/$v.nc
  rm var*_tmp.nc
done

