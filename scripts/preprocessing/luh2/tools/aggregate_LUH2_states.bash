#!/bin/bash
# JN 06.12.16, script to aggregate LUH2 states to LUH1 states ("CMIP6 way")
# (julia.nabel@mpimet.mpg.de)

#module load nco

echo "--> aggregate_LUH2_states.bash"

statesFile=${1}
outfolder=${2}

variablesToAggregateFor_gsecd=( secdf secdn range )
variablesToAggregateFor_gcrop=( c3ann c3per c4ann c4per c3nfx )

if [ ! -d ${outfolder} ]; then
    mkdir ${outfolder}
fi

cd ${outfolder}

#--- aggregate gother: primf and primn
cdo -s -chname,primf,gothr -add -selvar,primf ${statesFile} -selvar,primn ${statesFile} tmp_aggregated_gothr_selyear.nc
ncatted -a long_name,gothr,d,, tmp_aggregated_gothr_selyear.nc

#--- aggregate gpast: pastr and urban (TR decided urban should go into pastr!)
cdo -s -chname,pastr,gpast -add -selvar,pastr ${statesFile} -selvar,urban ${statesFile} tmp_aggregated_gpast_selyear.nc
ncatted -a long_name,gpast,d,, tmp_aggregated_gpast_selyear.nc

#--- aggregate gsecd: apparently what jsbach would define to be deserts are often defined as rangelands in LUH2v2h, therefore we assign rangelands to secondary lands 
cdoCommand=""
index=0
(( stopIndex = ${#variablesToAggregateFor_gsecd[@]} - 1 ))
while [ $index -lt ${stopIndex} ]; do
    cdoCommand="${cdoCommand} -add -selvar,${variablesToAggregateFor_gsecd[$index]} ${statesFile}"
   (( index = index + 1 ))
done
cdoCommand="${cdoCommand} -selvar,${variablesToAggregateFor_gsecd[$index]} ${statesFile}"
cdo -s -chname,${variablesToAggregateFor_gsecd[0]},gsecd ${cdoCommand} tmp_aggregated_gsecd_selyear.nc
ncatted -a long_name,gsecd,d,, tmp_aggregated_gsecd_selyear.nc

#--- aggregate gcrop
cdoCommand=""
index=0
(( stopIndex = ${#variablesToAggregateFor_gcrop[@]} - 1 ))
while [ $index -lt ${stopIndex} ]; do
    cdoCommand="${cdoCommand} -add -selvar,${variablesToAggregateFor_gcrop[$index]} ${statesFile}"
   (( index = index + 1 ))
done
cdoCommand="${cdoCommand} -selvar,${variablesToAggregateFor_gcrop[$index]} ${statesFile}"
cdo -s -chname,${variablesToAggregateFor_gcrop[0]},gcrop ${cdoCommand} tmp_aggregated_gcrop_selyear.nc
ncatted -a long_name,gcrop,d,, tmp_aggregated_gcrop_selyear.nc

# join the 4 files
cdo -s merge tmp_aggregated_?????_selyear.nc tmp_all_states.nc
rm tmp_aggregated_*_selyear.nc

# JSBACH and the JSBACH scripting expect one file per year
cdo -s -splityear tmp_all_states.nc LUH_states_
rm tmp_all_states.nc
