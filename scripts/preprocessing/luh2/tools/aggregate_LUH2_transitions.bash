#!/bin/bash
# JN 06.12.16, script to aggregate LUH2 transitions to LUH1 states ("CMIP6 way")
# (julia.nabel@mpimet.mpg.de)
# SF 18.05.22, add config file to load nco and cdo versions in MAIN
# (stefanie.falk@lmu.de)

echo "--> aggregate_LUH2_transitions.bash"

transFileWithPath=${1}
outfolder=${2}

org_states=( primf primn urban secdf secdn pastr range c3ann c3per c4ann c4per c3nfx )
aggregated_states=( primary primary pasture secondary secondary pasture secondary crop crop crop crop crop )

nonExistingTransitions=( primf_to_secdf primn_to_secdn )
# and there are obviously no transitions to primf or primn

if [ ! -d ${outfolder} ]; then
    mkdir ${outfolder}
fi

cd ${outfolder}

state_from_index=0
fileList=""
while [ ${state_from_index} -lt ${#org_states[@]} ]; do
  state_to_index=0
  state_from=${org_states[$state_from_index]}
  aggregated_state_from=${aggregated_states[$state_from_index]}
  while [ ${state_to_index} -lt ${#org_states[@]} ]; do
    state_to=${org_states[$state_to_index]}
    aggregated_state_to=${aggregated_states[$state_to_index]}

    thisVar="${state_from}_to_${state_to}"
    echo "current var: ${thisVar}"

    if [ ! ${aggregated_state_from} = ${aggregated_state_to} ]; then
      thisTargetVar="${aggregated_state_from}2${aggregated_state_to}"

      # there are no transitions to primf or primn
      if [ ${state_to} = primf ] || [ ${state_to} = primn ]; then
        (( state_to_index = state_to_index + 1 ))
        continue
      fi

      doesNotExist=false
      for nonTrans in ${nonExistingTransitions[@]}; do
        if [ ${nonTrans} = ${thisVar} ]; then
          doesNotExist=true
          break
        fi
      done
      if [ $doesNotExist = true ]; then
        (( state_to_index = state_to_index + 1 ))
        continue
      fi

      echo "... adding ${thisVar} to ${thisTargetVar}"

      # get the data
      cdo -s -selvar,${thisVar} ${transFileWithPath} tmp.nc

      # convert to double (currently they are float but the f90 mapping programm is only implemented for double)
      ${ncap} -s "${thisVar}=double(${thisVar})" tmp.nc tmp1.nc

      thisAggregationFile=tmp_${thisTargetVar}.nc
      if [ -f ${thisAggregationFile} ]; then
        cdo -s -add ${thisAggregationFile} tmp1.nc tmp2.nc
        mv tmp2.nc ${thisAggregationFile}
      else
        cdo -s -chname,${thisVar},${thisTargetVar} tmp1.nc ${thisAggregationFile}
        # remove old not suitable attributes
        ncatted -a standard_name,${thisTargetVar},d,, ${thisAggregationFile}
        fileList="${fileList} ${thisAggregationFile}"
      fi

      rm -f tmp.nc tmp1.nc tmp2.nc
    fi
    (( state_to_index = state_to_index + 1 ))
  done

  (( state_from_index = state_from_index + 1 ))
done

# join the files
cdo -s -merge ${fileList} tmp_all_transitions.nc
rm ${fileList}

# JSBACH and the JSBACH scripting expect one file per year
cdo -s -splityear tmp_all_transitions.nc LUH_transitions_
rm tmp_all_transitions.nc

