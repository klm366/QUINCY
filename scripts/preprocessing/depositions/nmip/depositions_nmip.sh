#!/bin/bash 

# ICON-Land
#
# ---------------------------------------
# Copyright (C) 2013-2026, MPI-M, MPI-BGC
#
# Contact: icon-model.org
# Authors: AUTHORS.md
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------

#!/usr/bin/env bash
#--------------------------------------------------------------------
### Batch Queuing System is SLURM
#SBATCH --output=depositions_nmip.o%j
#SBATCH --error=depositions_nmip.o%j
#SBATCH --account=mj0143
#SBATCH --partition=interactive
#SBATCH --mem-per-cpu=3000
#SBATCH --exclusive
#SBATCH --ntasks=1

# script is based on ncl script from Cheng Gong (/work/mj0143/b301072/ncl_scripts/pre-processing_ndep.ncl)

module load nco
module unload cdo
module load cdo/2.5.0-gcc-11.2.0
# set ncatted command to create no history
ncatted="ncatted -h"
# set cdo command to be silent, create no history and use double precision
cdo="cdo -s --no_history -b 64"

scriptName="depositions_nmip.sh"
#----------------------------------------------------------------------

# path and name of the target output folder
outdir=./YOUROUTPUTFOLDER
input_data_freq="monthly" # monthly/yearly

# path to the files
quincy_pool="/work/mj0143/icon_land_data_pool"
if [[ "$input_data_freq" == "monthly" ]]; then
  input_data_n="${quincy_pool}/trendy/2026/raw_data_20260804/"
elif [[ "$input_data_freq" == "yearly" ]]; then
  input_data_n="${quincy_pool}/trendy/2026/raw-data-20260702/ndep/"
else
  echo "The input frequency is not specified"
fi
input_path_p="${quincy_pool}/preprocessing/np_depositions/p_dep_files/"
fraction_files="${quincy_pool}/trendy/2026/preprocessesed-01/fraction_files/R2B5/0030/"
states="${quincy_pool}/trendy/2026/raw-data-20260702/luh2/states_LUH2_GCB2026.nc"

# grid
target_grid="/pool/data/ICON/grids/public/edzw/icon_grid_0030_R02B05_G.nc"
#target_grid=/pool/data/ICON/grids/public/mpim/0049/icon_grid_0049_R02B04_G.nc

#----------------------------------------------------------------------

temp_path="${outdir}"/temp

mkdir -p "${outdir}"
mkdir -p "${temp_path}"

# Define years to process depositions
start_year=1850  # start year from which on to use nmip deposition data 
n_dep_end_year=2025 # end year until which on to use nmip deposition data 
p_dep_end_year=2019
nhx_in="${input_data_n}/Ndep_NHx_monthly.nc"
noy_in="${input_data_n}/Ndep_NOy_monthly.nc"
mask="luh" # frac or luh, luh is preferred
dep_unit="kg/m2/s"

# Some global attributes to be set
contact="szaehle@bgc-jena.mpg.de"
reference="; Brahney, J., N. Mahowald, D. S. Ward, A. P. Ballantyne, and J. C. Neff (2015), Is atmospheric phosphorus pollution altering global alpine Lake stoichiometry?, Global Biogeochem. Cycles, 29, 1369–1383, doi:10.1002/2015GB005137."
title="N input data for NMIP3"
#----------------------------------------------------------------------

# create files with land from fraction files
if [[ "$mask" == "frac" ]]; then
  # derive land from fraction files
  for year in $(seq 1850 2025); do
    ${cdo} -expr,'land= notsea*fract_land' "${fraction_files}"/"bc_land_frac_iq_13pfts_${year}.nc" "${temp_path}/land_frac_file.nc"
    ${cdo} -selvar,land "${temp_path}/land_frac_file.nc" "${temp_path}"/"tmp_land_bc_land_frac_iq_13pfts_${year}.nc"
  done
  # merge land variable from fraction files
  ${cdo} -mergetime "${temp_path}"/tmp_land_bc_land_frac_iq_13pfts_*.nc "${temp_path}/land_frac_1850-2025_r2b5.nc"
  rm "${outdir}/tmp_land_bc_land_frac_iq_13pfts_*.nc"
  # remap to n depositions grid
  ${cdo} -remapcon,${nhx_in} "${temp_path}/land_frac_1850-2025_r2b5.nc" "${outdir}/land_frac_1850-2025.nc" 
  # weight grid-cell deposition by the effective land fraction
  # to retain only the deposition contribution over land 
  ${cdo} -mul ${nhx_in} "${outdir}/land_frac_1850-2025.nc" "${outdir}/nhx_land.nc"
  ${cdo} -mul ${noy_in} "${outdir}/land_frac_1850-2025.nc" "${outdir}/noy_land.nc"
else 
  # derive land from states file (luh data)
  ${cdo} -selyear,1850/2024 ${states} "${outdir}"/"states_1850_2025.nc"
  ${cdo} -expr,'land= primf + primn + secdf + secdn + urban + c3ann + c4ann + c3per + c4per + c3nfx + pastr + range' \
   "${outdir}"/"states_1850_2025.nc" "${outdir}"/"land_states_1850_2025_grid.nc"   
  # remap to n depositions grid
  ${cdo} -remapcon,"${noy_in}" "${outdir}"/"land_states_1850_2025_grid.nc" \
  "${outdir}"/"land_states_1850-2025.nc" 
  # weight grid-cell deposition by the effective land fraction
  # to retain only the deposition contribution over land 
  ${cdo} -mul ${nhx_in} "${outdir}"/"land_states_1850-2025.nc" "${outdir}/nhx_land.nc"
  ${cdo} -mul ${noy_in} "${outdir}"/"land_states_1850-2025.nc" "${outdir}/noy_land.nc"
fi

# convert units:
# gN m-2 land -> kg m-2 s-1
if [[ "$input_data_freq" == "monthly" ]]; then
  # divide by 1000 by seconds per month
  # 60 * 60 * 24 * days per month
  ${cdo} -mulc,0.001 -divc,86400 -divdpm "${outdir}/nhx_land.nc" "${outdir}/nhx_kg_s.nc"
  ${cdo} -mulc,0.001 -divc,86400 -divdpm "${outdir}/noy_land.nc" "${outdir}/noy_kg_s.nc"
else 
  # divide by 1000 and by seconds per year
  # sec_year=31536000 
  ${cdo} mulc,3.170979e-11 "${outdir}/nhx_land.nc" "${outdir}/nhx_kg_s.nc"
  ${cdo} mulc,3.170979e-11 "${outdir}/noy_land.nc" "${outdir}/noy_kg_s.nc"
fi

# change missing values to 0
${cdo} -O -setmissval,nan "${outdir}/nhx_kg_s.nc" "${temp_path}/nhx_kg_s_mvn.nc"
${cdo} -O -setmisstoc,0 "${temp_path}/nhx_kg_s_mvn.nc" "${temp_path}/nhx_kg_s_clean.nc"

${cdo} -O -setmissval,nan "${outdir}/noy_kg_s.nc" "${temp_path}/noy_kg_s_mvn.nc"
${cdo} -O -setmisstoc,0 "${temp_path}/noy_kg_s_mvn.nc" "${temp_path}/noy_kg_s_clean.nc"

# select only the NHx variable, rename it, and set attributes
ncks -O -h -v Ndep_NHx_monthly "${temp_path}/nhx_kg_s_clean.nc" ${temp_path}/foo2.nc
ncrename -O -h -v Ndep_NHx_monthly,NHx_deposition ${temp_path}/foo2.nc
${ncatted} -a long_name,NHx_deposition,o,c,"NHx_deposition" \
         -a standard_name,,d,c,"" \
         -a original_name,,d,c,"" \
         -a bounds,,d,c,"" \
         -a units,NHx_deposition,o,c,"${dep_unit}" ${temp_path}/foo2.nc

# select only the NOy variable, rename it, and set attributes
ncks -O -h -v Ndep_NOy_monthly "${temp_path}/noy_kg_s_clean.nc" ${temp_path}/foo3.nc
ncrename -O -h -v Ndep_NOy_monthly,NOy_deposition ${temp_path}/foo3.nc
${ncatted} -a long_name,NOy_deposition,o,c,"NOy_deposition" \
         -a standard_name,,d,c,"" \
         -a original_name,,d,c,"" \
         -a bounds,,d,c,"" \
         -a units,NOy_deposition,o,c,"${dep_unit}" ${temp_path}/foo3.nc

# split N depositions into years
${cdo} splityear ${temp_path}/foo2.nc "${temp_path}/nhx_"
${cdo} splityear ${temp_path}/foo3.nc "${temp_path}/noy_"

# loop over years and build the annual output
for year in $(seq ${start_year} ${n_dep_end_year}); do
  # crete an exeption for 2025 when the data are missing
  if [ "${year}" -eq 2025 ]; then
    last_year=2025
    ${cdo} -O settaxis,2025-01-01,00:00:00,1mon "${outdir}/bc_land_ndep_quincy_2024.nc" "${outdir}/bc_land_ndep_quincy_${last_year}.nc" 
  fi
  if [[ "$input_data_freq" == "yearly" ]]; then
    ${cdo} -O settaxis,${year}-01-01,00:00:00,1month -duplicate,12 \
        "${temp_path}/nhx_${year}.nc" "${temp_path}/nhx_${year}_12months.nc"
    ${cdo} -O settaxis,${year}-01-01,00:00:00,1month -duplicate,12 \
    "${temp_path}/noy_${year}.nc" "${temp_path}/noy_${year}_12months.nc"
    # merge yearly files with NOy and NHx
    ${cdo} merge "${temp_path}/nhx_${year}_12months.nc" "${temp_path}/noy_${year}_12months.nc" \
    "${temp_path}/bc_land_ndep_quincy_${year}_0.25.nc"
  elif [[ "$input_data_freq" == "monthly" ]]; then
    # merge yearly files with NOy and NHx
    ${cdo} merge "${temp_path}/nhx_${year}.nc" "${temp_path}/noy_${year}.nc" "${temp_path}/bc_land_ndep_quincy_${year}_0.25.nc"
  else 
    echo "The input frequency is not specified"
  fi
  ${cdo} --no_history -remapcon,${target_grid} "${temp_path}/bc_land_ndep_quincy_${year}_0.25.nc" \
  "${outdir}/bc_land_ndep_quincy_${year}.nc"
done

# loop through each year from start_year to end_year
year=${start_year}
while [ ${year} -le ${n_dep_end_year} ]; do
  # define input and output paths
  ndep_in="${outdir}/bc_land_ndep_quincy_${year}.nc"
  npdep_out_pre="${temp_path}/bc_land_npdep_quincy_${year}_pre.nc"
  npdep_out="${outdir}/bc_land_npdep_quincy_${year}.nc"
  # p depositions data are available till 2019
  if [ ${year} -gt ${p_dep_end_year} ]; then
    pdep_in="${input_path_p}/bc_land_npdep_quincy_2019.nc"
  else
    pdep_in="${input_path_p}/bc_land_npdep_quincy_${year}.nc"
  fi
  # select the 'pdep' variable from the files containing phosphorus depositions
  ${cdo} -selvar,pdep "${pdep_in}" ${temp_path}/pdep_${year}_temp.nc
  # remap the selected Pdep file to the target grid defined by grid_def_file
  ${cdo} -remapcon,${target_grid} ${temp_path}/pdep_${year}_temp.nc ${temp_path}/pdep_${year}.nc
  # merge yearly Ndep and Pdep
  ${cdo} -merge "${ndep_in}" "${temp_path}/pdep_${year}.nc" "${npdep_out_pre}"
  
  # set further global attributes
  ${ncatted} -O \
  -a references,global,a,c,"${reference}" \
  -a variable_id,global,d,, \
  -a title,global,d,,"${title}" \
  -a contact,global,o,c,"${contact}" \
  "${npdep_out_pre}" "${npdep_out}"
  # clean up temporary files
  rm -f "${temp_path}/pdep_${year}_temp.nc" "${temp_path}/pdep_${year}.nc"
  year=$((year + 1))
done
echo "Processing complete. Yearly NP-deposition files are located in ${output_path}."

# Use the 1850 NP-deposition file for years 1700–1850
for year in $(seq 1700 1850); do
    npdep_in="${outdir}/bc_land_npdep_quincy_1850.nc"
    npdep_out="${outdir}/bc_land_npdep_quincy_${year}.nc"

    if [ "${year}" -eq 1850 ]; then
        # The 1850 file already has the correct year
        cp -f "${npdep_in}" "${npdep_out}"
    else
        # Change the time axis to the requested year
        ${cdo} -settaxis,${year}-01-01,00:00:00,1month \
            "${npdep_in}" "${npdep_out}"
    fi
done
# Remove temporary directory
rm -r ${temp_path}

exit 0

