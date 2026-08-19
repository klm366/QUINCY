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

input_data_n="/work/mj0143/icon_land_data_pool/trendy/2026/raw_data_20260804/"
input_path_p="/work/mj0143/icon_land_data_pool/preprocessing/np_depositions/p_dep_files/"
grid_out="/work/mj0143/b301069/data/NMIP3_input/IQ_grid/landsea_grid.nc"
target_grid="/pool/data/ICON/grids/public/edzw/icon_grid_0030_R02B05_G.nc"
outdir="./Ndep_trendy2026_R2B5_MONTHLY_on_luh_01"
fraction_files="/work/mj0143/icon_land_data_pool/trendy/2026/preprocessesed-01/fraction_files/R2B5/0030/"
states="/work/mj0143/icon_land_data_pool/trendy/2026/raw-data-20260702/luh2/states_LUH2_GCB2026.nc"
temp_path="${outdir}"/temp

mkdir -p "${outdir}"
mkdir -p "${temp_path}"

# Define years to process depositions
start_year=1850
n_dep_end_year=2025
p_dep_end_year=2019
nhx_in="${input_data_n}/Ndep_NHx.nc"
noy_in="${input_data_n}/Ndep_NOy.nc"
mask="luh" # frac or luh
dep_unit="kg/m2/s"

module load nco
# Creating files with land from fraction files

if [[ "$mask" == "frac" ]]; then
  # Deriving land from fraction files

  for year in $(seq 1850 2025); do
    cdo -expr,'land= notsea*fract_land' "${fraction_files}"/"bc_land_frac_iq_13pfts_${year}.nc" "${temp_path}/land_frac_file.nc"
    cdo -selvar,land "${temp_path}/land_frac_file.nc" "${temp_path}"/"tmp_land_bc_land_frac_iq_13pfts_${year}.nc"
  done
  cdo -mergetime "${temp_path}"/tmp_land_bc_land_frac_iq_13pfts_*.nc "${temp_path}/land_frac_1850-2025_r2b5.nc"
  rm "${outdir}/tmp_land_bc_land_frac_iq_13pfts_*.nc"
# remaping to ndep grid

  cdo -remapcon,${nhx_in} "${temp_path}/land_frac_1850-2025_r2b5.nc" "${outdir}/land_frac_1850-2025.nc" 
  cdo -mul ${nhx_in} "${outdir}/land_frac_1850-2025.nc" "${outdir}/nhx_land.nc"
  cdo -mul ${noy_in} "${outdir}/land_frac_1850-2025.nc" "${outdir}/noy_land.nc"
else 

  # processing land created from states file
  cdo -selyear,1850/2024 ${states} "${outdir}"/"states_1850_2025.nc"
  cdo -expr,'land= primf + primn + secdf + secdn + urban + c3ann + c4ann + c3per + c4per + c3nfx + pastr + range' \
   "${outdir}"/"states_1850_2025.nc" "${outdir}"/"land_states_1850_2025_grid.nc"   

  cdo -remapcon,"${noy_in}" "${outdir}"/"land_states_1850_2025_grid.nc" \
  "${outdir}"/"land_states_1850-2025.nc" 

  cdo -mul ${nhx_in} "${outdir}"/"land_states_1850-2025.nc" "${outdir}/nhx_land.nc"
  cdo -mul ${noy_in} "${outdir}"/"land_states_1850-2025.nc" "${outdir}/noy_land.nc"
fi


# Convert units:
# gN m-2 land -> kg m-2 s-1
# divide by 1000 and by seconds per year
# sec_year=31536000
cdo mulc,3.170979e-11 "${outdir}/nhx_land.nc" "${outdir}/nhx_kg_s.nc"
cdo mulc,3.170979e-11 "${outdir}/noy_land.nc" "${outdir}/noy_kg_s.nc"

cdo -O -setmissval,nan "${outdir}/nhx_kg_s.nc" "${temp_path}/nhx_kg_s_mvn.nc"
cdo -O -setmisstoc,0 "${temp_path}/nhx_kg_s_mvn.nc" "${temp_path}/nhx_kg_s_clean.nc"

cdo -O -setmissval,nan "${outdir}/noy_kg_s.nc" "${temp_path}/noy_kg_s_mvn.nc"
cdo -O -setmisstoc,0 "${temp_path}/noy_kg_s_mvn.nc" "${temp_path}/noy_kg_s_clean.nc"

ncks -O -h -v Ndep_NHx "${temp_path}/nhx_kg_s_clean.nc" ${temp_path}/foo2.nc

ncrename -O -h -v Ndep_NHx,NHx_deposition ${temp_path}/foo2.nc
ncatted -h -a long_name,NHx_deposition,o,c,"NHx_deposition" \
        -h -a standard_name,,d,c,"" \
        -h -a original_name,,d,c,"" \
        -h -a bounds,,d,c,"" \
        -h -a units,NHx_deposition,o,c,"${dep_unit}" ${temp_path}/foo2.nc


ncks -O -h -v Ndep_NOy "${temp_path}/noy_kg_s_clean.nc" ${temp_path}/foo3.nc

ncrename -O -h -v Ndep_NOy,NOy_deposition ${temp_path}/foo3.nc
ncatted -h -a long_name,NOy_deposition,o,c,"NOy_deposition" \
        -h -a standard_name,,d,c,"" \
        -h -a original_name,,d,c,"" \
        -h -a bounds,,d,c,"" \
        -h -a units,NOy_deposition,o,c,"${dep_unit}" ${temp_path}/foo3.nc

cdo splityear ${temp_path}/foo2.nc "${temp_path}/nhx_"
cdo splityear ${temp_path}/foo3.nc "${temp_path}/noy_"

# 6) Loop over years and build the annual output
for year in $(seq ${start_year} ${n_dep_end_year}); do
  src_year=${year}
  if [ "${year}" -eq 2025 ]; then
    src_year=2024
  fi
  # make sure time dimension exists with 12 months
  cdo merge "${temp_path}/nhx_${year}.nc" "${temp_path}/noy_${year}.nc" "${temp_path}/bc_land_ndep_quincy_${year}_0.25.nc"
  cdo --no_history -remapcon,${target_grid} "${temp_path}/bc_land_ndep_quincy_${year}_0.25.nc" "${outdir}/bc_land_ndep_quincy_${year}.nc"

done



# Loop through each year from start_year to end_year
year=${start_year}
while [ ${year} -le ${n_dep_end_year} ]; do
  # Define input and output paths
  ndep_in="${outdir}/bc_land_ndep_quincy_${year}.nc"
  npdep_out_pre="${temp_path}/bc_land_npdep_quincy_${year}_pre.nc"
  npdep_out="${outdir}/bc_land_npdep_quincy_${year}.nc"
  if [ ${year} -gt ${p_dep_end_year} ]; then
    pdep_in="${input_path_p}/bc_land_npdep_quincy_2019.nc"
  else
    pdep_in="${input_path_p}/bc_land_npdep_quincy_${year}.nc"
  fi
  # Select the 'pdep' variable from the files containing phosphorus depositions
  cdo --no_history -selvar,pdep "${pdep_in}" ${temp_path}/pdep_${year}_temp.nc
  # Remap the selected Pdep file to the target grid defined by grid_def_file
  cdo --no_history -remapcon,${target_grid} ${temp_path}/pdep_${year}_temp.nc ${temp_path}/pdep_${year}.nc
  # Merge yearly Ndep and Pdep
  cdo --no_history -merge "${ndep_in}" "${temp_path}/pdep_${year}.nc" "${npdep_out_pre}"
  ncatted -O -h \
  -a references,global,a,c,"; Brahney, J., N. Mahowald, D. S. Ward, A. P. Ballantyne, and J. C. Neff (2015), Is atmospheric phosphorus pollution altering global alpine Lake stoichiometry?, Global Biogeochem. Cycles, 29, 1369–1383, doi:10.1002/2015GB005137." \
  -a variable_id,global,d,, \
  -a title,global,d,,"N input data for NMIP3" \
  "${npdep_out_pre}" "${npdep_out}"
  # Clean up temporary files
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
        cdo --no_history -settaxis,${year}-01-01,00:00:00,1month \
            "${npdep_in}" "${npdep_out}"
    fi
done
# Remove temporary directory
rm -r ${temp_path}

exit 0

