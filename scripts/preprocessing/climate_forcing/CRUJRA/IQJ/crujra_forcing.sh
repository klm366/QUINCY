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


#-----------------------------------------------------------------------------
# The script reads 6-hourly CRU-JRA v3.5 forcing over 1901- 2024 (longwave radiation,
# shortwave radiation, tmax, tmin, precipitation, wind speed, specific humidity),
# computes daily statistics, inserts Feb 29 for leap years by duplicating Feb 28,
# reorders longitude from [-180,180] to [0,360], merges all variables into a
# single file, and remap conservatively onto the ICON grid.
#----------------------------------------------------------------------------
#
###############################################################################
### Batch Queuing System is SLURM
#SBATCH --job-name=crujra_forcing
#SBATCH --output=crujra_forcing.o%j
#SBATCH --error=crujra_forcing.o%j
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=5120 
#SBATCH --ntasks=1
#SBATCH --account mj0143
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --time=07:30:00



module load cdo nco

DataPath="/work/mj0143/icon_land_data_pool/trendy/2026/raw-data-20260702/crujra3"
# Grid
#grid_def_r2b4=/pool/data/ICON/grids/public/mpim/0049/icon_grid_0049_R02B04_G.nc
grid_def_r2b5=/pool/data/ICON/grids/public/edzw/icon_grid_0030_R02B05_G.nc

outdir="./crujra_climate"
tmpdir="./tmp_crujra_climate"


mkdir -p "$outdir" "$tmpdir"
# Set cdo command to be silent, create no history and use double precision
cdo="cdo -s --no_history -b F64"
#
# Process CRUJRA data 1901-2025
for year in $(seq 2010 2025); do
  echo "Processing year $year"

  if (( year % 4 == 0 )); then leap=1; else leap=0; fi

  # 1. LONGWAVE
  ${cdo} daymean \
    "${DataPath}/crujra.v4.10.5d.dlwrf.${year}.365d.noc.nc" \
    "${tmpdir}/longwave_365_${year}.nc"
  echo "done longwave"

  # 2. SHORTWAVE
  ${cdo} daymean \
    "${DataPath}/../Radiation/Processed/tswrf_v15_${year}.nc" \
    "${tmpdir}/shortwave_365_${year}.nc"
  echo "done shortwave"

  # 3. TMAX
  ${cdo} subc,273.15 \
    -daymax "${DataPath}/crujra.v4.10.5d.tmax.${year}.365d.noc.nc" \
    "${tmpdir}/tmax_365_${year}.nc"
  echo "done tmax"

  # 4. TMIN
  ${cdo} subc,273.15 \
    -daymin "${DataPath}/crujra.v4.10.5d.tmin.${year}.365d.noc.nc" \
    "${tmpdir}/tmin_365_${year}.nc"
  echo "done tmin"

  # 5. PRECIP
  ${cdo} divc,86400.0 \
    -daysum "${DataPath}/crujra.v4.10.5d.pre.${year}.365d.noc.nc" \
    "${tmpdir}/precip_365_${year}.nc"
  echo "done precip"

  # 6. WIND SPEED
  ${cdo} daymean \
    -expr,'wspeed=sqrt(ugrd*ugrd+vgrd*vgrd)' \
    -merge \
      "${DataPath}/crujra.v4.10.5d.ugrd.${year}.365d.noc.nc" \
      "${DataPath}/crujra.v4.10.5d.vgrd.${year}.365d.noc.nc" \
    "${tmpdir}/wspeed_365_${year}.nc"
  echo "done wspeed"

  # 7. QAIR
  ${cdo} daymean \
    "${DataPath}/crujra.v4.10.5d.spfh.${year}.365d.noc.nc" \
    "${tmpdir}/qair_365_${year}.nc"
  echo "done qair"

  # 8. Add Feb 29 for leap years (duplicate Feb 28 = timestep 59)
  vars=(longwave shortwave tmax tmin precip wspeed qair)

  for var in "${vars[@]}"; do
    src="${tmpdir}/${var}_365_${year}.nc"
    dst="${tmpdir}/${var}_final_${year}.nc"
    if (( leap == 1 )); then
      ${cdo} cat \
        -seltimestep,1/59  "$src" \
        -seltimestep,59    "$src" \
        -seltimestep,60/365 "$src" \
        "$dst"
      ${cdo} settaxis,${year}-01-01,00:00:00,1day "$dst" "${dst%.nc}_fixed.nc"
      mv "${dst%.nc}_fixed.nc" "$dst"
    else
      ${cdo} settaxis,${year}-01-01,00:00:00,1day "$src" "$dst"
    fi
  done

  # 9. Rename variables + set metadata
  ${cdo} chname,dlwrf,longwave  "${tmpdir}/longwave_final_${year}.nc" "${tmpdir}/lw_r_${year}.nc"
  ${cdo} chname,tswrf,shortwave "${tmpdir}/shortwave_final_${year}.nc" "${tmpdir}/sw_r_${year}.nc"
  cp "${tmpdir}/tmax_final_${year}.nc"    "${tmpdir}/tmax_r_${year}.nc"
  cp "${tmpdir}/tmin_final_${year}.nc"    "${tmpdir}/tmin_r_${year}.nc"
  ${cdo} chname,pre,precip      "${tmpdir}/precip_final_${year}.nc"   "${tmpdir}/precip_r_${year}.nc"
  cp "${tmpdir}/wspeed_final_${year}.nc" "${tmpdir}/wspeed_r_${year}.nc"
  ${cdo} chname,spfh,qair       "${tmpdir}/qair_final_${year}.nc"     "${tmpdir}/qair_r_${year}.nc"

  ncatted -O \
    -a standard_name,longwave,o,c,"surface_downwelling_longwave_flux_in_air" \
    -a long_name,longwave,o,c,"Surface Downwelling Longwave Radiation" \
    -a units,longwave,o,c,"W m-2" \
    "${tmpdir}/lw_r_${year}.nc"
  ncatted -O \
    -a standard_name,shortwave,o,c,"surface_downwelling_shortwave_flux_in_air" \
    -a long_name,shortwave,o,c,"Surface Downwelling Shortwave Radiation" \
    -a units,shortwave,o,c,"W m-2" \
    -a _FillValue,shortwave,o,d,9.96920996838687e+36 \
    -a missing_value,shortwave,o,d,9.96920996838687e+36 \
    "${tmpdir}/sw_r_${year}.nc"
  ncatted -O \
    -a standard_name,tmax,o,c,"air_temperature" \
    -a long_name,tmax,o,c,"Daily Maximum Near-Surface Air Temperature" \
    -a units,tmax,o,c,"degC" \
    "${tmpdir}/tmax_r_${year}.nc"
  ncatted -O \
    -a standard_name,tmin,o,c,"air_temperature" \
    -a long_name,tmin,o,c,"Daily Minimum Near-Surface Air Temperature" \
    -a units,tmin,o,c,"degC" \
    "${tmpdir}/tmin_r_${year}.nc"
  ncatted -O \
    -a standard_name,precip,o,c,"precipitation_flux" \
    -a long_name,precip,o,c,"Precipitation" \
    -a units,precip,o,c,"kg m-2 s-1" \
    "${tmpdir}/precip_r_${year}.nc"
  ncatted -O \
    -a standard_name,wspeed,o,c,"wind_speed" \
    -a long_name,wspeed,o,c,"Near-Surface Wind Speed" \
    -a units,wspeed,o,c,"m s-1" \
    "${tmpdir}/wspeed_r_${year}.nc"
  ncatted -O \
    -a standard_name,qair,o,c,"specific_humidity" \
    -a long_name,qair,o,c,"Near-Surface Specific Humidity" \
    -a units,qair,o,c,"kg kg-1" \
    "${tmpdir}/qair_r_${year}.nc"

  # 10. Reorder longitude [-180,180] -> [0,360] and merge
  for var in lw sw tmax tmin precip wspeed qair; do
    ${cdo} sellonlatbox,0,360,-90,90 \
      "${tmpdir}/${var}_r_${year}.nc" \
      "${tmpdir}/${var}_lon360_${year}.nc"
  done
 
  ${cdo} merge \
    "${tmpdir}/lw_lon360_${year}.nc" \
    "${tmpdir}/sw_lon360_${year}.nc" \
    "${tmpdir}/tmax_lon360_${year}.nc" \
    "${tmpdir}/tmin_lon360_${year}.nc" \
    "${tmpdir}/precip_lon360_${year}.nc" \
    "${tmpdir}/wspeed_lon360_${year}.nc" \
    "${tmpdir}/qair_lon360_${year}.nc" \
    "${tmpdir}/climate_${year}_grid.nc"

  # 11. Remap to ICON grid
  for grid_tag in r2b4 r2b5; do
    # if [[ "$grid_tag" == "r2b4" ]]; then
    #   grid_def="$grid_def_r2b4"
    #   grid_outdir="${outdir}/r2b4"
    #   institute_tag="mpim"
    # else
    #   grid_def="$grid_def_r2b5"
    #   grid_outdir="${outdir}/r2b5"
    #   institute_tag="edzw"
    # fi
    grid_def="$grid_def_r2b5"
    grid_outdir="${outdir}/r2b5"
    institute_tag="edzw"
    mkdir -p "$grid_outdir"

    final_file="${grid_outdir}/climate_crujra_${grid_tag}_${institute_tag}_${year}.nc"
    
    # Change relative time axis to absolute
    ${cdo} -a remapcon,"$grid_def" \
      "${tmpdir}/climate_${year}_grid.nc" \
      "$final_file"   
    ${cdo} -a setcalendar,proleptic_gregorian \
      "$final_file" \
      "${tmpdir}/climate_crujra_${grid_tag}_${institute_tag}_${year}_calendar.nc"

    mv "${tmpdir}/climate_crujra_${grid_tag}_${institute_tag}_${year}_calendar.nc" \
      "$final_file"

    # Check for negative precipitation values in the final forcing file.
    precip_min="$(${cdo} output -fldmin -timmin -selname,precip "$final_file")"

    if awk "BEGIN {exit !($precip_min < 0)}"; then
      echo "ERROR: Negative precipitation detected in $final_file (min=${precip_min})." >&2
      exit 1
    fi

    # Remove NCO-generated history attributes and overwrite the global contact attribute
    ncatted -O -h \
      -a history,global,d,, \
      -a history_of_appended_files,global,d,, \
      -a title,global,o,c,"Daily climate forcing based on CRUJRA v3.5" \
      -a contact,global,o,c,"szaehle@bgc-jena.mpg.de" \
      "${grid_outdir}/climate_crujra_${grid_tag}_${institute_tag}_${year}.nc"   
   
    echo "Year $year remapped -> ${grid_outdir}/climate_crujra_${grid_tag}_${institute_tag}_${year}.nc"
  done

  # Cleanup
  rm -f "${tmpdir}"/*_${year}.nc

done 

# Fill 1850-1900 
# #echo "Filling years 1850-1900."

# #declare -A src_year=(
#   [1900]=1921 [1899]=1919 [1898]=1918 [1897]=1917 [1896]=1916
#   [1895]=1915 [1894]=1914 [1893]=1913 [1892]=1912 [1891]=1911
#   [1890]=1910 [1889]=1909 [1888]=1908 [1887]=1907 [1886]=1906
#   [1885]=1905 [1884]=1904 [1883]=1903 [1882]=1902 [1881]=1901
#   [1880]=1920 [1879]=1919 [1878]=1918 [1877]=1917 [1876]=1916
#   [1875]=1915 [1874]=1914 [1873]=1913 [1872]=1912 [1871]=1911
#   [1870]=1910 [1869]=1909 [1868]=1908 [1867]=1907 [1866]=1906
#   [1865]=1905 [1864]=1904 [1863]=1903 [1862]=1902 [1861]=1901
#   [1860]=1920 [1859]=1919 [1858]=1918 [1857]=1917 [1856]=1916
#   [1855]=1915 [1854]=1914 [1853]=1913 [1852]=1912 [1851]=1911
#   [1850]=1910
# )

# #for grid_tag in r2b4 r2b5; do
#   if [[ "$grid_tag" == "r2b4" ]]; then
#     grid_outdir="${outdir}/r2b4"
#     institute_tag="mpim"
#   else
#     grid_outdir="${outdir}/r2b5"
#     institute_tag="edzw"
#   fi

#   for y in $(seq 1900 -1 1850); do
#     src=${src_year[$y]}
#     echo "Copying $src -> $y"
#     ${cdo} setyear,${y} \
#       "${grid_outdir}/climate_crujra_${grid_tag}_${institute_tag}_${src}.nc" \
#       "${grid_outdir}/climate_crujra_${grid_tag}_${institute_tag}_${y}.nc"
#   done

# done

# echo "All years processed."

