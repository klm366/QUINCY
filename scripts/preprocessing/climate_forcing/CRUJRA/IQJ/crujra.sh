#!/bin/bash

diri="/work/mj0143/b301069/data/NMIP3_input/climate"
diri2="/work/mj0143/icon_land_data_pool/mpim/0049/climate_forcing/GSWP3/rev002"
grid="/pool/data/ICON/grids/public/mpim/0049/icon_grid_0049_R02B04_G.nc"
outdir="./climate"
tmpdir="./tmp_climate"

mkdir -p "$outdir" "$tmpdir"

for year in $(seq 2020 2024); do
  echo "Processing year $year"

  # ---- Determine leap year ----
  if (( year % 4 == 0 )); then
    leap=1
  else
    leap=0
  fi

  # ---- Determine time axis source file ----
  # (used only to get the time variable in the original script;
  #  CDO will generate correct time axes automatically via settaxis)
  if (( year >= 2020 )); then
    if (( leap == 1 )); then
      ndays=366
    else
      ndays=365
    fi
  fi

  # ================================================================
  # 1. LONGWAVE: 6-hourly -> daily mean
  # ================================================================
  cdo daymean \
    "${diri}/dlwrf/crujra.v3.5d.dlwrf.${year}.365d.noc.nc" \
    "${tmpdir}/longwave_365_${year}.nc"
  echo "done longwave"

  # ================================================================
  # 2. SHORTWAVE: 6-hourly -> daily mean
  # ================================================================
  cdo daymean \
    "${diri}/Radiation/tswrf_gcb2025_${year}.nc" \
    "${tmpdir}/shortwave_365_${year}.nc"
  echo "done shortwave"

  # ================================================================
  # 3. TMAX: 6-hourly -> daily max, convert K -> degC
  # ================================================================
  cdo subc,273.15 \
    -daymax "${diri}/tmax/crujra.v3.5d.tmax.${year}.365d.noc.nc" \
    "${tmpdir}/tmax_365_${year}.nc"
  echo "done tmax"

  # ================================================================
  # 4. TMIN: 6-hourly -> daily min, convert K -> degC
  # ================================================================
  cdo subc,273.15 \
    -daymin "${diri}/tmin/crujra.v3.5d.tmin.${year}.365d.noc.nc" \
    "${tmpdir}/tmin_365_${year}.nc"
  echo "done tmin"

  # ================================================================
  # 5. PRECIP: 6-hourly sum -> daily sum, convert mm/6h -> kg m-2 s-1
  #    Original: sum of 4 timesteps (mm/6h each) = mm/day
  #    Then divided by 86400 -> kg m-2 s-1
  # ================================================================
  cdo divc,86400.0 \
    -daysum "${diri}/pre/crujra.v3.02.5d.pre.${year}.365d.noc.nc" \
    "${tmpdir}/precip_365_${year}.nc"
  echo "done precip"

  # ================================================================
  # 6. WIND SPEED: sqrt(u^2 + v^2), then daily mean
  # ================================================================
  cdo daymean \
    -expr,'wsspeed=sqrt(ugrd*ugrd+vgrd*vgrd)' \
    -merge \
      "${diri}/ugrd/crujra.v3.5d.ugrd.${year}.365d.noc.nc" \
      "${diri}/vgrd/crujra.v3.5d.vgrd.${year}.365d.noc.nc" \
    "${tmpdir}/wsspeed_365_${year}.nc"
  echo "done wsspeed"

  # ================================================================
  # 7. QAIR (specific humidity): 6-hourly -> daily mean
  # ================================================================
  cdo daymean \
    "${diri}/spfh/crujra.v3.5d.spfh.${year}.365d.noc.nc" \
    "${tmpdir}/qair_365_${year}.nc"
  echo "done qair"

  # ================================================================
  # 8. INSERT FEB 29 FOR LEAP YEARS
  #    Duplicate Feb 28 (day 59) and insert as Feb 29 (day 60)
  # ================================================================
  vars=(longwave shortwave tmax tmin precip wsspeed qair)

  for var in "${vars[@]}"; do
    src="${tmpdir}/${var}_365_${year}.nc"
    dst="${tmpdir}/${var}_final_${year}.nc"

    if (( leap == 1 )); then
      # Extract days 1-59, duplicate day 59 as Feb 29, then days 60-365
      cdo cat \
        -selday,1/59   -settaxis,${year}-01-01,00:00:00,1day "$src" \
        -selday,59      -settaxis,${year}-02-29,00:00:00,1day "$src" \
        -selday,60/365  -settaxis,${year}-03-01,00:00:00,1day "$src" \
        "$dst"
    else
      cp "$src" "$dst"
    fi
  done

  # ================================================================
  # 9. RENAME VARIABLES & SET METADATA, MERGE INTO ONE FILE
  # ================================================================

  # Rename variables to target names if needed (CRU-JRA names -> output names)
  cdo chname,dlwrf,longwave   "${tmpdir}/longwave_final_${year}.nc"  "${tmpdir}/lw_r_${year}.nc"
  cdo chname,tswrf,shortwave  "${tmpdir}/shortwave_final_${year}.nc" "${tmpdir}/sw_r_${year}.nc"
  cdo chname,tmax,tmax        "${tmpdir}/tmax_final_${year}.nc"      "${tmpdir}/tmax_r_${year}.nc"
  cdo chname,tmin,tmin        "${tmpdir}/tmin_final_${year}.nc"      "${tmpdir}/tmin_r_${year}.nc"
  cdo chname,pre,precip       "${tmpdir}/precip_final_${year}.nc"    "${tmpdir}/precip_r_${year}.nc"
  # wsspeed and qair already have correct names from expr/spfh rename:
  cdo chname,spfh,qair        "${tmpdir}/qair_final_${year}.nc"      "${tmpdir}/qair_r_${year}.nc"
  cp "${tmpdir}/wsspeed_final_${year}.nc" "${tmpdir}/wsspeed_r_${year}.nc"

  # Set attributes via ncatted (requires NCO) — or skip if not needed for ICON ingest
  # If NCO is available and attributes are required:
  if command -v ncatted &>/dev/null; then
    ncatted -O \
      -a standard_name,longwave,o,c,"surface_downwelling_longwave_flux_in_air" \
      -a long_name,longwave,o,c,"Surface Downwelling Longwave Radiation" \
      -a units,longwave,o,c,"W m-2" \
      "${tmpdir}/lw_r_${year}.nc"

    ncatted -O \
      -a standard_name,shortwave,o,c,"surface_downwelling_shortwave_flux_in_air" \
      -a long_name,shortwave,o,c,"Surface Downwelling Shortwave Radiation" \
      -a units,shortwave,o,c,"W m-2" \
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
      -a standard_name,wsspeed,o,c,"wind_speed" \
      -a long_name,wsspeed,o,c,"Near-Surface Wind Speed" \
      -a units,wsspeed,o,c,"m s-1" \
      "${tmpdir}/wsspeed_r_${year}.nc"

    ncatted -O \
      -a standard_name,qair,o,c,"specific_humidity" \
      -a long_name,qair,o,c,"Near-Surface Specific Humidity" \
      -a units,qair,o,c,"kg kg-1" \
      "${tmpdir}/qair_r_${year}.nc"
  fi

  # ================================================================
  # 10. REORDER LONGITUDE [-180,180] -> [0,360] & MERGE VARIABLES
  # ================================================================
  # CDO sellonlatbox / invertlat won't reorder lon this way, but
  # "cdo sellonlatbox,0,360,-90,90" on a [-180,180] grid will
  # wrap it correctly to [0,360]:
  for var in lw sw tmax tmin precip wsspeed qair; do
    cdo sellonlatbox,0,360,-90,90 \
      "${tmpdir}/${var}_r_${year}.nc" \
      "${tmpdir}/${var}_lon360_${year}.nc"
  done

  # Merge all variables into one grid file
  cdo merge \
    "${tmpdir}/lw_lon360_${year}.nc" \
    "${tmpdir}/sw_lon360_${year}.nc" \
    "${tmpdir}/tmax_lon360_${year}.nc" \
    "${tmpdir}/tmin_lon360_${year}.nc" \
    "${tmpdir}/precip_lon360_${year}.nc" \
    "${tmpdir}/wsspeed_lon360_${year}.nc" \
    "${tmpdir}/qair_lon360_${year}.nc" \
    "${outdir}/climate_${year}_grid.nc"

  # ================================================================
  # 11. REMAP TO ICON GRID (conservative remapping)
  # ================================================================
  cdo remapcon,"$grid" \
    "${outdir}/climate_${year}_grid.nc" \
    "${outdir}/NMIP3_climate_${year}.nc"

  echo "Year $year done -> ${outdir}/NMIP3_climate_${year}.nc"

  # ================================================================
  # 12. PATCH LAND FRACTION FILE (once, using first year's output)
  # ================================================================
  bc_out="./bc_land_frac_iq_13pfts_1850_NMIP.nc"
  bc_src="/pool/data/ICON/grids/public/mpim/0049/land/r0103/bc_land_frac_iq_13pfts_1850.nc"

  if [[ ! -f "$bc_out" ]]; then
    echo "Creating patched land fraction file..."
    # Find cells where sea < 1 AND tmin is missing in remapped output
    # This logic requires Python/NCO since CDO can't do conditional cell patching easily
    python3 - <<'PYEOF'
import netCDF4 as nc
import numpy as np

bc_src  = "/pool/data/ICON/grids/public/mpim/0049/land/r0103/bc_land_frac_iq_13pfts_1850.nc"
clim    = "./climate/NMIP3_climate_2020.nc"   # use first available year
bc_out  = "./bc_land_frac_iq_13pfts_1850_NMIP.nc"

with nc.Dataset(bc_src)  as f_src, \
     nc.Dataset(clim)    as f_clim, \
     nc.Dataset(bc_out, "w") as f_out:

    sea    = f_src["sea"][:]
    notsea = f_src["notsea"][:]
    tmin   = f_clim["tmin"][0, :]          # first timestep

    # cells where sea<1 but tmin is missing -> force to sea
    mask = (sea < 1.0) & tmin.mask
    sea[mask]    = 1.0
    notsea[mask] = 0.0

    # copy dimensions
    for name, dim in f_src.dimensions.items():
        f_out.createDimension(name, len(dim) if not dim.isunlimited() else None)

    # copy all variables, replacing sea/notsea
    for name, var in f_src.variables.items():
        out_var = f_out.createVariable(name, var.datatype, var.dimensions)
        out_var.setncatts({k: var.getncattr(k) for k in var.ncattrs()})
        if name == "sea":
            out_var[:] = sea
        elif name == "notsea":
            out_var[:] = notsea
        else:
            out_var[:] = var[:]

print("Land fraction file written.")
PYEOF
  fi

  # ================================================================
  # 13. CLEANUP INTERMEDIATES (optional)
  # ================================================================
  rm -f "${tmpdir}"/*_${year}.nc

done

echo "All years processed."

