#!/bin/ksh
#-------------------------------------------------------------------------------
# Generation of soil parameter arrays for JSBACH in different resolutions 
# 
# Based on 
# soil_para.com  Vs. 3.0 - October 2011 - Stefan Hagemann, MPI-M
#
#
# Veronika Gayler, August 2013
#-------------------------------------------------------------------------------
set -e

echam=true
icon=true

#----------------------------
pool=/pool/land/jsbach/INI/soil_parameters/0.5deg
[[ -d 0.5deg ]] || mkdir 0.5deg
#----------------------------
# input data on 0.5 deg grid
#----------------------------

cd 0.5deg

varlist="soildepth hydcond matrixpot bclapp volporosity \
         volfc volpwp heatcapacity heatcond psi_lambda"

# grid description for 0.5 deg grid
cat > 0.5grid <<EOF
gridtype=lonlat
xsize=720
ysize=360
xfirst=-179.75
xinc=0.5
yfirst=89.75
yinc=-0.5
EOF

# convert data from service format to netcdf
for var in ${varlist}; do
  if [[ ${var} = soildepth ]]; then
    ifile=bin_${var}_ext.srv
  else
    ifile=ber_05_${var}.srv
  fi
  [[ -f ${ifile} ]] || ln -s ${pool}/${ifile} .
  cdo -s -f nc setgrid,0.5grid -setvar,${var} ${ifile} ${var}_0.5grid.tmp
  cdo -s setmissval,-9999 ${var}_0.5grid.tmp ${var}_0.5grid.nc

  case ${var} in
    soildepth)
      longname="soil depth until bedrock"
      unit="m"
      reference="?" ;;
    hydcond)
      # Unit conversion cm --> m
      cdo -s mulc,0.01 ${var}_0.5grid.nc ${var}_0.5grid.tmp
      mv ${var}_0.5grid.tmp ${var}_0.5grid.nc
      longname="saturated hydraulic conductivity"
      unit="m/s"
      reference="Beringer et al. (2001)" ;;
    matrixpot)
      longname="saturated matrix potential"
      unit="m"
      reference="Beringer et al. (2001)" ;;
    bclapp)
      longname="Clapp and Hornberger (1978) exponent b"
      unit=" "
      reference="Beringer et al. (2001)" ;;
    volporosity)
      longname="volumetric soil porosity"
      unit="m/m"
      reference="Beringer et al. (2001)" ;;
    volfc)
      longname="volumetric soil field capacity"
      unit="m/m"
      reference="Patterson (1990)" ;;
    volpwp)
      longname="volumetric wilting point"
      unit="m/m"
      reference="Patterson (1990)" ;;
    heatcapacity)
      longname="heat capacity of dry soil"
      unit="J m-3 K-1"
      reference="Beringer et al. (2001)" ;;
    heatcond)
      longname="heat conductivity of dry soil"
      unit="J m-1 s-1 K-1"
      reference="Beringer et al. (2001)" ;;
    psi_lambda)
      longname="soil pore size distribution index"
      unit=" "
      reference="William and Ahuja (2003)" ;;
  esac

  code=$(cdo -s showcode ${ifile} | tr -d " ")       # code number in srv format

  ncatted -h -a long_name,${var},o,c,"${longname}" \
             -a units,${var},o,c,"${unit}" \
             -a code,${var},o,c,${code} \
             -a references,${var},o,c,"${reference}"   ${var}_0.5grid.nc
done
rm -f *.tmp

# soil available water wcap (T. Stacke)
#    global attributes copied to variable attributes 

[[ -f wcap_lsp3.nc ]] || ln -s ${pool}/wcap_lsp3.nc .
cdo -s -f nc -b f32 \
   setmissval,-9999 -setmisstoc,-9999 -setmissval,0. wcap_lsp3.nc wcap_0.5grid.nc
ncatted -h -a contact,wcap,o,c,"tobias.stacke@zmaw.de"  \
  -a references,wcap,o,c,"Based on an optimized root depth data set (Kleidon 2004)" \
  -a title,wcap,o,c,"Distribution of maximum soil water capacity in LSP3" \
     wcap_0.5grid.nc

[[ -f wcap_lsp3_dynveg.nc ]] || ln -s ${pool}/wcap_lsp3_dynveg.nc .
cdo -s -f nc -b f32 \
   setmissval,-9999 -setmisstoc,-9999 -setmissval,0. -setvar,wcap_dynveg \
   wcap_lsp3_dynveg.nc wcap_dynveg_0.5grid.nc
ncatted -h -a contact,wcap_dynveg,o,c,"tobias.stacke@zmaw.de"  \
  -a references,wcap_dynveg,o,c,"Based on an optimized root depth data set (Kleidon 2004)" \
  -a title,wcap_dynveg,o,c,"Distribution of maximum soil water capacity in LSP3" \
  -a comment,wcap_dynveg,o,c,"Desert type values are substituted with semi desert shrubs for dynveg simulations" \
  wcap_dynveg_0.5grid.nc


# interpolation to Gaussian grid
#--------------------------------

if [[ ${echam} = true ]]; then
#for res in 21 31 42 63 85 106 127 159 255; do
for res in 63; do
  mergelist=""
  for var in ${varlist} wcap wcap_dynveg; do
    if [[ ! -f t${res}weights ]]; then 
      cdo -s gencon,t${res}grid ${var}_0.5grid.nc t${res}weights
    fi
    cdo -s setctomiss,0 ${var}_0.5grid.nc ${var}_0.5grid.tmp  # mask out zeros on glaciers 
    cdo -s remap,t${res}grid,t${res}weights ${var}_0.5grid.tmp ${var}_T${res}.nc
    mergelist="${mergelist} ${var}_T${res}.nc"
  done
  cdo -s merge ${mergelist} ../soil_parameters_T${res}.nc

  # define global attributes
  date=$(date)
  history="${date}: generated with $0"
  source="base data on 0.5deg grid provided by Stefan Hagemann, MPI-M"
  ncatted -h -a history,global,o,c,"${history}" \
             -a source,global,o,c,"${source}" ../soil_parameters_T${res}.nc
  rm ${mergelist}
done
fi

# interpolation to ICON grid
#--------------------------------

if [[ ${icon} = true ]]; then
for grid in r2b4_amip r2b6; do

  case ${grid} in
    r2b4_amip)
      gridfile=/pool/data/ICON/grids/private/${grid}/${grid}.nc ;;
    r2b6)
      gridfile=/pool/data/ICON/grids/private/leonidas/gcc_globalGrids/Earth_Icos_0039km_etopo05_landlakes.nc ;;
  esac

  mergelist=""
  for var in ${varlist} wcap wcap_dynveg; do
    if [[ ! -f ${grid}_weights ]]; then 
      cdo -s gencon,${gridfile} ${var}_0.5grid.nc ${grid}_weights
    fi
    cdo -s setctomiss,0 ${var}_0.5grid.nc ${var}_0.5grid.tmp  # mask out zeros on glaciers 
    cdo -s fillmiss -setmissval,-9999 ${var}_0.5grid.tmp ${var}_0.5grid.tmp1
    cdo -s remap,${gridfile},${grid}_weights ${var}_0.5grid.tmp1 ${var}_${grid}.nc
    mergelist="${mergelist} ${var}_${grid}.nc"
  done
  cdo -s merge ${mergelist} ../soil_parameters_${grid}.nc

  # define global attributes
  date=$(date)
  history="${date}: generated with $0"
  source="base data on 0.5deg grid provided by Stefan Hagemann, MPI-M"
  uuidOfHGrid=$(ncdump -h ${gridfile} | grep ':uuidOfHGrid = ' | cut -f2 -d'"')
  ncatted -h -a history,global,o,c,"${history}" \
             -a uuidOfHGrid,global,o,c,${uuidOfHGrid} \
             -a source,global,o,c,"${source}" ../soil_parameters_${grid}.nc
  rm ${mergelist}

done
fi

# clean up
rm -f *weights 0.5grid *.tmp *.tmp1
