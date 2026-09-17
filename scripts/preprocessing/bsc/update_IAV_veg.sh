#!/bin/bash

#=============================================
#
#    Script update_IAV_veg.ksh
#    ------------------                                                     
#
#    extract 10-day climate/NRT LAI and albedo data
#
#    Usage : update_IAV_veg.ksh [-n LYYYY ] [ -t YYYY ] [-d ECVEG ] [-g GRID ]
#
#    Parameters :  [-n LYYYY ] [ -t YYYY ] [-d ECVEG ] [-g GRID ]
#
#    Externals:                                                       
#    ---------   
#
#    S. Boussetta (ECMWF) July 2023 update generation of IAV veg in the model task from 49R1
#    ----------     
#
#=============================================

set -ex

#Default 
YYYY=2019
LYYYY=clim_2019
#ECVEG=/pa1/VEG_IAV/$CLIMVERSION/${RESOL}${GTYPE}
GRID=lsm_regular.nc
errflg=0
REGULAR=true


USAGE="Usage: update_IAV_veg.ksh [-n LYYYY ] [ -t YYYY ] [-d ECVEG ] [-g GRID ]"
while getopts n:t:d:g: option
do
   case $option in
   n)  LYYYY=$OPTARG;;  # YYYY for a particular year or clim_YYYY  for LAI climate 2010_2019 and CCI YYYY, e.g. clim_2015
   t)  YYYY=$OPTARG;;   # YYYY for a particular year 
   d)  ECVEG=$OPTARG;;  # ECFS path for storing IAV data
   g)  GRID=$OPTARG;;   # GRID file with grid and mask information
   \?)  errflg=1;;
   esac
done
if [ $errflg -ne 0 ]  ; then
  echo $USAGE >&2
  exit 2
fi

OUTDIR=${LYYYY}

##================================================================================
## 1. Prepare paths and filenames
##================================================================================
#STATICDIR="/ec/res5/project/land_surface/ESACCI/scripts/land_cover_veg/static/"
#XDATA_DIR=$XDATA/climate/$CLIMVERSION/${RESOL}${GTYPE}

# Paths for LAI
if [ ${LYYYY:0:4} = clim ]; then
   YYYY=${LYYYY:5:4} # default year for CCI "clim"
   #   FCGLSINPUT="/ec/res5/project/land_surface/CGLS/LAI_v2/corrected/clim/2010_2019/netcdf/"
   FINPUTstart="CGLS_LAI_v2_2010_2019_"
   FINPUTend="15_1km_clim.nc"
   FCGLSINPUT="${LC_LAI_DATA}/CGLS/LAI_v2/corrected/clim/2010_2019/netcdf/"
   TAG="clim_CGLS_LAI"
else
   #   FCGLSINPUT="/ec/res5/project/land_surface/CONFESS/LAI_v1_1km/${LYYYY}/"
   FINPUTstart="CGLS_LAI_v2_${LYYYY}"
   FINPUTend="15_1km.nc"
   FCGLSINPUT="${LC_LAI_DATA}/CGLS/LAI_v2/corrected/netcdf/${LYYYY}/"
   TAG="CGLS_LAI"
   YYYY=${LYYYY}
fi

# Paths for LC
# CCI netcdf original file at 300m 
ccibasename="ESACCI-LC-L4-LCCS-Map-300m-P1Y-"
if [ "${YYYY}" -lt "2016" ]; then
     endname="-v2.0.7cds.nc"
else
     endname="-v2.1.1.nc"
fi
#FCCI="/ec/res5/project/land_surface/ESACCI/data/"${ccibasename}${YYYY}${endname}
FCCI=${LC_LAI_DATA}/ESACCI/${ccibasename}${YYYY}${endname}  #ESACCI-LC-L4-LCCS-Map-300m-P1Y-2019-v2.1.1.nc

python3 -u cut_region.py -r ${GRID} -i ${FCCI} -o ${ccibasename}${YYYY}${endname}_${REGION}.nc

##================================================================================
## 2. Generate Land cover from ESA CCI
##================================================================================

if [[ $REGULAR ]]; then
# Creates an output file called esacci_lc_agg_REG_LL.nc that contains tvl, tvh, cvl, cvh and lsm at lsm_regular.nc resolution
  python3 -u agg_esacci_lc.py -i ${ccibasename}${YYYY}${endname}_${REGION}.nc -c esacci_cwt_ifs.csv -g ${GRID} -s 1 -o ${OUTDIR} -t REG_LL

else
  python3 -u agg_esacci_lc.py -i $FCCI -c esacci_cwt_ifs.csv -g lsm_global -s 0 -o ./ -t ${RESOL}${GTYPE} #Should be changed when gaussian grid option is added
fi

##================================================================================
## 3. Generate LAI from CGLS
##================================================================================

##================================================================================
## 3.1 Generate land cover on LAI grid
##================================================================================

python3 -u cut_region.py -r ${GRID} -i ${FCGLSINPUT}/${FINPUTstart}01${FINPUTend} -o ${FINPUTstart}01${FINPUTend}_${REGION}

#Creates an output file called esacci_lc_agg_ifs_${TAG}.nc that contains tvl, tvh, cvl, cvh and lsm at LAI resolution
python3 -u agg_esacci_lc.py -i ${ccibasename}${YYYY}${endname}_${REGION}.nc  -c esacci_cwt_ifs.csv -g ${FINPUTstart}01${FINPUTend}_${REGION} -s 1 -o ${OUTDIR} -t $TAG

##================================================================================
## 3.2 Process LAI - Aggregation to target grid
##================================================================================

FLCOVER=${OUTDIR}/esacci_lc_agg_ifs_${TAG}.nc # file generated in previous step 4.1 

## loop in months 
for MM in 01 02 03 04 05 06 07 08 09 10 11 12
do 
  FINPUT=${FCGLSINPUT}/${FINPUTstart}${MM}${FINPUTend}

  python3 -u cut_region.py -r ${GRID} -i ${FINPUT} -o ${FINPUTstart}${MM}${FINPUTend}_${REGION}

  SDATE=${YYYY}${MM}15
  FOUTPUT=${OUTDIR}/xlai_${MM} 

 # python3 -u agg_cgls.py -i $FINPUT \
 #                            -d $SDATE \
 #                            -cci $FLCOVER \
 #                            -g $XDATA_DIR/lsm \
 #                            -o $FOUTPUT

  python3 -u agg_cgls.py -i ${FINPUTstart}${MM}${FINPUTend}_${REGION} \
                             -d $SDATE \
                             -cci $FLCOVER \
                             -g $GRID \
                             -s 1 \
                             -o $FOUTPUT
done 

cdo mergetime ${OUTDIR}/xlai_* ${OUTDIR}/month_lai_v0

rm -f ${OUTDIR}/xlai_??

##================================================================================
## 3.3 Process LAI - Process LAI - Dissaggregation
##================================================================================

## note that the python script, expects to find in the -i directory (current directory):
## lsm cvl cvh tvl tvh month_lai month_lail_v0 month_laih_v0 esacci_cwt_ifs.csv
## and will create: month_lail , month_laih  

#cp $XDATA_DIR/lsm . 

cp $GRID ${OUTDIR}/lsm
ncrename -v landsea,lsm ${OUTDIR}/lsm
cp esacci_cwt_ifs.csv ${OUTDIR}

python3 -u disagg_cgls.py -i ${OUTDIR} -t nc -d ${YYYY}

cdo -O merge ${OUTDIR}/esacci_lc_agg_ifs_REG_LL.nc ${OUTDIR}/month_lai_${YYYY}.nc ${OUTDIR}/surfviav_${YYYY}

echo "DONE esacci"

#=============================================================================

echo "finished update_IAV_veg"
