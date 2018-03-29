#!/bin/bash
set +ex

#module load  nco

##############################################################################################
### Dieses Programm holt eine Liste von Variablen z.B. c-Pools aus einem ICON Restartfile
### und packt sie in ein eigenes File zusammen ###
###
### Note: Tested with nco version 4.6.7 ... older versions might not work.
###       As of March 2018, it doesn't work on MPI workstations because nco is too old.
###       It does work on mistral.
##############################################################################################

### User Interface
INFILE="JS045_restart_atm_19870201T000000Z.nc"
OUTFILE="ic_land_carbon.nc"

PROCESS="carbon"
VARNAME_LIST="c_green   c_woods   c_reserve 
              c_acid_ag1   c_water_ag1  c_ethanol_ag1  c_nonsoluble_ag1
              c_acid_bg1   c_water_bg1  c_ethanol_bg1  c_nonsoluble_bg1  c_humus_1
              c_acid_ag2   c_water_ag2  c_ethanol_ag2  c_nonsoluble_ag2
              c_acid_bg2   c_water_bg2  c_ethanol_bg2  c_nonsoluble_bg2  c_humus_2"
TILE_LIST="box  land  veg  pft01  pft02  pft03  pft04  pft05  pft06  pft07  pft08  pft09  pft10  pft11 "

### Evtl. löschen des vorhandenen Output-Files
if [ -f  ${OUTFILE} ];then
   echo "Das Outputfile ${OUTFILE} ist bereits vorhanden. Soll es überschrieben werden?"
   echo "RETURN =File wird überschrieben."
   echo "XXX    =Beliebige andere Eingabe fuert zum Abbruch."
   read y
   if [ ! ${y} ] ; then 
      rm ${OUTFILE}
      echo "Das bereits vorhandene File ${OUTFILE} wurde gelöscht."
      echo ""
   else
      exit
   fi
fi

### Variablen aus Restartfile holen
for VARNAME in ${VARNAME_LIST} ; do
   echo -e "\033[01;31;40m"
   echo "=== ${VARNAME} ==="
   tput sgr0

   for TILE in ${TILE_LIST} ; do
      STREAM_NAME="${PROCESS}_${VARNAME}_${TILE}"
      echo "${STREAM_NAME}"

      cdo selvar,${STREAM_NAME}                      ${INFILE}    ZWERG1
      #cdo chname,${STREAM_NAME},${VARNAME}_${TILE}   ZWERG1       ZWERG2  # workaround for long varnames
      ncks -A                                        ZWERG1       ZWERG2
   done
done

ncks -C -x -v time            ZWERG2       ZWERG3      # Delete variable time
ncwa -O -a time               ZWERG3       ${OUTFILE}  # Delete dimension time
ncrename -d ncells_2,ncells ${OUTFILE}   ! This needs a recent nco version!

### Clean up
rm ZWERG1  ZWERG2  ZWERG3

### Ende
echo -e "\033[00;31;48m" 
echo ""
echo "Normales Shellende erreicht !"
echo ""
tput sgr0
exit

