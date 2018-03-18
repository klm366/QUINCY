#!/bin/bash
set +ex

#module load  nco

########################################################################################################################################
### Dieses Programm holt eine Liste von Variablen z.B. c-Pools aus einem ICON Restartfile und packt sie in ein eigenes File zusammen ###
########################################################################################################################################

### User Interface
INFILE="JS045_restart_atm_19870201T000000Z.nc"
OUTFILE="ic_land_Cpools.nc"

PROCESS="carbon"
VARNAME_LIST="c_green_ta   c_woods_ta   c_reserve_ta 
              c_acid_ag1_ta   c_water_ag1_ta  c_ethanol_ag1_ta  c_nonsoluble_ag1_ta
              c_acid_bg1_ta   c_water_bg1_ta  c_ethanol_bg1_ta  c_nonsoluble_bg1_ta  c_humus_1_ta
              c_acid_ag2_ta   c_water_ag2_ta  c_ethanol_ag2_ta  c_nonsoluble_ag2_ta
              c_acid_bg2_ta   c_water_bg2_ta  c_ethanol_bg2_ta  c_nonsoluble_bg2_ta  c_humus_2_ta"
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

### Clean up
rm ZWERG1  ZWERG2  ZWERG3

### Ende
echo -e "\033[00;31;48m" 
echo ""
echo "Normales Shellende erreicht !"
echo ""
tput sgr0
exit

