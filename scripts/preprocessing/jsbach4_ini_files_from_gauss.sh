#!/bin/bash
#_____________________________________________________________________________
# Script to construct for a given ICON grid file initial fields and boundary conditions
# for the land model JSBACH4+ and the ICON atmosphere by interpolating existing files
# from ECHAM6 on Gaussian grid
#
# Reiner Schnur, MPI-M, 2017-04-10 initial version
# 
# Based on scripts from Veronika Gayler, Thomas Raddatz and Marco Giorgetta
#
set -eu
#_____________________________________________________________________________

rootdir=/pool/data/ICON/grids

basedir=${rootdir}/private/jsbach
basedir=/scratch/local1/m212005/ini_files/jsbach4/icon

workdir=${basedir}/workdir

grid_path=/pool/data/ICON/grids/private/mpim/icon_preprocessing/source/grids  # todo
grid_path=/pool/data/ICON/grids/private/r2b4_amip/r0002

bisections="0 1 2 3 4 5 6 7 8 9 10 11"
bisections="4"

remap_path=${rootdir}/private/mpim/icon_preprocessing/source/remap_weights    # todo
remap_path=/pool/data/ICON/grids/private/r2b4_amip/r0002

# Currently, only 0 or 11
npfts=11
npfts=0
[[ $npfts -gt 0 ]] && pft_tag="${npfts}pfts_" || pft_tag=""

year=1976

input_path=/pool/data/JSBACH/input/r0005/T255
[[ npfts -gt 0 ]] && \
    input_file=jsbach_T255TP6M_fractional_${npfts}tiles_5layers_${year}.nc || \
    input_file=jsbach_T255TP6M_fractional_11tiles_5layers_${year}.nc
input_path=/pool/data/JSBACH/input/r0005/T63
[[ npfts -gt 0 ]] && \
    input_file=jsbach_T63TP04_fractional_${npfts}tiles_5layers_${year}.nc || \
    input_file=jsbach_T63TP04_fractional_11tiles_5layers_${year}.nc
oro_input_file=/pool/data/ECHAM6/input/r0005/T63/T63TP04_jan_surf.nc

# Note: glac must be at beginning of varlist
varlist="glac lake init_moist snow roughness_length albedo elevation orography_std_dev 
fao forest_fract maxmoist lai_clim veg_fract surf_temp veg_ratio_max 
albedo_veg_vis albedo_veg_nir albedo_soil_vis albedo_soil_nir roughness_length_oro
soil_depth soil_porosity pore_size_index soil_field_cap heat_capacity heat_conductivity 
moisture_pot hyd_cond_sat wilting_point bclapp root_depth layer_moist
OROMEA OROSTD OROSIG OROGAM OROTHE OROPIC OROVAL"
[[ ${npfts} > 0 ]] && varlist="${varlist} cover_fract"

coupled=false    # land sea mask for coupled experiment?
hd_file=/pool/data/ICON/coupled_input_temp/aloy_pre04/Land/hd_para_icon_r2b4_pre04_vs1_5.nc

output_path=${rootdir}/public/mpim
output_path=/scratch/local1/m212005/data/icon

debug=false             # true/false
debug_grid=t63grid
debug_remap_scheme=nn

dryrun=false  # Process files in workdir, but don't copy them to the destination folder (true/false)
clean=true    # Delete workdir after finished (true/false)

#_____________________________________________________________________________

set +u
. ${MODULESHOME}/init/bash
module unload cdo nco netcdf
module load cdo/1.8.0-gccsys
module load nco/4.5.2
module load netcdf
set -u

#_____________________________________________________________________________
#
# Function definitions
#_____________________________________________________________________________
#
function generate_land_sea_mask {

    echo "Generating land sea mask ..."

    if [[ ${coupled} == "true" ]]
    then
        if [[ ${bisection} -ne 4 ]]
        then
            echo "Generation of files for coupled experiments currently only supported for r2b4"
            exit
        else
            ${cdo} chname,lsf,notsea /pool/data/ICON/coupled_input_temp/aloy_pre04/lsf.remapycon.pre04_0158km_2016.nc \
                    ${grid_name}_notsea.nc
        fi
    else
        ${cdo} setname,notsea -gtc,0 -selvar,cell_sea_land_mask ${icon_grid_file} ${grid_name}_notsea.nc
    fi

    #[[ ${debug} == "true" ]] && ${cdo} remapnn,t255grid ${grid_name}_notsea.nc debug_${grid_name}_notsea_t255.nc
    debug_file ${grid_name}_notsea.nc
}

function extrapolate_var {

    local var=$1     # Variable name

    local cdo_files="gauss_${var}.nc gauss_${var}_ext.nc"

    case ${var} in

        init_moist | layer_moist | roughness_length | root_depth | fao )
            echo "Extrapolating $var ..."
            ${cdo} -fillmiss -setmissval,0 ${cdo_files}
            ;;
        albedo )
            echo "Extrapolating $var ..."
            ${cdo} -fillmiss -setmissval,0.07 ${cdo_files}
            ;;
        albedo_veg_vis | albedo_veg_nir | albedo_soil_vis | albedo_soil_nir | surf_temp )
            echo "Extrapolating $var ..."
            ${cdo} copy ${cdo_files}
            ;;
        maxmoist )
            echo "Extrapolating $var ..."
            ${cdo} -fillmiss -setmissval,1.e-13 ${cdo_files}
            ;;
        forest_fract| snow | glac | lake )
            echo "Extrapolating $var ..."
            ${cdo} setvar,${var} -fillmiss -ifthen gauss_slf.nc ${cdo_files}
            ;;
        lai_clim | veg_fract | veg_ratio_max | roughness_length_oro | \
        soil_depth | bclapp | soil_field_cap | heat_capacity | heat_conductivity | hyd_cond_sat | \
        moisture_pot | pore_size_index | soil_porosity | wilting_point )
            echo "Extrapolating $var ..."
            ${cdo} setvar,${var} -fillmiss -sub -setmissval,-9e33 -setmisstoc,-9.e33 -setmissval,0 \
                    -add gauss_${var}.nc gauss_slf.nc gauss_slf.nc gauss_${var}_ext.nc
            ;;
        cover_fract )
            echo "Extrapolating $var ..."
            temp=$(tempfile -d .)
            ${cdo} splitlevel gauss_${var}.nc cfract
            ${cdo} sub cfract000001.nc gauss_glac.nc cfract_tmp.nc && mv cfract_tmp.nc cfract000001.nc
            ${cdo} merge cfract0000??.nc ${temp} && rm cfract0000??.nc
            ${cdo} setvar,${var} -fillmiss -ifnotthen gauss_glac.nc -ifthen gauss_slf.nc ${temp} gauss_${var}_ext.nc
            rm ${temp}
            ;;
        #* ) echo "Not extrapolating $var ..."
        #   ;;
    esac

    return 0

}

function remap_var {

    local var=$1

    case ${var} in

        fao )
            echo "Remapping $var ..."
            ${cdo} setmissval,-9e33 -remapnn,${icon_grid_file} gauss_${var}_ext.nc ${grid_name}_${var}_remap.nc
            ;;
        albedo )
            echo "Remapping $var ..."
            temp=$(tempfile -d .)
            ${cdo} remap,${icon_grid_file},${icon_remap_weights} gauss_${var}_ext.nc ${temp}
            ${cdo} setvar,${var} -setmisstoc,0.07 -ifthenelse ${grid_name}_glac_remap.nc \
                -mulc,0.7 ${grid_name}_glac_remap.nc ${temp} ${grid_name}_${var}_remap.nc
            rm ${temp}
            ;;
        albedo_veg_vis | albedo_veg_nir | albedo_soil_vis | albedo_soil_nir )
            echo "Remapping $var ..."
            temp=$(tempfile -d .)
            ${cdo} remap,${icon_grid_file},${icon_remap_weights} gauss_${var}_ext.nc ${temp}
            ${cdo} setvar,${var} -setmisstoc,0.07 -ifthen ${grid_name}_notsea.nc ${temp} ${grid_name}_${var}_remap.nc
            rm ${temp}
            ;;
        elevation )
            echo "Extracting elevation from ICON grid file"
            ${cdo} setvar,elevation -selvar,cell_elevation ${icon_grid_file} ${grid_name}_${var}_remap.nc
            ;;
        ORO* )
            echo "Remapping $var ..."
            temp=$(tempfile -d .)
            ${cdo} remap,${icon_grid_file},${icon_remap_weights} gauss_${var}.nc ${temp}
            ${cdo} chname,${var},$(echo $var | tr 'A-Z', 'a-z') ${temp} ${grid_name}_${var}_remap.nc
            rm ${temp}
            ;;
        cover_fract )
            # npft should be consistent with cover_fract from the JSBACH3 input file
            echo "Remapping $var ..."
            temp=$(tempfile -d .)
            ${cdo} remap,${icon_grid_file},${icon_remap_weights} gauss_${var}_ext.nc ${temp}
            ${cdo} setvar,${var} -setmissval,-9.e33 -ifthen ${grid_name}_notsea.nc ${temp} ${grid_name}_${var}_remap.nc

            ${cdo} splitlevel ${grid_name}_${var}_remap.nc cfract
            local n=1
            while [[ ${n} -le ${npfts} ]]
            do
                nt=$(printf "%2.2i" $n)
                ${cdo} mul cfract0000${nt}.nc -mulc,-1 -subc,1 ${grid_name}_glac_remap.nc cfract${nt}.nc && rm cfract0000${nt}.nc
                (( n = n + 1 ))
            done
            ${cdo} add cfract01.nc ${grid_name}_glac_remap.nc ${temp} && mv ${temp} cfract01.nc
            rm ${grid_name}_${var}_remap.nc
            ${cdo} merge cfract??.nc ${grid_name}_${var}_remap.nc
            ;;
        roughness_length | roughness_length_oro )
            echo "Remapping $var ..."
            temp=$(tempfile -d .)
            ${cdo} remap,${icon_grid_file},${icon_remap_weights} gauss_${var}_ext.nc ${temp}
            ${cdo} setvar,${var} -setmisstoc,0.005 -ifthen ${grid_name}_notsea.nc ${temp} ${grid_name}_${var}_remap.nc
            ;;
        glac )
            echo "Remapping $var ..."
            temp=$(tempfile -d .)
            ${cdo} remap,${icon_grid_file},${icon_remap_weights} gauss_${var}_ext.nc ${temp}
            ${cdo} setvar,${var} -gec,0.5 -setmissval,-9.e33 -ifthen ${grid_name}_notsea.nc ${temp} ${grid_name}_${var}_remap.nc
            ;;
        surf_temp )
            echo "Remapping $var ..."
            temp=$(tempfile -d .)
            ${cdo} remap,${icon_grid_file},${icon_remap_weights} gauss_${var}_ext.nc ${temp}
            ${cdo} setvar,${var} -setmissval,-9.e33 -ifthen ${grid_name}_notsea.nc ${temp} ${grid_name}_${var}_remap.nc
            ;;
        * )
            if [[ -f gauss_${var}_ext.nc ]]
            then
                echo "Remapping $var ..."
                temp=$(tempfile -d .)
                ${cdo} remap,${icon_grid_file},${icon_remap_weights} gauss_${var}_ext.nc ${temp}
                ${cdo} setvar,${var} -setmissval,-9.e33 -ifthen ${grid_name}_notsea.nc ${temp} ${grid_name}_${var}_remap.nc
                ${cdo} mul ${grid_name}_${var}_remap.nc -mulc,-1 -subc,1 ${grid_name}_glac_remap.nc ${temp}
                mv ${temp} ${grid_name}_${var}_remap.nc
                #rm ${temp}
            fi
            ;;
    esac

    # if [[ ${var} == "glac" ]]
    # then
    #     temp=$(tempfile -d .)
    #     ${cdo} gec,0.5 ${grid_name}_glac_remap.nc ${temp}
    #     mv ${temp} ${grid_name}_glac_remap.nc
    # fi

    [[ -f ${grid_name}_${var}_remap.nc ]] && debug_file ${grid_name}_${var}_remap.nc

    return 0
}

function crosscheck_vars {

    temp=$(tempfile -d .)
    ${cdo} ifthenelse ${grid_name}_landnoglacier.nc -setrtoc,0,0.1,0.1 ${grid_name}_soil_depth_remap.nc ${grid_name}_soil_depth_remap.nc ${temp}
    mv ${temp} ${grid_name}_soil_depth_remap.nc

    temp=$(tempfile -d .)
    ${cdo} mul ${grid_name}_soil_field_cap_remap.nc ${grid_name}_soil_depth_remap.nc max_soilwater.nc
    ${cdo} min ${grid_name}_maxmoist_remap.nc max_soilwater.nc ${temp}
    mv ${temp} ${grid_name}_maxmoist_remap.nc
    rm max_soilwater.nc

    temp=$(tempfile -d .)
    cdo setmisstoc,280 ${grid_name}_surf_temp_remap.nc ${temp}
    mv ${temp} ${grid_name}_surf_temp_remap.nc
    debug_file ${grid_name}_surf_temp_remap.nc
}

function generate_fractions {

    ${cdo} setvar,fract_glac -setmisstoc,0 -gec,0.5 ${grid_name}_glac_remap.nc ${grid_name}_fract_glac.nc
    debug_file ${grid_name}_fract_glac.nc

    ${cdo} setvar,landnoglacier -mul ${grid_name}_notsea.nc -ltc,0.5 ${grid_name}_fract_glac.nc ${grid_name}_landnoglacier.nc
    debug_file ${grid_name}_landnoglacier.nc

    ${cdo} setvar,sea -mulc,-1 -subc,1 ${grid_name}_notsea.nc ${grid_name}_sea.nc
    debug_file ${grid_name}_sea.nc

    ${cdo} setvar,fract_lake -setmisstoc,0 ${grid_name}_lake_remap.nc ${grid_name}_fract_lake.nc
    debug_file ${grid_name}_fract_lake.nc

    ${cdo} setvar,fract_land -mulc,-1 -subc,1 ${grid_name}_fract_lake.nc ${grid_name}_fract_land.nc
    debug_file ${grid_name}_fract_land.nc

    [[ ${npfts} -eq 0 ]] && return 0

    # Generate PFT fractions
    local i=1
    while [[ ${i} -le ${npfts} ]]
    do
        ipft=$(printf "%2.2i" $i)
        ${cdo} sellevel,$i ${grid_name}_cover_fract_remap.nc pft${ipft}.nc
        if [[ ${i} -eq 1 ]]
        then
            ${cdo} setrtoc,1,1.1,0 pft${ipft}.nc zpft && mv zpft pft${ipft}.nc
        fi
        ncrename -v cover_fract,fract_pft${ipft} pft${ipft}.nc
        ncwa -a ntiles pft${ipft}.nc zpft 
        ${cdo} setmisstoc,0 zpft pft${ipft}.nc && rm zpft
        (( i = i + 1 ))
    done
    ${cdo} enssum pft*.nc fract_sum.nc
    ${cdo} mulc,-1 -subc,1 -gtc,0 -setmisstoc,0 fract_sum.nc fract_sum2.nc && mv fract_sum2.nc fract_sum.nc
    ${cdo} add pft03.nc fract_sum.nc zpft && mv zpft pft03.nc
    for i in $(seq 1 $npfts)
    do
        ipft=$(printf "%2.2i" $i)
        ${cdo} setname,fract_pft${ipft} pft${ipft}.nc ${grid_name}_fract_pft${ipft}_remap.nc && rm pft${ipft}.nc
        debug_file ${grid_name}_fract_pft${ipft}_remap.nc
    done        

}

function copy_uuid {

    local file=$1

    uuidOfHGrid=$(ncdump -h ${icon_grid_file} | grep ':uuidOfHGrid = ' | cut -f2 -d'"')
    ncatted -a uuidOfHGrid,global,o,c,${uuidOfHGrid} ${file}

}

function generate_output_files {

    local outfile=""

    outfile="bc_land_frac_${pft_tag}${year}.nc"
    echo "Generating output file ${outfile} ..."
    if [[ ${npfts} -gt 0 ]]
    then
        ${cdo} merge \
            ${grid_name}_notsea.nc \
            ${grid_name}_fract_glac.nc \
            ${grid_name}_sea.nc \
            ${grid_name}_fract_lake.nc \
            ${grid_name}_fract_land.nc \
            ${grid_name}_landnoglacier.nc \
            ${grid_name}_veg_ratio_max_remap.nc \
            ${grid_name}_fract_pft*_remap.nc \
            ${outfile}
    else
        ${cdo} merge \
            ${grid_name}_notsea.nc \
            ${grid_name}_fract_glac.nc \
            ${grid_name}_sea.nc \
            ${grid_name}_fract_lake.nc \
            ${grid_name}_fract_land.nc \
            ${grid_name}_landnoglacier.nc \
            ${grid_name}_veg_ratio_max_remap.nc \
            ${outfile}
    fi
    ncrename -d cell,ncells ${outfile} >& /dev/null
    copy_uuid ${outfile}
    debug_file ${outfile}

    outfile="bc_land_phys_${pft_tag}${year}.nc"
    echo "Generating output file ${outfile} ..."
    ${cdo} merge \
        ${grid_name}_lai_clim_remap.nc \
        ${grid_name}_veg_fract_remap.nc \
        ${grid_name}_roughness_length_remap.nc \
        ${grid_name}_roughness_length_oro_remap.nc \
        ${grid_name}_albedo_remap.nc \
        ${grid_name}_albedo_veg_vis_remap.nc \
        ${grid_name}_albedo_veg_nir_remap.nc \
        ${grid_name}_albedo_soil_vis_remap.nc \
        ${grid_name}_albedo_soil_nir_remap.nc \
        ${grid_name}_forest_fract_remap.nc \
        ${outfile}
    ncrename -d cell,ncells ${outfile} >& /dev/null
    copy_uuid ${outfile}
    debug_file ${outfile}

    outfile="bc_land_soil_${pft_tag}${year}.nc"
    echo "Generating output file ${outfile} ..."
    ${cdo} merge \
        ${grid_name}_soil_depth_remap.nc \
        ${grid_name}_root_depth_remap.nc \
        ${grid_name}_bclapp_remap.nc \
        ${grid_name}_soil_field_cap_remap.nc \
        ${grid_name}_heat_capacity_remap.nc \
        ${grid_name}_heat_conductivity_remap.nc \
        ${grid_name}_hyd_cond_sat_remap.nc \
        ${grid_name}_moisture_pot_remap.nc \
        ${grid_name}_pore_size_index_remap.nc \
        ${grid_name}_soil_porosity_remap.nc \
        ${grid_name}_wilting_point_remap.nc \
        ${grid_name}_fao_remap.nc \
        ${grid_name}_maxmoist_remap.nc \
        ${outfile}
    ncrename -d cell,ncells ${outfile} >& /dev/null
    copy_uuid ${outfile}
    debug_file ${outfile}

    outfile="bc_land_sso_${year}.nc"
    echo "Generating output file ${outfile} ..."
    ${cdo} merge \
        ${grid_name}_elevation_remap.nc \
        ${grid_name}_OROMEA_remap.nc \
        ${grid_name}_OROSTD_remap.nc \
        ${grid_name}_OROSIG_remap.nc \
        ${grid_name}_OROGAM_remap.nc \
        ${grid_name}_OROTHE_remap.nc \
        ${grid_name}_OROPIC_remap.nc \
        ${grid_name}_OROVAL_remap.nc \
        ${outfile}
    ncrename -d cell,ncells ${outfile} >& /dev/null
    copy_uuid ${outfile}
    debug_file ${outfile}

    outfile="ic_land_soil_${pft_tag}${year}.nc"
    echo "Generating output file ${outfile} ..."
    ${cdo} merge \
        ${grid_name}_init_moist_remap.nc \
        ${grid_name}_surf_temp_remap.nc \
        ${grid_name}_layer_moist_remap.nc \
        ${grid_name}_snow_remap.nc \
        ${outfile}
    ncrename -d cell,ncells ${outfile} >& /dev/null
    copy_uuid ${outfile}
    debug_file ${outfile}

}

function debug_file {

    local file=$1

    if [[ ${debug} == "true" ]]
    then
        if [[ ! -f debug_weights.nc ]]
        then
            #${cdo} -f nc gen${debug_remap_scheme},${debug_grid} -random,debug_griddes.txt debug_weights.nc
            ${cdo} gen${debug_remap_scheme},${debug_grid} ${grid_name}_notsea.nc debug_weights.nc
        fi
        ${cdo} remap,${debug_grid},debug_weights.nc ${file} debug_${file}
    fi

}
#_____________________________________________________________________________
#
# Main script
#_____________________________________________________________________________
#
mkdir -p ${workdir} || true
cd ${workdir}
rm -rf *

cdo -b F64 splitname ${input_path}/${input_file} gauss_
cdo -b F64 splitname ${oro_input_file} gauss_

for bisection in $bisections
do
  np=  
  n0=$(printf "%2.2i" $bisection)
  (( n1 = bisection + 1 ))
  n1=$(printf "%2.2i" $n1)
  if [[ $bisection -ge 8  ]]
  then
      np="-P 8"
  fi
  cdo="cdo -s -b F64 $np"
  #cdo="cdo -b F64 $np"

  grid_name=R02B${n0}

  icon_grid=icon_grid_00${n1}_R02B${n0}_G.nc           # todo
  remap_weights=wgtdis_T255_to_R02B${n0}.nc            # todo
  icon_grid=r2b4_amip.nc
  #remap_source=r2b4_amip_cell.nc
  remap_weights=remapdis_t63grid_r2b4_amip_cell.nc

  icon_grid_file=${grid_path}/${icon_grid}
  icon_remap_weights=${remap_path}/${remap_weights}

  echo Processing $icon_grid ...

  generate_land_sea_mask

  for var in $varlist
  do
    extrapolate_var $var 
  done

  for var in $varlist
  do
    remap_var $var
  done

  generate_fractions

  crosscheck_vars

  generate_output_files

  output_dir="${output_path}/R2B${bisection}/land"
  mkdir -p ${output_dir} || true

  if [[ ${dryrun} != "true" ]]
  then
    echo "Copying output files to ${output_dir}"
    cp -p bc_land_frac_${pft_tag}$year.nc ${output_dir}
    cp -p bc_land_phys_${pft_tag}$year.nc ${output_dir}
    cp -p bc_land_soil_${pft_tag}$year.nc ${output_dir}
    cp -p bc_land_sso_$year.nc            ${output_dir}
    cp -p ic_land_soil_${pft_tag}$year.nc ${output_dir}

    [[ ${coupled} == "true" ]] && cp -p ${hd_file} ${output_dir}/bc_land_hd.nc
  fi

  if [[ ${clean} == "true" ]]
  then
    echo "Cleaning up ..."
    cd ${workdir}/..
    rm -rf ${workdir}
  fi

done

