#!/bin/ksh
#------------------------------------------------------------------------------
# generate jsbach initial boundary conditions and initial conditions for icon grid
# Veronika Gayler, MPI-M, Hamburg, 2013
#------------------------------------------------------------------------------
set -e

#cdo="cdo -f nc2"
cdo="cdo -f nc"

year=${year_cf}
[[ ${echam_fractional} = true ]] && fractional=fractional_ || fractional=""

# input files
jsbachfile=${srcdir}/input/jsbach/${grid}/jsbach_${res_atm}${res_oce}_${fractional}${ntiles}tiles_5layers_${year}.nc
ifile_topo=${srcdir}/input/echam6/${res_atm}${res_oce}_jan_surf.nc
filename=${grid}_jsbach_${ntiles}tiles_${year}

#${cdo} copy ${pool_soil}/${soilparamfile} ${soilparamfile}
[[ -f ${soilparamfile} ]] || ln -s ${pool_soil}/${soilparamfile} .

# preparations
#--------------

${cdo} copy ${jsbachfile} ${jsbachfile}_tmp && mv ${jsbachfile}_tmp ${jsbachfile}
${cdo} copy ${ifile_topo} ${ifile_topo}_tmp && mv ${ifile_topo}_tmp ${ifile_topo}

# generate remapping matrices if not yet available
[[ -f rmp_${res_atm}_to_${grid} ]] || ${cdo} gendis,${gridfile} ${jsbachfile} rmp_${res_atm}_to_${grid}

# land sea mask on Gaussian grid
${cdo} selvar,slf ${jsbachfile} slm_${res_atm}.nc

# glacier mask on Gaussian grido
${cdo} selvar,glac ${jsbachfile} glac_${res_atm}.nc

# land sea mask on icon grid
${cdo} chname,cell_sea_land_mask,slf -gtc,0 -selvar,cell_sea_land_mask ${gridfile} slf_${grid}.nc
#${cdo} chname,cell_sea_land_mask,slf -selvar,cell_sea_land_mask ${srcdir}/lsf.remapycon.tri0039km.editslo-amip_r2b4.nc slf_${grid}.nc
#${cdo} chname,lsf,slf -selvar,lsf ${srcdir}/lsf.remapycon.pre04_0158km_2016.nc slf_${grid}.nc

rm -f ${filename}_*.tmp
varlist="init_moist snow slf slm roughness_length albedo elevation glac fao \
         maxmoist forest_fract lai_clim veg_fract surf_temp cover_fract cover_type veg_ratio_max \
         albedo_veg_vis albedo_veg_nir albedo_soil_vis albedo_soil_nir roughness_length_oro \
         natural_veg bclapp soil_field_cap heat_capacity heat_conductivity hyd_cond_sat moisture_pot pore_size_index soil_porosity \
         soil_depth wilting_point layer_moist root_depth"

for var in ${varlist}; do

    # extrapolate land values to the ocean on the source grid
    case ${var} in

        init_moist | layer_moist | roughness_length | root_depth | fao )
            ${cdo} -fillmiss -setmissval,0 -selvar,${var} ${jsbachfile} gauss_${var}.tmp
            rmp=rmp_${res_atm}_to_${grid} ;;
        albedo )
            ${cdo} -fillmiss -setmissval,0.07 -selvar,${var} ${jsbachfile} gauss_${var}.tmp
            rmp=rmp_${res_atm}_to_${grid} ;;
        maxmoist )
            ${cdo} -fillmiss -setmissval,1.e-13 -selvar,${var} ${jsbachfile} gauss_${var}.tmp
            rmp=rmp_${res_atm}_to_${grid} ;;
        cover_fract | natural_veg )
            ${cdo} splitlevel -selvar,${var} ${jsbachfile} ${var}
            ${cdo} sub ${var}000001.nc glac_${res_atm}.nc ${var}000001.tmp
            mv ${var}000001.tmp ${var}000001.nc
            rm -f gauss_${var}.tmp1
            ${cdo} merge ${var}0000??.nc gauss_${var}.tmp1
            rm ${var}0000??.nc
            ${cdo} setvar,${var} -fillmiss -ifnotthen glac_${res_atm}.nc -ifthen slm_${res_atm}.nc \
                gauss_${var}.tmp1 gauss_${var}.tmp
            rm gauss_${var}.tmp1
            rmp=rmp_${res_atm}_to_${grid} ;;
        cover_type )
            ${cdo} splitlevel -selvar,${var} ${jsbachfile} ${var}
            ${cdo} add ${var}000001.nc glac_${res_atm}.nc ${var}000001.tmp
            mv ${var}000001.tmp ${var}000001.nc
            rm -f gauss_${var}.tmp1
            ${cdo} merge ${var}0000??.nc gauss_${var}.tmp1
            rm ${var}0000??.nc
            ${cdo} setvar,${var} -fillmiss -ifthen slm_${res_atm}.nc gauss_${var}.tmp1 gauss_${var}.tmp
            rm gauss_${var}.tmp1
            rmp=rmp_${res_atm}_to_${grid} ;;
        forest_fract | glac | snow )
            ${cdo} setvar,${var} -fillmiss -ifthen slm_${res_atm}.nc -selvar,${var} ${jsbachfile} \
                gauss_${var}.tmp
            rmp=rmp_${res_atm}_to_${grid} ;;
        lai_clim | veg_fract | veg_ratio_max | roughness_length_oro )
            ${cdo} setvar,${var} -fillmiss -sub -setmissval,-9e33 -setmisstoc,-9e33 -setmissval,0 \
                         -add slm_${res_atm}.nc -selvar,${var} ${jsbachfile} \
                         slm_${res_atm}.nc gauss_${var}.tmp
            rmp=rmp_${res_atm}_to_${grid} ;;
        albedo_veg_vis | albedo_veg_nir | albedo_soil_vis | albedo_soil_nir | surf_temp )
            ${cdo} selvar,${var} ${jsbachfile} gauss_${var}.tmp
            rmp=rmp_${res_atm}_to_${grid} ;;
    esac

    # do the remapping 
    # Soil parameter fields are just read in (and renamed) from already remapped soil parameter file
    case ${var} in

 	elevation )
            #${cdo} setvar,elevation -ifthen slf_${grid}.nc -selvar,cell_elevation ${gridfile} \
            ${cdo} setvar,elevation -selvar,cell_elevation ${gridfile} \
                       ${filename}_elevation.tmp ;;
	slf ) 
            #rm -f ${filename}_${var}.tmp1 ${filename}_${var}.tmp ;;
            cp -p slf_${grid}.nc ${filename}_slf.tmp ;;
	slm ) 
            ${cdo} setvar,slm -gtc,0 slf_${grid}.nc ${filename}_slm.tmp ;;
        fao )
            ${cdo} remapnn,${gridfile} gauss_${var}.tmp ${filename}_${var}.tmp1
            #${cdo} setvar,${var} -setmissval,-9e+33 -ifthen slf_${grid}.nc ${filename}_${var}.tmp1 \
            ${cdo} setvar,${var} -setmissval,-9e+33 ${filename}_${var}.tmp1 \
                ${filename}_${var}.tmp ;;
        maxmoist )
            [[ ${dynveg} = true ]] && \
              #${cdo} setvar,${var} -setmissval,-9e+33 -ifthen slf_${grid}.nc -selvar,wcap_dynveg \
              ${cdo} setvar,${var} -setmissval,-9e+33 -selvar,wcap_dynveg \
                ${soilparamfile} ${filename}_${var}.tmp || \
              #${cdo} setvar,${var} -setmissval,-9e+33 -ifthen slf_${grid}.nc -selvar,wcap \
              ${cdo} setvar,${var} -setmissval,-9e+33 -selvar,wcap \
                ${soilparamfile} ${filename}_${var}.tmp
            ;;
	soil_depth )
            #${cdo} setvar,${var} -setmissval,-9e+33 -ifthen slf_${grid}.nc -selvar,soildepth ${soilparamfile} ${filename}_${var}.tmp ;;
            ${cdo} setvar,${var} -setmissval,-9e+33 -selvar,soildepth ${soilparamfile} ${filename}_${var}.tmp ;;
        bclapp )
            #${cdo} setvar,${var} -setmissval,-9e+33 -ifthen slf_${grid}.nc -selvar,bclapp ${soilparamfile} ${filename}_${var}.tmp ;;
            ${cdo} setvar,${var} -setmissval,-9e+33 -selvar,bclapp ${soilparamfile} ${filename}_${var}.tmp ;;
        soil_field_cap )
            #${cdo} setvar,${var} -setmissval,-9e+33 -ifthen slf_${grid}.nc -selvar,volfc ${soilparamfile} ${filename}_${var}.tmp ;;
            ${cdo} setvar,${var} -setmissval,-9e+33 -selvar,volfc ${soilparamfile} ${filename}_${var}.tmp ;;
        heat_capacity )
            #${cdo} setvar,${var} -setmissval,-9e+33 -ifthen slf_${grid}.nc -selvar,heatcapacity ${soilparamfile} ${filename}_${var}.tmp ;;
            ${cdo} setvar,${var} -setmissval,-9e+33 -selvar,heatcapacity ${soilparamfile} ${filename}_${var}.tmp ;;
        heat_conductivity )
            #${cdo} setvar,${var} -setmissval,-9e+33 -ifthen slf_${grid}.nc -selvar,heatcond ${soilparamfile} ${filename}_${var}.tmp ;;
            ${cdo} setvar,${var} -setmissval,-9e+33 -selvar,heatcond ${soilparamfile} ${filename}_${var}.tmp ;;
        hyd_cond_sat )
            #${cdo} setvar,${var} -setmissval,-9e+33 -ifthen slf_${grid}.nc -selvar,hydcond ${soilparamfile} ${filename}_${var}.tmp ;;
            ${cdo} setvar,${var} -setmissval,-9e+33 -selvar,hydcond ${soilparamfile} ${filename}_${var}.tmp ;;
        moisture_pot )
            #${cdo} setvar,${var} -setmissval,-9e+33 -ifthen slf_${grid}.nc -selvar,matrixpot ${soilparamfile} ${filename}_${var}.tmp ;;
            ${cdo} setvar,${var} -setmissval,-9e+33 -selvar,matrixpot ${soilparamfile} ${filename}_${var}.tmp ;;
        pore_size_index )
            #${cdo} setvar,${var} -setmissval,-9e+33 -ifthen slf_${grid}.nc -selvar,psi_lambda ${soilparamfile} ${filename}_${var}.tmp ;;
            ${cdo} setvar,${var} -setmissval,-9e+33 -selvar,psi_lambda ${soilparamfile} ${filename}_${var}.tmp ;;
        soil_porosity )
            #${cdo} setvar,${var} -setmissval,-9e+33 -ifthen slf_${grid}.nc -selvar,volporosity ${soilparamfile} ${filename}_${var}.tmp ;;
            ${cdo} setvar,${var} -setmissval,-9e+33 -selvar,volporosity ${soilparamfile} ${filename}_${var}.tmp ;;
        wilting_point )
            #${cdo} setvar,${var} -setmissval,-9e+33 -ifthen slf_${grid}.nc -selvar,volpwp ${soilparamfile} ${filename}_${var}.tmp ;;
            ${cdo} setvar,${var} -setmissval,-9e+33 -selvar,volpwp ${soilparamfile} ${filename}_${var}.tmp ;;
        * )
            ${cdo} remap,${gridfile},${rmp} gauss_${var}.tmp ${filename}_${var}.tmp1
 	        ${cdo} setvar,${var} -setmissval,-9e+33 -ifthen slf_${grid}.nc ${filename}_${var}.tmp1 ${filename}_${var}.tmp ;;
    esac
done

# treatment of glaciers
${cdo} gec,0.5 ${filename}_glac.tmp glac_${grid}.nc
${cdo} ltc,0.5 ${filename}_glac.tmp glac_inv_${grid}.nc
 
for var in ${varlist}; do
    case ${var} in

        cover_fract | natural_veg )
            ${cdo} splitlevel ${filename}_${var}.tmp ${filename}_${var}_
            typeset -Z2 nt=01
            while [[ ${nt} -le ${ntiles} ]]; do
                ${cdo} mul ${filename}_${var}_0000${nt}.nc -selvar,glac glac_inv_${grid}.nc \
                ${filename}_${var}_${nt}.tmp
                rm ${filename}_${var}_0000${nt}.nc
                (( nt = nt + 1 ))
            done
            ${cdo} add ${filename}_${var}_01.tmp -selvar,glac glac_${grid}.nc ${filename}_${var}_01.tmp1
            mv ${filename}_${var}_01.tmp1 ${filename}_${var}_01.tmp
            rm -f ${filename}_${var}.tmp
            ${cdo} merge ${filename}_${var}_??.tmp ${filename}_${var}.tmp
            ${cdo} setvar,${var} ${filename}_${var}.tmp ${filename}_${var}.tmp2
            rm ${filename}_${var}_??.tmp ;; 
        cover_type )
            ${cdo} splitlevel ${filename}_${var}.tmp ${filename}_${var}_
            typeset -Z2 nt=01
            ${cdo} sub ${filename}_${var}_000001.nc -selvar,glac glac_${grid}.nc ${filename}_${var}_01.tmp1
            mv ${filename}_${var}_01.tmp1 ${filename}_${var}_000001.nc
            ${cdo} nint ${filename}_${var}_000011.nc ${filename}_${var}_11.tmp1
            mv ${filename}_${var}_11.tmp1 ${filename}_${var}_000011.nc
            rm -f ${filename}_${var}.tmp
            ${cdo} merge ${filename}_${var}_0000??.nc ${filename}_${var}.tmp2
            rm ${filename}_${var}_0000??.nc ;; 
        glac )
            ${cdo} copy glac_${grid}.nc ${filename}_${var}.tmp2 ;;
        lai_clim | bclapp | soil_field_cap | forest_fract | heat_capacity | heat_conductivity | root_depth | \
        hyd_cond_sat | init_moist | layer_moist | maxmoist | moisture_pot | \
        pore_size_index | soil_porosity | snow | soil_depth | veg_fract | veg_ratio_max | wilting_point )
            ${cdo} mul ${filename}_${var}.tmp glac_inv_${grid}.nc ${filename}_${var}.tmp2 ;;
        albedo )
            ${cdo} setvar,${var} -setmisstoc,0.07 -ifthenelse glac_inv_${grid}.nc ${filename}_${var}.tmp \
            -mulc,0.7 glac_${grid}.nc ${filename}_${var}.tmp2 ;;
        albedo_veg_vis | albedo_veg_nir | albedo_soil_vis | albedo_soil_nir )
            ${cdo} setmisstoc,0.07 ${filename}_${var}.tmp ${filename}_${var}.tmp2 ;;
    esac
    [[ -f ${filename}_${var}.tmp2 ]] && mv ${filename}_${var}.tmp2  ${filename}_${var}.tmp
done

${cdo} mul slf_${grid}.nc -ltc,1 glac_${grid}.nc ${filename}_landnoglacier.tmp

${cdo} ifthenelse ${filename}_landnoglacier.tmp -setrtoc,0,0.1,0.1 ${filename}_soil_depth.tmp ${filename}_soil_depth.tmp ${filename}_soil_depth.tmp2
mv ${filename}_soil_depth.tmp2 ${filename}_soil_depth.tmp

${cdo} mul ${filename}_soil_field_cap.tmp ${filename}_soil_depth.tmp ${filename}_max_soilwater.tmp

${cdo} min ${filename}_maxmoist.tmp ${filename}_max_soilwater.tmp ${filename}_maxmoist.tmp2
mv ${filename}_maxmoist.tmp2 ${filename}_maxmoist.tmp
rm ${filename}_max_soilwater.tmp

#------------------------------------------------------------------------------
# generate boundary conditions depending on topography for icon grid
# Marco Giorgetta, MPI-M, Hamburg 2013-03-08
#------------------------------------------------------------------------------

# 1. remapping ECHAM data from the gaussian grid to the ICON grid at cell centers
# ------------------------------------------------------------------------------
# method = remapnn = nearest neighbor
#

${cdo} remapnn,${gridfile} ${ifile_topo} jan_surf_tmp.nc

# 2. land sea mask on the ICON grid
# ---------------------------------

# integer sea land mask from the icon grid file
# - ocean cell without land edge : -2
# - ocean cell with land edge    : -1
# - land cell with ocean edge    :  1
# - land cell without ocean edge :  2
#
#${cdo} selname,cell_sea_land_mask ${gridfile} integer_slm.nc
#
# floating point sea mask
# - at land : 0.
# - at sea  : 1.
#
${cdo} -b F64 chname,slf,sea -mulc,-1 -subc,1 slf_${grid}.nc ${filename}_sea.tmp
#
# floating point land mask
# - at land : 1.
# - at sea  : 0.
#
#${cdo} -b F64 chname,cell_sea_land_mask,notsea -gtc,0 integer_slm.nc ${filename}_notsea.tmp
${cdo} -b F64 chname,slf,notsea slf_${grid}.nc ${filename}_notsea.tmp

# 3. set all land parameters to 0 where at sea
# --------------------------------------------
# set surface land variables to zero at sea
# - at land : keep value
# - at sea  : set to 0
#
#${cdo} mul jan_surf_tmp.nc ${filename}_notsea.tmp jan_surf.nc
#${cdo} ifthen ${filename}_notsea.tmp jan_surf_tmp.nc jan_surf.nc
cp jan_surf_tmp.nc jan_surf.nc

# 4. separate variables transfered to new files
# ---------------------------------------------
#
# lake surface fractions
#
${cdo} -b F64 chname,ALAKE,lake -setmisstoc,0 -selname,ALAKE jan_surf.nc ${filename}_lake.tmp
#${cdo} chname,ALAKE,lake -setmisstoc,0 -setrtoc,-1,2,0 -selname,ALAKE jan_surf.nc ${filename}_lake.tmp
#
# land fraction (not sea and not lake, including glacier)
#
#${cdo} chname,notsea,land -setmisstoc,0 -sub ${filename}_notsea.tmp ${filename}_lake.tmp ${filename}_land.tmp
${cdo} -b F64 setname,land -setmisstoc,0 -sub slf_${grid}.nc ${filename}_lake.tmp ${filename}_land.tmp
#
# land fraction (not sea, not lake, not glacier)
#
${cdo} -b F64 setname,landnoglacier -setmisstoc,0 -add ${filename}_lake.tmp -sub ${filename}_land.tmp ${filename}_glac.tmp ${filename}_landnoglacier.tmp
#
# surface physical properties
#
${cdo} chname,AZ0,z0  -setmisstoc,0.005 -selname,AZ0 jan_surf.nc ${filename}_z0.tmp
#
# orography
#
${cdo} chname,OROSTD,orostd -selname,OROSTD jan_surf.nc ${filename}_orostd.tmp
${cdo} chname,OROSIG,orosig -selname,OROSIG jan_surf.nc ${filename}_orosig.tmp
${cdo} chname,OROGAM,orogam -selname,OROGAM jan_surf.nc ${filename}_orogam.tmp
${cdo} chname,OROTHE,orothe -selname,OROTHE jan_surf.nc ${filename}_orothe.tmp
${cdo} chname,OROPIC,oropic -selname,OROPIC jan_surf.nc ${filename}_oropic.tmp
${cdo} chname,OROVAL,oroval -selname,OROVAL jan_surf.nc ${filename}_oroval.tmp
${cdo} chname,OROMEA,oromea -selname,OROMEA jan_surf.nc ${filename}_oromea.tmp

# generate four separate files for boundary conditions (bc) and one file for initial conditions (ic)

rm -f ${grid}_jsbach_clim.nc ${grid}_jsbach_${ntiles}tiles_${year}.nc ${grid}_jsbach.nc

rm -f bc_land_frac_${grid}_${year}.nc bc_land_phys_${grid}_${year}.nc bc_land_soil_${grid}_${year}.nc \
    bc_land_sso_${grid}_${year}.nc ic_land_soil_${grid}_${year}.nc
# define surfaces
${cdo} merge ${filename}_slf.tmp ${filename}_glac.tmp ${filename}_sea.tmp ${filename}_notsea.tmp ${filename}_lake.tmp \
    ${filename}_land.tmp ${filename}_cover_fract.tmp ${filename}_cover_type.tmp ${filename}_natural_veg.tmp \
    ${filename}_veg_ratio_max.tmp ${filename}_landnoglacier.tmp bc_land_frac_${grid}_${year}.nc
ncrename -d cell,ncells bc_land_frac_${grid}_${year}.nc
rm ${filename}_slf.tmp ${filename}_glac.tmp ${filename}_sea.tmp ${filename}_notsea.tmp ${filename}_lake.tmp \
    ${filename}_land.tmp ${filename}_cover_fract.tmp ${filename}_cover_type.tmp ${filename}_natural_veg.tmp \
    ${filename}_veg_ratio_max.tmp ${filename}_landnoglacier.tmp
 physical properties
${cdo} merge ${filename}_lai_clim.tmp ${filename}_veg_fract.tmp ${filename}_roughness_length.tmp \
    ${filename}_roughness_length_oro.tmp ${filename}_z0.tmp ${filename}_albedo.tmp ${filename}_albedo_veg_vis.tmp \
    ${filename}_albedo_veg_nir.tmp ${filename}_albedo_soil_vis.tmp ${filename}_albedo_soil_nir.tmp \
    ${filename}_forest_fract.tmp bc_land_phys_${grid}_${year}.nc
ncrename -d cell,ncells bc_land_phys_${grid}_${year}.nc
rm ${filename}_lai_clim.tmp ${filename}_veg_fract.tmp ${filename}_roughness_length.tmp ${filename}_roughness_length_oro.tmp \
    ${filename}_z0.tmp ${filename}_albedo.tmp ${filename}_albedo_veg_vis.tmp ${filename}_albedo_veg_nir.tmp \
    ${filename}_albedo_soil_vis.tmp ${filename}_albedo_soil_nir.tmp ${filename}_forest_fract.tmp
# soil properties
${cdo} merge ${filename}_fao.tmp ${filename}_maxmoist.tmp ${filename}_bclapp.tmp ${filename}_soil_field_cap.tmp \
    ${filename}_heat_capacity.tmp ${filename}_heat_conductivity.tmp ${filename}_hyd_cond_sat.tmp ${filename}_moisture_pot.tmp \
    ${filename}_pore_size_index.tmp ${filename}_soil_porosity.tmp ${filename}_soil_depth.tmp ${filename}_wilting_point.tmp \
    ${filename}_root_depth.tmp bc_land_soil_${grid}_${year}.nc
${cdo} -O setgrid,${filename}_fao.tmp bc_land_soil_${grid}_${year}.nc z.nc && mv z.nc bc_land_soil_${grid}_${year}.nc
ncrename -d cell,ncells bc_land_soil_${grid}_${year}.nc
rm ${filename}_fao.tmp ${filename}_maxmoist.tmp ${filename}_bclapp.tmp ${filename}_soil_field_cap.tmp \
    ${filename}_heat_capacity.tmp ${filename}_heat_conductivity.tmp ${filename}_hyd_cond_sat.tmp ${filename}_moisture_pot.tmp \
    ${filename}_pore_size_index.tmp ${filename}_soil_porosity.tmp ${filename}_soil_depth.tmp ${filename}_wilting_point.tmp \
    ${filename}_root_depth.tmp
# topography
${cdo} merge ${filename}_elevation.tmp ${filename}_orostd.tmp ${filename}_orosig.tmp ${filename}_orogam.tmp ${filename}_orothe.tmp \
    ${filename}_oropic.tmp ${filename}_oroval.tmp ${filename}_oromea.tmp bc_land_sso_${grid}_${year}.nc
ncrename -d cell,ncells bc_land_sso_${grid}_${year}.nc
rm ${filename}_elevation.tmp ${filename}_orostd.tmp ${filename}_orosig.tmp ${filename}_orogam.tmp ${filename}_orothe.tmp \
    ${filename}_oropic.tmp ${filename}_oroval.tmp ${filename}_oromea.tmp
# initial soil conditions
${cdo} setmisstoc,280 ${filename}_surf_temp.tmp z.nc ; mv z.nc ${filename}_surf_temp.tmp
${cdo} merge ${filename}_surf_temp.tmp ${filename}_init_moist.tmp ${filename}_layer_moist.tmp ${filename}_snow.tmp ic_land_soil_${grid}_${year}.nc
ncrename -d cell,ncells ic_land_soil_${grid}_${year}.nc
rm ${filename}_surf_temp.tmp ${filename}_init_moist.tmp ${filename}_layer_moist.tmp ${filename}_snow.tmp
# uuid as an attribute
uuidOfHGrid=$(ncdump -h ${gridfile} | grep ':uuidOfHGrid = ' | cut -f2 -d'"')
ncatted -a uuidOfHGrid,global,o,c,${uuidOfHGrid} bc_land_frac_${grid}_${year}.nc
ncatted -a uuidOfHGrid,global,o,c,${uuidOfHGrid} bc_land_phys_${grid}_${year}.nc
ncatted -a uuidOfHGrid,global,o,c,${uuidOfHGrid} bc_land_soil_${grid}_${year}.nc
ncatted -a uuidOfHGrid,global,o,c,${uuidOfHGrid} bc_land_sso_${grid}_${year}.nc
ncatted -a uuidOfHGrid,global,o,c,${uuidOfHGrid} ic_land_soil_${grid}_${year}.nc

# clean up
rm ${jsbachfile}
rm gauss_*.tmp
rm ${filename}_*.tmp1
rm slf_${grid}.nc  slm_${res_atm}.nc
rm glac_${res_atm}.nc glac_${grid}.nc glac_inv_${grid}.nc
rm rmp_${res_atm}_to_${grid}
rm jan_surf_tmp.nc jan_surf.nc
rm ${soilparamfile}
#rm integer_slm.nc

# test: remap to regular grid to check data with ncview
[[ -f rmp_${grid}_to_t255 ]] || ${cdo} gencon,t255grid bc_land_frac_${grid}_${year}.nc rmp_${grid}_to_t255
${cdo} remap,t255grid,rmp_${grid}_to_t255 bc_land_frac_${grid}_${year}.nc test_bc_land_frac.t255.nc
${cdo} remap,t255grid,rmp_${grid}_to_t255 bc_land_phys_${grid}_${year}.nc test_bc_land_phys.t255.nc
${cdo} remap,t255grid,rmp_${grid}_to_t255 bc_land_soil_${grid}_${year}.nc test_bc_land_soil.t255.nc
${cdo} remap,t255grid,rmp_${grid}_to_t255 bc_land_sso_${grid}_${year}.nc test_bc_land_sso.t255.nc
${cdo} remap,t255grid,rmp_${grid}_to_t255 ic_land_soil_${grid}_${year}.nc test_ic_land_soil.t255.nc
rm rmp_${grid}_to_t255

[[ -f rmp_${grid}_to_t63 ]] || ${cdo} gencon,t63grid bc_land_frac_${grid}_${year}.nc rmp_${grid}_to_t63
${cdo} remap,t63grid,rmp_${grid}_to_t63 bc_land_frac_${grid}_${year}.nc test_bc_land_frac.t63.nc
${cdo} remap,t63grid,rmp_${grid}_to_t63 bc_land_phys_${grid}_${year}.nc test_bc_land_phys.t63.nc
${cdo} remap,t63grid,rmp_${grid}_to_t63 bc_land_soil_${grid}_${year}.nc test_bc_land_soil.t63.nc
${cdo} remap,t63grid,rmp_${grid}_to_t63 bc_land_sso_${grid}_${year}.nc test_bc_land_sso.t63.nc
${cdo} remap,t63grid,rmp_${grid}_to_t63 ic_land_soil_${grid}_${year}.nc test_ic_land_soil.t63.nc
rm rmp_${grid}_to_t63

exit
