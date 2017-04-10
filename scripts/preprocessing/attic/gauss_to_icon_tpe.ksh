#!/bin/ksh
#------------------------------------------------------------------------------
# generate jsbach initial boundary conditions and initial conditions for icon grid
# Veronika Gayler, MPI-M, Hamburg, 2013
#------------------------------------------------------------------------------
set -e

year=${year_cf}
[[ ${echam_fractional} = true ]] && fractional=fractional_ || fractional=""

# input files
jsbachfile=${srcdir}/input/jsbach/${grid}/jsbach_${res_atm}${res_oce}_${fractional}${ntiles}tiles_5layers_${year}.nc
ifile_topo=${srcdir}/input/echam6/${res_atm}${res_oce}_jan_surf.nc
filename=${grid}_jsbach_${ntiles}tiles_${year}

cdo -f nc2 copy ${pool_soil}/${soilparamfile} ${soilparamfile}

# preparations
#--------------

cdo -f nc2 copy ${jsbachfile} ${jsbachfile}_tmp && mv ${jsbachfile}_tmp ${jsbachfile}
cdo -f nc2 copy ${ifile_topo} ${ifile_topo}_tmp && mv ${ifile_topo}_tmp ${ifile_topo}

# generate remapping matrices if not yet available
[[ -f rmp_${res_atm}_to_${grid} ]] || cdo gendis,${gridfile} ${jsbachfile} rmp_${res_atm}_to_${grid}

# land sea mask on Gaussian grid
cdo selvar,slm ${jsbachfile} slm_${res_atm}_tmp.nc
# modification for terra planet
cdo addc,1. -mulc,0. slm_${res_atm}_tmp.nc slm_${res_atm}.nc
rm slm_${res_atm}_tmp.nc

# glacier mask on Gaussian grido
cdo selvar,glac ${jsbachfile} glac_${res_atm}.nc

# land sea mask on icon grid
cdo gtc,-3 -selvar,cell_sea_land_mask ${gridfile} slm_${grid}.nc

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
            cdo -fillmiss -setmissval,0 -selvar,${var} ${jsbachfile} gauss_${var}.tmp
            rmp=rmp_${res_atm}_to_${grid} ;;
        albedo )
            cdo -fillmiss -setmissval,0.07 -selvar,${var} ${jsbachfile} gauss_${var}.tmp
            rmp=rmp_${res_atm}_to_${grid} ;;
        maxmoist )
            cdo -fillmiss -setmissval,1.e-13 -selvar,${var} ${jsbachfile} gauss_${var}.tmp
            rmp=rmp_${res_atm}_to_${grid} ;;
        cover_fract | natural_veg )
            cdo splitlevel -selvar,${var} ${jsbachfile} ${var}
            cdo sub ${var}000001.nc2 glac_${res_atm}.nc ${var}000001.tmp
            mv ${var}000001.tmp ${var}000001.nc2
            rm -f gauss_${var}.tmp1
            cdo merge ${var}0000??.nc2 gauss_${var}.tmp1
            rm ${var}0000??.nc2
            cdo setvar,${var} -fillmiss -ifnotthen glac_${res_atm}.nc -ifthen slm_${res_atm}.nc \
                gauss_${var}.tmp1 gauss_${var}.tmp
            rm gauss_${var}.tmp1
            rmp=rmp_${res_atm}_to_${grid} ;;
        cover_type )
            cdo splitlevel -selvar,${var} ${jsbachfile} ${var}
            cdo add ${var}000001.nc2 glac_${res_atm}.nc ${var}000001.tmp
            mv ${var}000001.tmp ${var}000001.nc2
            rm -f gauss_${var}.tmp1
            cdo merge ${var}0000??.nc2 gauss_${var}.tmp1
            rm ${var}0000??.nc2
            cdo setvar,${var} -fillmiss -ifthen slm_${res_atm}.nc gauss_${var}.tmp1 gauss_${var}.tmp
            rm gauss_${var}.tmp1
            rmp=rmp_${res_atm}_to_${grid} ;;
        forest_fract | glac | snow )
            cdo setvar,${var} -fillmiss -ifthen slm_${res_atm}.nc -selvar,${var} ${jsbachfile} \
		gauss_${var}.tmp
            rmp=rmp_${res_atm}_to_${grid} ;;
        lai_clim | veg_fract | veg_ratio_max | roughness_length_oro )
            cdo setvar,${var} -fillmiss -sub -setmissval,-9e33 -setmisstoc,-9e33 -setmissval,0 \
                         -add slm_${res_atm}.nc -selvar,${var} ${jsbachfile} \
                         slm_${res_atm}.nc gauss_${var}.tmp
            rmp=rmp_${res_atm}_to_${grid} ;;
        albedo_veg_vis | albedo_veg_nir | albedo_soil_vis | albedo_soil_nir | surf_temp )
            cdo selvar,${var} ${jsbachfile} gauss_${var}.tmp
            rmp=rmp_${res_atm}_to_${grid} ;;
    esac

    # do the remapping 
    # Soil parameter fields are just read in (and renamed) from already remapped soil parameter file
    case ${var} in

 	elevation )
            cdo setvar,elevation -ifthen slm_${grid}.nc -selvar,cell_elevation ${gridfile} \
                       ${filename}_elevation.tmp ;;
	slf ) 
            rm -f ${filename}_${var}.tmp1 ${filename}_${var}.tmp ;;
	slm ) 
            cdo setvar,slm slm_${grid}.nc ${filename}_slm.tmp ;;
        fao )
            cdo remapnn,${gridfile} gauss_${var}.tmp ${filename}_${var}.tmp1
            cdo setvar,${var} -setmissval,-9e+33 -ifthen slm_${grid}.nc ${filename}_${var}.tmp1 \
                ${filename}_${var}.tmp ;;
        maxmoist )
            [[ ${dynveg} = true ]] && \
              cdo setvar,${var} -setmissval,-9e+33 -ifthen slm_${grid}.nc -selvar,wcap_dynveg \
                ${soilparamfile} ${filename}_${var}.tmp || \
	      cdo setvar,${var} -setmissval,-9e+33 -ifthen slm_${grid}.nc -selvar,wcap \
                ${soilparamfile} ${filename}_${var}.tmp
            ;;
	soil_depth )
            cdo setvar,${var} -setmissval,-9e+33 -ifthen slm_${grid}.nc -selvar,soildepth ${soilparamfile} ${filename}_${var}.tmp ;;
        bclapp )
            cdo setvar,${var} -setmissval,-9e+33 -ifthen slm_${grid}.nc -selvar,bclapp ${soilparamfile} ${filename}_${var}.tmp ;;
        soil_field_cap )
            cdo setvar,${var} -setmissval,-9e+33 -ifthen slm_${grid}.nc -selvar,volfc ${soilparamfile} ${filename}_${var}.tmp ;;
        heat_capacity )
            cdo setvar,${var} -setmissval,-9e+33 -ifthen slm_${grid}.nc -selvar,heatcapacity ${soilparamfile} ${filename}_${var}.tmp ;;
        heat_conductivity )
            cdo setvar,${var} -setmissval,-9e+33 -ifthen slm_${grid}.nc -selvar,heatcond ${soilparamfile} ${filename}_${var}.tmp ;;
        hyd_cond_sat )
            cdo setvar,${var} -setmissval,-9e+33 -ifthen slm_${grid}.nc -selvar,hydcond ${soilparamfile} ${filename}_${var}.tmp ;;
        moisture_pot )
            cdo setvar,${var} -setmissval,-9e+33 -ifthen slm_${grid}.nc -selvar,matrixpot ${soilparamfile} ${filename}_${var}.tmp ;;
        pore_size_index )
            cdo setvar,${var} -setmissval,-9e+33 -ifthen slm_${grid}.nc -selvar,psi_lambda ${soilparamfile} ${filename}_${var}.tmp ;;
        soil_porosity )
            cdo setvar,${var} -setmissval,-9e+33 -ifthen slm_${grid}.nc -selvar,volporosity ${soilparamfile} ${filename}_${var}.tmp ;;
        wilting_point )
            cdo setvar,${var} -setmissval,-9e+33 -ifthen slm_${grid}.nc -selvar,volpwp ${soilparamfile} ${filename}_${var}.tmp ;;
        * )
            cdo remap,${gridfile},${rmp} gauss_${var}.tmp ${filename}_${var}.tmp1
 	    cdo setvar,${var} -setmissval,-9e+33 -ifthen slm_${grid}.nc ${filename}_${var}.tmp1 \
		${filename}_${var}.tmp ;;
    esac
done

# treatment of glaciers
cdo gec,0.5 ${filename}_glac.tmp glac_${grid}.nc
cdo ltc,0.5 ${filename}_glac.tmp glac_inv_${grid}.nc
 
for var in ${varlist}; do
    case ${var} in

       cover_fract | natural_veg )
           cdo splitlevel ${filename}_${var}.tmp ${filename}_${var}_
	   typeset -Z2 nt=01
	   while [[ ${nt} -le ${ntiles} ]]; do
	       cdo mul ${filename}_${var}_0000${nt}.nc2 glac_inv_${grid}.nc \
		   ${filename}_${var}_${nt}.tmp
	       rm ${filename}_${var}_0000${nt}.nc2
	       (( nt = nt + 1 ))
	   done
	   cdo add ${filename}_${var}_01.tmp glac_${grid}.nc ${filename}_${var}_01.tmp1
           mv ${filename}_${var}_01.tmp1 ${filename}_${var}_01.tmp
           rm -f ${filename}_${var}.tmp
           cdo merge ${filename}_${var}_??.tmp ${filename}_${var}.tmp
           cdo -f nc2 setvar,${var} ${filename}_${var}.tmp ${filename}_${var}.tmp2
           rm ${filename}_${var}_??.tmp ;; 
       cover_type )
           cdo splitlevel ${filename}_${var}.tmp ${filename}_${var}_
	   typeset -Z2 nt=01
           cdo sub ${filename}_${var}_000001.nc2 glac_${grid}.nc ${filename}_${var}_01.tmp1
           mv ${filename}_${var}_01.tmp1 ${filename}_${var}_000001.nc2
           cdo nint ${filename}_${var}_000011.nc2 ${filename}_${var}_11.tmp1
           mv ${filename}_${var}_11.tmp1 ${filename}_${var}_000011.nc2
           rm -f ${filename}_${var}.tmp
           cdo merge ${filename}_${var}_0000??.nc2 ${filename}_${var}.tmp2
           rm ${filename}_${var}_0000??.nc2 ;; 
       glac )
	    cp glac_${grid}.nc ${filename}_${var}.tmp2 ;;
       lai_clim | bclapp | soil_field_cap | forest_fract | heat_capacity | heat_conductivity | root_depth | \
	   hyd_cond_sat | init_moist | layer_moist | maxmoist | moisture_pot | \
           pore_size_index | soil_porosity | snow | soil_depth | veg_fract | veg_ratio_max | wilting_point )
            cdo mul ${filename}_${var}.tmp glac_inv_${grid}.nc ${filename}_${var}.tmp2 ;;
       albedo )
            cdo setvar,${var} -ifthenelse glac_inv_${grid}.nc ${filename}_${var}.tmp \
		-mulc,0.7 glac_${grid}.nc ${filename}_${var}.tmp2 ;;
    esac
    [[ -f ${filename}_${var}.tmp2 ]] && mv ${filename}_${var}.tmp2  ${filename}_${var}.tmp
done

cdo mul slm_${grid}.nc -ltc,1 glac_${grid}.nc ${filename}_landnoglacier.tmp

cdo ifthenelse ${filename}_landnoglacier.tmp -setrtoc,0,0.1,0.1 ${filename}_soil_depth.tmp ${filename}_soil_depth.tmp ${filename}_soil_depth.tmp2
mv ${filename}_soil_depth.tmp2 ${filename}_soil_depth.tmp

cdo mul ${filename}_soil_field_cap.tmp ${filename}_soil_depth.tmp ${filename}_max_soilwater.tmp

cdo min ${filename}_maxmoist.tmp ${filename}_max_soilwater.tmp ${filename}_maxmoist.tmp2
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

cdo remapnn,${gridfile} ${ifile_topo} jan_surf_tmp.nc

# 2. land sea mask on the ICON grid
# ---------------------------------

# integer sea land mask from the icon grid file
# - ocean cell without land edge : -2
# - ocean cell with land edge    : -1
# - land cell with ocean edge    :  1
# - land cell without ocean edge :  2
#
cdo -f nc2 selname,cell_sea_land_mask ${gridfile} integer_slm.nc
#
# floating point sea mask
# - at land : 0.
# - at sea  : 1.
#
cdo -f nc2 -b F64 chname,cell_sea_land_mask,sea -ltc,0 integer_slm.nc ${filename}_sea.tmp
#
# floating point land mask
# - at land : 1.
# - at sea  : 0.
#
cdo -f nc2 -b F64 chname,cell_sea_land_mask,notsea -gtc,0 integer_slm.nc ${filename}_notsea.tmp

# 3. set all land parameters to 0 where at sea
# --------------------------------------------
# set surface land variables to zero at sea
# - at land : keep value
# - at sea  : set to 0
#
cdo mul jan_surf_tmp.nc ${filename}_notsea.tmp jan_surf.nc

# 4. separate variables transfered to new files
# ---------------------------------------------
#
# lake surface fractions
#
cdo chname,ALAKE,lake -selname,ALAKE jan_surf.nc ${filename}_lake.tmp
#
# land fraction (not sea and not lake)
#
cdo chname,notsea,land -sub ${filename}_notsea.tmp ${filename}_lake.tmp ${filename}_land.tmp
#
# surface physical properties
#
cdo chname,AZ0,z0  -selname,AZ0 jan_surf.nc ${filename}_z0.tmp
#
# orography
#
cdo chname,OROSTD,orostd -selname,OROSTD jan_surf.nc ${filename}_orostd.tmp
cdo chname,OROSIG,orosig -selname,OROSIG jan_surf.nc ${filename}_orosig.tmp
cdo chname,OROGAM,orogam -selname,OROGAM jan_surf.nc ${filename}_orogam.tmp
cdo chname,OROTHE,orothe -selname,OROTHE jan_surf.nc ${filename}_orothe.tmp
cdo chname,OROPIC,oropic -selname,OROPIC jan_surf.nc ${filename}_oropic.tmp
cdo chname,OROVAL,oroval -selname,OROVAL jan_surf.nc ${filename}_oroval.tmp
cdo chname,OROMEA,oromea -selname,OROMEA jan_surf.nc ${filename}_oromea.tmp

# generate four separate files for boundary conditions (bc) and one file for initial conditions (ic)

rm -f ${grid}_jsbach_clim.nc ${grid}_jsbach_${ntiles}tiles_${year}.nc ${grid}_jsbach.nc

# modifications for terra planet
# surface definitions
cdo mulc,0. ${filename}_glac.tmp ${filename}_glac.tmp2
cdo mulc,0. ${filename}_sea.tmp ${filename}_sea.tmp2
cdo addc,1. -mulc,0. ${filename}_notsea.tmp ${filename}_notsea.tmp2
cdo mulc,0. ${filename}_lake.tmp ${filename}_lake.tmp2
cdo addc,1. -mulc,0. ${filename}_land.tmp ${filename}_land.tmp2
cdo splitlevel ${filename}_cover_fract.tmp ${filename}_cover_fract
cdo addc,1.e-10 -mulc,0. ${filename}_cover_fract000001.nc2 ${filename}_cover_fract000001_tmp.nc
cdo addc,1.e-10 -mulc,0. ${filename}_cover_fract000002.nc2 ${filename}_cover_fract000002_tmp.nc
cdo addc,1.e-10 -mulc,0. ${filename}_cover_fract000003.nc2 ${filename}_cover_fract000003_tmp.nc
cdo addc,1.e-10 -mulc,0. ${filename}_cover_fract000004.nc2 ${filename}_cover_fract000004_tmp.nc
cdo addc,1.e-10 -mulc,0. ${filename}_cover_fract000005.nc2 ${filename}_cover_fract000005_tmp.nc
cdo addc,1.e-10 -mulc,0. ${filename}_cover_fract000006.nc2 ${filename}_cover_fract000006_tmp.nc
cdo subc,1.e-09 -addc,1. -mulc,0. ${filename}_cover_fract000007.nc2 ${filename}_cover_fract000007_tmp.nc
cdo addc,1.e-10 -mulc,0. ${filename}_cover_fract000008.nc2 ${filename}_cover_fract000008_tmp.nc
cdo addc,1.e-10 -mulc,0. ${filename}_cover_fract000009.nc2 ${filename}_cover_fract000009_tmp.nc
cdo addc,1.e-10 -mulc,0. ${filename}_cover_fract000010.nc2 ${filename}_cover_fract000010_tmp.nc
cdo addc,1.e-10 -mulc,0. ${filename}_cover_fract000011.nc2 ${filename}_cover_fract000011_tmp.nc
cdo merge ${filename}_cover_fract0000??_tmp.nc ${filename}_cover_fract.tmp2
rm ${filename}_cover_fract0000??.nc2 ${filename}_cover_fract0000??_tmp.nc
cdo splitlevel ${filename}_cover_type.tmp ${filename}_cover_type
cdo addc,2. -mulc,0. ${filename}_cover_type000001.nc2 ${filename}_cover_type000001_tmp.nc
cdo addc,3. -mulc,0. ${filename}_cover_type000002.nc2 ${filename}_cover_type000002_tmp.nc
cdo addc,4. -mulc,0. ${filename}_cover_type000003.nc2 ${filename}_cover_type000003_tmp.nc
cdo addc,5. -mulc,0. ${filename}_cover_type000004.nc2 ${filename}_cover_type000004_tmp.nc
cdo addc,10. -mulc,0. ${filename}_cover_type000005.nc2 ${filename}_cover_type000005_tmp.nc
cdo addc,11. -mulc,0. ${filename}_cover_type000006.nc2 ${filename}_cover_type000006_tmp.nc
cdo addc,12. -mulc,0. ${filename}_cover_type000007.nc2 ${filename}_cover_type000007_tmp.nc
cdo addc,13. -mulc,0. ${filename}_cover_type000008.nc2 ${filename}_cover_type000008_tmp.nc
cdo addc,15. -mulc,0. ${filename}_cover_type000009.nc2 ${filename}_cover_type000009_tmp.nc
cdo addc,16. -mulc,0. ${filename}_cover_type000010.nc2 ${filename}_cover_type000010_tmp.nc
cdo addc,20. -mulc,0. ${filename}_cover_type000011.nc2 ${filename}_cover_type000011_tmp.nc
cdo merge ${filename}_cover_type0000??_tmp.nc ${filename}_cover_type.tmp2
rm ${filename}_cover_type0000??.nc2 ${filename}_cover_type0000??_tmp.nc
cdo splitlevel ${filename}_natural_veg.tmp ${filename}_natural_veg
cdo addc,1.e-10 -mulc,0. ${filename}_natural_veg000001.nc2 ${filename}_natural_veg000001_tmp.nc
cdo addc,1.e-10 -mulc,0. ${filename}_natural_veg000002.nc2 ${filename}_natural_veg000002_tmp.nc
cdo addc,1.e-10 -mulc,0. ${filename}_natural_veg000003.nc2 ${filename}_natural_veg000003_tmp.nc
cdo addc,1.e-10 -mulc,0. ${filename}_natural_veg000004.nc2 ${filename}_natural_veg000004_tmp.nc
cdo addc,1.e-10 -mulc,0. ${filename}_natural_veg000005.nc2 ${filename}_natural_veg000005_tmp.nc
cdo addc,1.e-10 -mulc,0. ${filename}_natural_veg000006.nc2 ${filename}_natural_veg000006_tmp.nc
cdo subc,1.e-09 -addc,1. -mulc,0. ${filename}_natural_veg000007.nc2 ${filename}_natural_veg000007_tmp.nc
cdo addc,1.e-10 -mulc,0. ${filename}_natural_veg000008.nc2 ${filename}_natural_veg000008_tmp.nc
cdo addc,1.e-10 -mulc,0. ${filename}_natural_veg000009.nc2 ${filename}_natural_veg000009_tmp.nc
cdo addc,1.e-10 -mulc,0. ${filename}_natural_veg000010.nc2 ${filename}_natural_veg000010_tmp.nc
cdo addc,1.e-10 -mulc,0. ${filename}_natural_veg000011.nc2 ${filename}_natural_veg000011_tmp.nc
cdo merge ${filename}_natural_veg0000??_tmp.nc ${filename}_natural_veg.tmp2
rm ${filename}_natural_veg0000??.nc2 ${filename}_natural_veg0000??_tmp.nc
cdo addc,1. -mulc,0. ${filename}_veg_ratio_max.tmp ${filename}_veg_ratio_max.tmp2
rm ${filename}_glac.tmp ${filename}_sea.tmp ${filename}_notsea.tmp ${filename}_lake.tmp ${filename}_land.tmp ${filename}_cover_fract.tmp ${filename}_cover_type.tmp ${filename}_natural_veg.tmp ${filename}_veg_ratio_max.tmp 
# physical properties
cdo addc,3. -mulc,0. ${filename}_lai_clim.tmp ${filename}_lai_clim.tmp2
cdo addc,0.9 -mulc,0. ${filename}_veg_fract.tmp ${filename}_veg_fract.tmp2
cdo addc,0.05 -mulc,0. ${filename}_roughness_length.tmp ${filename}_roughness_length.tmp2
cdo addc,0.04 -mulc,0. ${filename}_roughness_length_oro.tmp ${filename}_roughness_length_oro.tmp2
cdo addc,0.05 -mulc,0. ${filename}_z0.tmp ${filename}_z0.tmp2
cdo addc,0.12 -mulc,0. ${filename}_albedo.tmp ${filename}_albedo.tmp2
cdo addc,0.12 -mulc,0. ${filename}_albedo_veg_vis.tmp ${filename}_albedo_veg_vis.tmp2
cdo addc,0.12 -mulc,0. ${filename}_albedo_veg_nir.tmp ${filename}_albedo_veg_nir.tmp2
cdo addc,0.12 -mulc,0. ${filename}_albedo_soil_vis.tmp ${filename}_albedo_soil_vis.tmp2
cdo addc,0.12 -mulc,0. ${filename}_albedo_soil_nir.tmp ${filename}_albedo_soil_nir.tmp2
cdo -mulc,0. ${filename}_forest_fract.tmp ${filename}_forest_fract.tmp2
rm ${filename}_lai_clim.tmp ${filename}_veg_fract.tmp ${filename}_roughness_length.tmp ${filename}_roughness_length_oro.tmp ${filename}_z0.tmp ${filename}_albedo.tmp ${filename}_albedo_veg_vis.tmp ${filename}_albedo_veg_nir.tmp ${filename}_albedo_soil_vis.tmp ${filename}_albedo_soil_nir.tmp ${filename}_forest_fract.tmp
# soil properties
cdo addc,3. -mulc,0. ${filename}_fao.tmp ${filename}_fao.tmp2
cdo addc,0.4721088 -mulc,0. ${filename}_maxmoist.tmp ${filename}_maxmoist.tmp2
cdo addc,8.32 -mulc,0. ${filename}_bclapp.tmp ${filename}_bclapp.tmp2
cdo addc,0.3216 -mulc,0. ${filename}_soil_field_cap.tmp ${filename}_soil_field_cap.tmp2
cdo addc,2.26039e+06 -mulc,0. ${filename}_heat_capacity.tmp ${filename}_heat_capacity.tmp2
cdo addc,5.7709 -mulc,0. ${filename}_heat_conductivity.tmp ${filename}_heat_conductivity.tmp2
cdo addc,2.845e-06 -mulc,0. ${filename}_hyd_cond_sat.tmp ${filename}_hyd_cond_sat.tmp2
cdo addc,0.28893 -mulc,0. ${filename}_moisture_pot.tmp ${filename}_moisture_pot.tmp2
cdo addc,0.194 -mulc,0. ${filename}_pore_size_index.tmp ${filename}_pore_size_index.tmp2
cdo addc,0.4487 -mulc,0. ${filename}_soil_porosity.tmp ${filename}_soil_porosity.tmp2
cdo addc,6.0 -mulc,0. ${filename}_soil_depth.tmp ${filename}_soil_depth.tmp2
cdo addc,0.187 -mulc,0. ${filename}_wilting_point.tmp ${filename}_wilting_point.tmp2
cdo addc,1.468 -mulc,0. ${filename}_root_depth.tmp ${filename}_root_depth.tmp2
rm ${filename}_fao.tmp ${filename}_maxmoist.tmp ${filename}_bclapp.tmp ${filename}_soil_field_cap.tmp ${filename}_heat_capacity.tmp ${filename}_heat_conductivity.tmp ${filename}_hyd_cond_sat.tmp ${filename}_moisture_pot.tmp ${filename}_pore_size_index.tmp ${filename}_soil_porosity.tmp ${filename}_soil_depth.tmp ${filename}_wilting_point.tmp ${filename}_root_depth.tmp
# topographie
cdo gtc,0. ${filename}_elevation.tmp ${filename}_elevation_mask.tmp
cdo mul ${filename}_elevation.tmp ${filename}_elevation_mask.tmp ${filename}_elevation.tmp2
#cdo mulc,0. ${filename}_elevation.tmp ${filename}_elevation.tmp2
#cdo addc,1. -mulc,0. ${filename}_orostd.tmp ${filename}_orostd.tmp2
#cdo addc,0.001 -mulc,0. ${filename}_orosig.tmp ${filename}_orosig.tmp2
#cdo addc,0.5 -mulc,0. ${filename}_orogam.tmp ${filename}_orogam.tmp2
#cdo mulc,0. ${filename}_orothe.tmp ${filename}_orothe.tmp2
#cdo addc,1. -mulc,0. ${filename}_oropic.tmp ${filename}_oropic.tmp2
#cdo addc,-1. -mulc,0. ${filename}_oroval.tmp ${filename}_oroval.tmp2
#cdo mulc,0. ${filename}_oromea.tmp ${filename}_oromea.tmp2
rm ${filename}_elevation.tmp ${filename}_elevation_mask.tmp 
#rm ${filename}_elevation.tmp ${filename}_orostd.tmp ${filename}_orosig.tmp ${filename}_orogam.tmp ${filename}_orothe.tmp ${filename}_oropic.tmp ${filename}_oroval.tmp ${filename}_oromea.tmp
# initial soil conditions
cdo addc,0.2460544 -mulc,0. ${filename}_init_moist.tmp ${filename}_init_moist.tmp2
cdo splitlevel ${filename}_layer_moist.tmp ${filename}_layer_moist
cdo addc,0.010452 -mulc,0. ${filename}_layer_moist00.065.nc2 ${filename}_layer_moist000001_tmp.nc
cdo addc,0.0408432 -mulc,0. ${filename}_layer_moist00.319.nc2 ${filename}_layer_moist000002_tmp.nc
cdo addc,0.1468104 -mulc,0. ${filename}_layer_moist01.232.nc2 ${filename}_layer_moist000003_tmp.nc
cdo addc,0.0379488 -mulc,0. ${filename}_layer_moist04.134.nc2 ${filename}_layer_moist000004_tmp.nc
cdo addc,0.01 -mulc,0. ${filename}_layer_moist09.834.nc2 ${filename}_layer_moist000005_tmp.nc
cdo merge ${filename}_layer_moist00000?_tmp.nc ${filename}_layer_moist.tmp2
rm ${filename}_layer_moist0?.???.nc2 ${filename}_layer_moist00000?_tmp.nc
cdo mulc,0. ${filename}_snow.tmp ${filename}_snow.tmp2
rm ${filename}_init_moist.tmp ${filename}_layer_moist.tmp ${filename}_snow.tmp

# define surfaces
cdo merge ${filename}_slm.tmp ${filename}_glac.tmp2 ${filename}_sea.tmp2 ${filename}_notsea.tmp2 ${filename}_lake.tmp2 ${filename}_land.tmp2 ${filename}_cover_fract.tmp2 ${filename}_cover_type.tmp2 ${filename}_natural_veg.tmp2 ${filename}_veg_ratio_max.tmp2 bc_land_frac_${grid}_${year}.nc
rm ${filename}_slm.tmp ${filename}_glac.tmp2 ${filename}_sea.tmp2 ${filename}_notsea.tmp2 ${filename}_lake.tmp2 ${filename}_land.tmp2 ${filename}_cover_fract.tmp2 ${filename}_cover_type.tmp2 ${filename}_natural_veg.tmp2 ${filename}_veg_ratio_max.tmp2 
# physical properties
cdo merge ${filename}_lai_clim.tmp2 ${filename}_veg_fract.tmp2 ${filename}_roughness_length.tmp2 ${filename}_roughness_length_oro.tmp2 ${filename}_z0.tmp2 ${filename}_albedo.tmp2 ${filename}_albedo_veg_vis.tmp2 ${filename}_albedo_veg_nir.tmp2 ${filename}_albedo_soil_vis.tmp2 ${filename}_albedo_soil_nir.tmp2 ${filename}_forest_fract.tmp2 bc_land_phys_${grid}_${year}.nc
rm ${filename}_lai_clim.tmp2 ${filename}_veg_fract.tmp2 ${filename}_roughness_length.tmp2 ${filename}_roughness_length_oro.tmp2 ${filename}_z0.tmp2 ${filename}_albedo.tmp2 ${filename}_albedo_veg_vis.tmp2 ${filename}_albedo_veg_nir.tmp2 ${filename}_albedo_soil_vis.tmp2 ${filename}_albedo_soil_nir.tmp2 ${filename}_forest_fract.tmp2
# soil properties
cdo merge ${filename}_fao.tmp2 ${filename}_maxmoist.tmp2 ${filename}_bclapp.tmp2 ${filename}_soil_field_cap.tmp2 ${filename}_heat_capacity.tmp2 ${filename}_heat_conductivity.tmp2 ${filename}_hyd_cond_sat.tmp2 ${filename}_moisture_pot.tmp2 ${filename}_pore_size_index.tmp2 ${filename}_soil_porosity.tmp2 ${filename}_soil_depth.tmp2 ${filename}_wilting_point.tmp2 ${filename}_root_depth.tmp2 bc_land_soil_${grid}_${year}.nc
rm ${filename}_fao.tmp2 ${filename}_maxmoist.tmp2 ${filename}_bclapp.tmp2 ${filename}_soil_field_cap.tmp2 ${filename}_heat_capacity.tmp2 ${filename}_heat_conductivity.tmp2 ${filename}_hyd_cond_sat.tmp2 ${filename}_moisture_pot.tmp2 ${filename}_pore_size_index.tmp2 ${filename}_soil_porosity.tmp2 ${filename}_soil_depth.tmp2 ${filename}_wilting_point.tmp2 ${filename}_root_depth.tmp2
# topography
cdo merge ${filename}_elevation.tmp2 ${filename}_orostd.tmp ${filename}_orosig.tmp ${filename}_orogam.tmp ${filename}_orothe.tmp ${filename}_oropic.tmp ${filename}_oroval.tmp ${filename}_oromea.tmp bc_land_sso_${grid}_${year}.nc
rm ${filename}_elevation.tmp2 ${filename}_orostd.tmp ${filename}_orosig.tmp ${filename}_orogam.tmp ${filename}_orothe.tmp ${filename}_oropic.tmp ${filename}_oroval.tmp ${filename}_oromea.tmp
# initial soil conditions
cdo merge ${filename}_surf_temp.tmp ${filename}_init_moist.tmp2 ${filename}_layer_moist.tmp2 ${filename}_snow.tmp2 ic_land_soil_${grid}_${year}.nc
rm ${filename}_surf_temp.tmp ${filename}_init_moist.tmp2 ${filename}_layer_moist.tmp2 ${filename}_snow.tmp2
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
rm ${filename}_*.tmp1 ${filename}_landnoglacier.tmp
rm slm_${grid}.nc  slm_${res_atm}.nc
rm glac_${res_atm}.nc glac_${grid}.nc glac_inv_${grid}.nc
rm rmp_${res_atm}_to_${grid}
rm jan_surf_tmp.nc jan_surf.nc
rm ${soilparamfile}
rm integer_slm.nc

# test: remap to regular grid to check data with ncview
[[ -f rmp_${grid}_to_t255 ]] || cdo gencon,t255grid bc_land_frac_${grid}_${year}.nc rmp_${grid}_to_t255
cdo remap,t255grid,rmp_${grid}_to_t255 bc_land_frac_${grid}_${year}.nc test_bc_land_frac.t255.nc
cdo remap,t255grid,rmp_${grid}_to_t255 bc_land_phys_${grid}_${year}.nc test_bc_land_phys.t255.nc
cdo remap,t255grid,rmp_${grid}_to_t255 bc_land_soil_${grid}_${year}.nc test_bc_land_soil.t255.nc
cdo remap,t255grid,rmp_${grid}_to_t255 bc_land_sso_${grid}_${year}.nc test_bc_land_sso.t255.nc
cdo remap,t255grid,rmp_${grid}_to_t255 ic_land_soil_${grid}_${year}.nc test_ic_land_soil.t255.nc
rm rmp_${grid}_to_t255

[[ -f rmp_${grid}_to_t63 ]] || cdo gencon,t63grid bc_land_frac_${grid}_${year}.nc rmp_${grid}_to_t63
cdo remap,t63grid,rmp_${grid}_to_t63 bc_land_frac_${grid}_${year}.nc test_bc_land_frac.t63.nc
cdo remap,t63grid,rmp_${grid}_to_t63 bc_land_phys_${grid}_${year}.nc test_bc_land_phys.t63.nc
cdo remap,t63grid,rmp_${grid}_to_t63 bc_land_soil_${grid}_${year}.nc test_bc_land_soil.t63.nc
cdo remap,t63grid,rmp_${grid}_to_t63 bc_land_sso_${grid}_${year}.nc test_bc_land_sso.t63.nc
cdo remap,t63grid,rmp_${grid}_to_t63 ic_land_soil_${grid}_${year}.nc test_ic_land_soil.t63.nc
rm rmp_${grid}_to_t63

exit
