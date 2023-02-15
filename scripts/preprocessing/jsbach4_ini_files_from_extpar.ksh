set -x

# file names and paths
path_extpar=/work/mj0060/m212070/test_extpar4jsbach/extpar4jsbach/work
file_extpar_orig=icon_extpar__20220907
file_extpar=${file_extpar_orig}_grid

path_grid=/pool/data/ICON/grids/public/mpim/0015
file_grid=icon_grid_0015_R02B09_G

path_bc=/work/mj0060/m212070/extpar4jsbach/work/R2B9
file_bc_frac=bc_land_frac_2015
file_bc_phys=bc_land_phys_2015
file_bc_soil=bc_land_soil_2015

# parameters
value_albedo_max=60. # maximum value of background albedo accepted [%]
fillvalue_albedo_glac=12.
fillvalue_albedo_coast=12.
lower_bound_NDVI_veg_fract=0.12
upper_bound_NDVI_veg_fract=0.92
upper_bound_veg_fract=0.98
fact_NDVI_veg_fract=1.225 # should be: upper_bound_veg_fract / (upper_bound_NDVI_veg_fract - lower_bound_NDVI_veg_fract)
lower_bound_NDVI_lai_clim=0.12
upper_bound_NDVI_lai_clim=0.8
upper_bound_lai_clim=6.
fact_NDVI_lai_clim=8.8235294 # should be: upper_bound_lai_clim / (upper_bound_NDVI_lai_clim - lower_bound_NDVI_lai_clim)
minimum_root_depth=0.15
delta=0.0000001 # small number

cdo=/sw/spack-levante/mambaforge-4.11.0-0-Linux-x86_64-sobz6z/bin/cdo
rm=/usr/bin/rm
cp=/usr/bin/cp
expr=/usr/bin/expr

# copy extpar file and set grid
#${cdo} setgrid,${path_grid}/${file_grid}.nc ${path_extpar}/${file_extpar_orig}.nc ${file_extpar}.nc

# copy boundary condition files and initial condition files
${cp} -p ${path_bc}/${file_bc_frac}.nc ${file_bc_frac}_orig.nc
${cp} -p ${path_bc}/${file_bc_phys}.nc ${file_bc_phys}_orig.nc
${cp} -p ${path_bc}/${file_bc_soil}.nc ${file_bc_soil}_orig.nc

# modify lake fraction
# select extpar fields of lake fraction (from flake) and land fraction (Here it is assumed that the land fraction is the part of the grid cell, which is land (i.e. not ocean).)
${cdo} selvar,FR_LAKE ${file_extpar}.nc ${file_extpar}_fr_lake.nc
${cdo} selvar,FR_LAND ${file_extpar}.nc ${file_extpar}_fr_land.nc
# select extpar field of lake fraction from Globcover (lake fractions for Caspian Sean and Aral Lake area, missing value elsewhere)
${cdo} masklonlatbox,46.,62.,36.,48. -sellevel,21 -selvar,LU_CLASS_FRACTION ${file_extpar}.nc ${file_extpar}_LU_CLASS_FRACTION_lev21_casp.nc
# mask: 1 for Caspian Sea and Aral Lake, missing value elsewhere
${cdo} gtc,-1. ${file_extpar}_LU_CLASS_FRACTION_lev21_casp.nc ${file_extpar}_LU_CLASS_FRACTION_mask_casp.nc
# mask: 0 for Caspian Sea and Aral Lake area, 1 elsewhere
${cdo} subc,1. -setmisstoc,2. ${file_extpar}_LU_CLASS_FRACTION_mask_casp.nc ${file_extpar}_LU_CLASS_FRACTION_mask_mul_casp.nc
# lake fraction (from flake) excluding lakes in the Caspian Sea and Aral Lake area
${cdo} mul ${file_extpar}_fr_lake.nc ${file_extpar}_LU_CLASS_FRACTION_mask_mul_casp.nc ${file_extpar}_fr_lake_clean.nc
# add Caspian Sea and Aral Lake from Globcover to lake fractions from flake elsewhere
${cdo} add ${file_extpar}_fr_lake_clean.nc -setmisstoc,0. ${file_extpar}_LU_CLASS_FRACTION_lev21_casp.nc ${file_extpar}_fr_lake_casp.nc
# mask pure land grid boxes with no ocean
${cdo} gec,1. -selvar,notsea ${file_bc_frac}_orig.nc ${file_extpar}_fr_land_mask.nc
# lake fraction just for grid boxes with no ocean ==> this is the new lake fraction
${cdo} mul ${file_extpar}_fr_lake_casp.nc ${file_extpar}_fr_land_mask.nc ${file_extpar}_fr_lake_nosea.nc
# new land fraction (This is the part of the surface that is handled by JSBACH4 and that is not lake (,i.e. this 'land' fraction is 1 over the ocean). It is not the land part of the grid cell.)
${cdo} sub -addc,1. -mulc,0. ${file_extpar}_fr_land.nc ${file_extpar}_fr_lake_nosea.nc ${file_extpar}_fr_land_no_lake.nc
# change variable names
${cdo} chname,FR_LAKE,fract_lake ${file_extpar}_fr_lake_nosea.nc ${file_extpar}_fr_lake_name.nc
${cdo} chname,FR_LAKE,lake ${file_extpar}_fr_lake_nosea.nc ${file_extpar}_lake_name.nc
${cdo} chname,FR_LAND,fract_land ${file_extpar}_fr_land_no_lake.nc ${file_extpar}_fr_land_name.nc
# replace old lake fraction and land fraction by extpar fields in bc_frac file
${cdo} replace ${file_bc_frac}_orig.nc ${file_extpar}_fr_lake_name.nc ${file_bc_frac}_lake.nc
${cdo} replace ${file_bc_frac}_lake.nc ${file_extpar}_fr_land_name.nc ${file_bc_frac}_lake_land.nc
# replace old lake fraction also for atmosphere input in bc_frac file
${cdo} replace ${file_bc_frac}_lake_land.nc ${file_extpar}_lake_name.nc ${file_bc_frac}_lake_land_lake.nc
# clean
${rm} ${file_extpar}_fr_lake.nc ${file_extpar}_fr_land.nc ${file_extpar}_LU_CLASS_FRACTION_lev21_casp.nc ${file_extpar}_LU_CLASS_FRACTION_mask_casp.nc ${file_extpar}_LU_CLASS_FRACTION_mask_mul_casp.nc ${file_extpar}_fr_lake_clean.nc ${file_extpar}_fr_lake_casp.nc ${file_bc_frac}_lake.nc ${file_extpar}_fr_lake_name.nc ${file_extpar}_fr_land_mask.nc ${file_extpar}_fr_lake_nosea.nc ${file_extpar}_lake_name.nc ${file_extpar}_fr_land_no_lake.nc ${file_extpar}_fr_land_name.nc ${file_bc_frac}_lake_land.nc

# modify glacier fraction
# select extpar field of glacier fraction (from Globcover)
${cdo} selvar,ICE ${file_extpar}.nc ${file_extpar}_ICE.nc
# calculate new vegetation fraction (1 - glacier)
${cdo} addc,1. -mulc,-1. ${file_extpar}_ICE.nc ${file_extpar}_1-ICE.nc
# change variable names
${cdo} chname,ICE,fract_glac ${file_extpar}_ICE.nc ${file_extpar}_ICE_fract_glac.nc
${cdo} chname,ICE,glac ${file_extpar}_ICE.nc ${file_extpar}_ICE_glac.nc
${cdo} chname,ICE,fract_veg ${file_extpar}_1-ICE.nc ${file_extpar}_ICE_fract_veg.nc
# replace old glacier fraction and vegetation fraction by extpar fields in bc_frac file
${cdo} replace ${file_bc_frac}_lake_land_lake.nc ${file_extpar}_ICE_fract_glac.nc ${file_bc_frac}_lake_land_lake_glac.nc
${cdo} replace ${file_bc_frac}_lake_land_lake_glac.nc ${file_extpar}_ICE_fract_veg.nc ${file_bc_frac}_lake_land_lake_glac_veg.nc
# replace old glacier fraction also for atmosphere input in bc_frac file
${cdo} replace ${file_bc_frac}_lake_land_lake_glac_veg.nc ${file_extpar}_ICE_glac.nc ${file_bc_frac}_extpar.nc
# clean
${rm} ${file_extpar}_ICE.nc ${file_extpar}_1-ICE.nc ${file_extpar}_ICE_fract_glac.nc ${file_extpar}_ICE_glac.nc ${file_extpar}_ICE_fract_veg.nc ${file_bc_frac}_lake_land_lake.nc ${file_bc_frac}_lake_land_lake_glac.nc ${file_bc_frac}_lake_land_lake_glac_veg.nc

# soil texture
# % ==> fraction, fill missing value, and add soil texture to bc_soil-file
${cdo} selvar,soil_depth ${file_bc_soil}_orig.nc ${file_bc_soil}_soil_depth.nc

${cdo} chname,soil_depth,FR_SAND -add -mulc,0. ${file_bc_soil}_soil_depth.nc -mulc,0.01 -selvar,FR_SAND ${file_extpar}.nc ${file_extpar}_FR_SAND.nc
${cdo} chname,soil_depth,FR_SILT -add -mulc,0. ${file_bc_soil}_soil_depth.nc -mulc,0.01 -selvar,FR_SILT ${file_extpar}.nc ${file_extpar}_FR_SILT.nc
${cdo} chname,soil_depth,FR_CLAY -add -mulc,0. ${file_bc_soil}_soil_depth.nc -mulc,0.01 -selvar,FR_CLAY ${file_extpar}.nc ${file_extpar}_FR_CLAY.nc
${cdo} gtc,-0.001 ${file_extpar}_FR_SAND.nc ${file_extpar}_FR_SAND_MASK.nc
${cdo} gtc,-0.001 ${file_extpar}_FR_SILT.nc ${file_extpar}_FR_SILT_MASK.nc
${cdo} gtc,-0.001 ${file_extpar}_FR_CLAY.nc ${file_extpar}_FR_CLAY_MASK.nc
${cdo} mul ${file_extpar}_FR_SAND_MASK.nc ${file_extpar}_FR_SILT_MASK.nc ${file_extpar}_SAND_SILT_MASK.nc
${cdo} mul ${file_extpar}_SAND_SILT_MASK.nc ${file_extpar}_FR_CLAY_MASK.nc ${file_extpar}_SAND_SILT_CLAY_MASK.nc
${cdo} mul ${file_extpar}_FR_SAND.nc ${file_extpar}_SAND_SILT_CLAY_MASK.nc ${file_extpar}_FR_SAND_0.nc
${cdo} mul ${file_extpar}_FR_SILT.nc ${file_extpar}_SAND_SILT_CLAY_MASK.nc ${file_extpar}_FR_SILT_0.nc
${cdo} mul ${file_extpar}_FR_CLAY.nc ${file_extpar}_SAND_SILT_CLAY_MASK.nc ${file_extpar}_FR_CLAY_0.nc
${cdo} add ${file_extpar}_FR_SAND_0.nc -mulc,-0.5 -addc,-1. ${file_extpar}_SAND_SILT_CLAY_MASK.nc ${file_extpar}_FR_SAND_CORRECTED.nc
${cdo} add ${file_extpar}_FR_SILT_0.nc -mulc,-0.25 -addc,-1. ${file_extpar}_SAND_SILT_CLAY_MASK.nc ${file_extpar}_FR_SILT_CORRECTED.nc
${cdo} add ${file_extpar}_FR_CLAY_0.nc -mulc,-0.25 -addc,-1. ${file_extpar}_SAND_SILT_CLAY_MASK.nc ${file_extpar}_FR_CLAY_CORRECTED.nc

${cdo} chname,soil_depth,SUB_FR_SAND -add -mulc,0. ${file_bc_soil}_soil_depth.nc -mulc,0.01 -selvar,SUB_FR_SAND ${file_extpar}.nc ${file_extpar}_SUB_FR_SAND.nc
${cdo} chname,soil_depth,SUB_FR_SILT -add -mulc,0. ${file_bc_soil}_soil_depth.nc -mulc,0.01 -selvar,SUB_FR_SILT ${file_extpar}.nc ${file_extpar}_SUB_FR_SILT.nc
${cdo} chname,soil_depth,SUB_FR_CLAY -add -mulc,0. ${file_bc_soil}_soil_depth.nc -mulc,0.01 -selvar,SUB_FR_CLAY ${file_extpar}.nc ${file_extpar}_SUB_FR_CLAY.nc
${cdo} gtc,-0.001 ${file_extpar}_SUB_FR_SAND.nc ${file_extpar}_SUB_FR_SAND_MASK.nc
${cdo} gtc,-0.001 ${file_extpar}_SUB_FR_SILT.nc ${file_extpar}_SUB_FR_SILT_MASK.nc
${cdo} gtc,-0.001 ${file_extpar}_SUB_FR_CLAY.nc ${file_extpar}_SUB_FR_CLAY_MASK.nc
${cdo} mul ${file_extpar}_SUB_FR_SAND_MASK.nc ${file_extpar}_SUB_FR_SILT_MASK.nc ${file_extpar}_SUB_SAND_SILT_MASK.nc
${cdo} mul ${file_extpar}_SUB_SAND_SILT_MASK.nc ${file_extpar}_SUB_FR_CLAY_MASK.nc ${file_extpar}_SUB_SAND_SILT_CLAY_MASK.nc
${cdo} mul ${file_extpar}_SUB_FR_SAND.nc ${file_extpar}_SUB_SAND_SILT_CLAY_MASK.nc ${file_extpar}_SUB_FR_SAND_0.nc
${cdo} mul ${file_extpar}_SUB_FR_SILT.nc ${file_extpar}_SUB_SAND_SILT_CLAY_MASK.nc ${file_extpar}_SUB_FR_SILT_0.nc
${cdo} mul ${file_extpar}_SUB_FR_CLAY.nc ${file_extpar}_SUB_SAND_SILT_CLAY_MASK.nc ${file_extpar}_SUB_FR_CLAY_0.nc
${cdo} add ${file_extpar}_SUB_FR_SAND_0.nc -mulc,-0.5 -addc,-1. ${file_extpar}_SUB_SAND_SILT_CLAY_MASK.nc ${file_extpar}_SUB_FR_SAND_CORRECTED.nc
${cdo} add ${file_extpar}_SUB_FR_SILT_0.nc -mulc,-0.25 -addc,-1. ${file_extpar}_SUB_SAND_SILT_CLAY_MASK.nc ${file_extpar}_SUB_FR_SILT_CORRECTED.nc
${cdo} add ${file_extpar}_SUB_FR_CLAY_0.nc -mulc,-0.25 -addc,-1. ${file_extpar}_SUB_SAND_SILT_CLAY_MASK.nc ${file_extpar}_SUB_FR_CLAY_CORRECTED.nc

${cdo} chname,soil_depth,FR_OC -add -mulc,0. ${file_bc_soil}_soil_depth.nc -mulc,0.01 -selvar,FR_OC ${file_extpar}.nc ${file_extpar}_FR_OC.nc
${cdo} gtc,-0.001 ${file_extpar}_FR_OC.nc ${file_extpar}_FR_OC_MASK.nc
${cdo} mul ${file_extpar}_FR_OC.nc ${file_extpar}_FR_OC_MASK.nc ${file_extpar}_FR_OC_CORRECTED.nc

${cdo} chname,soil_depth,SUB_FR_OC -add -mulc,0. ${file_bc_soil}_soil_depth.nc -mulc,0.01 -selvar,SUB_FR_OC ${file_extpar}.nc ${file_extpar}_SUB_FR_OC.nc
${cdo} gtc,-0.001 ${file_extpar}_SUB_FR_OC.nc ${file_extpar}_SUB_FR_OC_MASK.nc
${cdo} mul ${file_extpar}_SUB_FR_OC.nc ${file_extpar}_SUB_FR_OC_MASK.nc ${file_extpar}_SUB_FR_OC_CORRECTED.nc

${cdo} merge ${file_extpar}_FR_SAND_CORRECTED.nc ${file_extpar}_FR_SILT_CORRECTED.nc ${file_extpar}_FR_CLAY_CORRECTED.nc ${file_extpar}_SUB_FR_SAND_CORRECTED.nc ${file_extpar}_SUB_FR_SILT_CORRECTED.nc ${file_extpar}_SUB_FR_CLAY_CORRECTED.nc ${file_extpar}_FR_OC_CORRECTED.nc ${file_extpar}_SUB_FR_OC_CORRECTED.nc ${file_extpar}_texture.nc
${cdo} merge ${file_bc_soil}_orig.nc ${file_extpar}_texture.nc ${file_bc_soil}_texture.nc
${rm} ${file_extpar}_FR_SAND.nc ${file_extpar}_FR_CLAY.nc ${file_extpar}_FR_SILT.nc ${file_extpar}_SUB_FR_SAND.nc ${file_extpar}_SUB_FR_SILT.nc ${file_extpar}_SUB_FR_CLAY.nc ${file_extpar}_FR_OC.nc ${file_extpar}_SUB_FR_OC.nc ${file_extpar}_FR_OC_MASK.nc ${file_extpar}_SUB_FR_OC_MASK.nc ${file_extpar}_FR_OC_CORRECTED.nc ${file_extpar}_SUB_FR_OC_CORRECTED.nc
${rm} ${file_extpar}_FR_SAND_MASK.nc ${file_extpar}_FR_SILT_MASK.nc ${file_extpar}_FR_CLAY_MASK.nc ${file_extpar}_SUB_FR_SAND_MASK.nc ${file_extpar}_SUB_FR_SILT_MASK.nc ${file_extpar}_SUB_FR_CLAY_MASK.nc ${file_extpar}_SAND_SILT_MASK.nc ${file_extpar}_SUB_SAND_SILT_MASK.nc
${rm} ${file_extpar}_SAND_SILT_CLAY_MASK.nc ${file_extpar}_SUB_SAND_SILT_CLAY_MASK.nc
${rm} ${file_extpar}_FR_????_0.nc ${file_extpar}_SUB_FR_????_0.nc ${file_extpar}_FR_????_CORRECTED.nc ${file_extpar}_SUB_FR_????_CORRECTED.nc
${rm} ${file_extpar}_texture.nc ${file_bc_soil}_soil_depth.nc

# modify root depth
# select extpar field of root depth
${cdo} selvar,ROOTDP ${file_extpar}.nc ${file_extpar}_ROOTDP.nc
# avoid new root depth larger than old soil depth
# select old bc field of soil depth
${cdo} selvar,soil_depth ${file_bc_soil}_texture.nc ${file_bc_soil}_soil_depth.nc
# masks where extpar root depth is either larger or smaller than old soil depth
${cdo} gt ${file_extpar}_ROOTDP.nc ${file_bc_soil}_soil_depth.nc ${file_extpar}_ROOTDP_gt_soil_depth.nc
${cdo} mulc,-1. -addc,-1. ${file_extpar}_ROOTDP_gt_soil_depth.nc ${file_extpar}_ROOTDP_lt_soil_depth.nc
# set extpar root depth to old soil depth, where it is larger than old soil depth
${cdo} add -mul ${file_extpar}_ROOTDP.nc ${file_extpar}_ROOTDP_lt_soil_depth.nc -subc,${delta} -mul ${file_bc_soil}_soil_depth.nc ${file_extpar}_ROOTDP_gt_soil_depth.nc ${file_extpar}_ROOTDP_corrected.nc
# set minimum root depth
${cdo} ltc,${minimum_root_depth} ${file_extpar}_ROOTDP_corrected.nc ${file_extpar}_ROOTDP_minimum.nc
${cdo} add -mul ${file_extpar}_ROOTDP_corrected.nc -mulc,-1. -addc,-1. ${file_extpar}_ROOTDP_minimum.nc -mulc,${minimum_root_depth} ${file_extpar}_ROOTDP_minimum.nc ${file_extpar}_ROOTDP_corrected2.nc
# change variable name
${cdo} chname,ROOTDP,root_depth ${file_extpar}_ROOTDP_corrected2.nc ${file_extpar}_root_depth.nc
# replace old root_depth by extpar root depth in bc_soil file
${cdo} replace ${file_bc_soil}_texture.nc ${file_extpar}_root_depth.nc ${file_bc_soil}_extpar.nc
# clean
${rm} ${file_extpar}_ROOTDP.nc ${file_bc_soil}_texture.nc ${file_bc_soil}_soil_depth.nc ${file_extpar}_ROOTDP_gt_soil_depth.nc ${file_extpar}_ROOTDP_lt_soil_depth.nc ${file_extpar}_ROOTDP_corrected.nc ${file_extpar}_ROOTDP_minimum.nc ${file_extpar}_ROOTDP_corrected2.nc ${file_extpar}_root_depth.nc

# modify bc_phys-file
# roughness length
# select extpar field of roughness length
${cdo} selvar,Z0 ${file_extpar}.nc ${file_extpar}_Z0.nc
# change variable name
${cdo} chname,Z0,roughness_length ${file_extpar}_Z0.nc ${file_extpar}_roughness_length.nc
# replace old roughness length by extpar roughness length in bc_phys file
${cdo} replace ${file_bc_phys}_orig.nc ${file_extpar}_roughness_length.nc ${file_bc_phys}_z0.nc
# clean
${rm} ${file_extpar}_Z0.nc ${file_extpar}_roughness_length.nc

# albedo
# remove very high albedo of mistaken glacier points
${cdo} selvar,ALB ${file_extpar}.nc ${file_extpar}_ALB.nc
${cdo} ltc,${value_albedo_max} ${file_extpar}_ALB.nc ${file_extpar}_ALB_mask_glac.nc
${cdo} mul ${file_extpar}_ALB.nc ${file_extpar}_ALB_mask_glac.nc ${file_extpar}_ALB_glac0.nc
${cdo} setmisstoc,${fillvalue_albedo_glac} -setctomiss,0. ${file_extpar}_ALB_glac0.nc ${file_extpar}_ALB_noglac.nc
# remove very low albedo at the coast
${cdo} gtc,6. ${file_extpar}_ALB_noglac.nc ${file_extpar}_ALB_noglac_mask_coast.nc
${cdo} mul ${file_extpar}_ALB_noglac.nc ${file_extpar}_ALB_noglac_mask_coast.nc ${file_extpar}_ALB_noglac_coast0.nc
${cdo} setmisstoc,${fillvalue_albedo_coast} -setctomiss,0. ${file_extpar}_ALB_noglac_coast0.nc ${file_extpar}_ALB_noglac_nocoast.nc
# average in time (weighted by days per month) and divide by 100 ([% reflection] ==> [albedo])
${cdo} -divc,100. -yearmonmean ${file_extpar}_ALB_noglac_nocoast.nc ${file_extpar}_ALB_noglac_nocoast_mean.nc
# change variable name
${cdo} chname,ALB,albedo ${file_extpar}_ALB_noglac_nocoast_mean.nc ${file_extpar}_albedo_noglac_nocoast_mean.nc
# replace old albedo by extpar albedo in bc_phys file
${cdo} replace ${file_bc_phys}_z0.nc ${file_extpar}_albedo_noglac_nocoast_mean.nc ${file_bc_phys}_z0_albedo.nc
# clean albedo
${rm} ${file_extpar}_ALB.nc ${file_extpar}_ALB_mask_glac.nc ${file_extpar}_ALB_glac0.nc ${file_extpar}_ALB_noglac.nc ${file_extpar}_ALB_noglac_mask_coast.nc ${file_extpar}_ALB_noglac_coast0.nc ${file_extpar}_ALB_noglac_nocoast.nc ${file_extpar}_ALB_noglac_nocoast_mean.nc ${file_extpar}_albedo_noglac_nocoast_mean.nc

# forest fraction
${cdo} selvar,FOR_D ${file_extpar}.nc ${file_extpar}_FOR_D.nc
${cdo} selvar,FOR_E ${file_extpar}.nc ${file_extpar}_FOR_E.nc
${cdo} add ${file_extpar}_FOR_D.nc ${file_extpar}_FOR_E.nc ${file_extpar}_FOR.nc
${cdo} chname,FOR_D,forest_fract ${file_extpar}_FOR.nc ${file_extpar}_forest.nc
# replace old forest fraction by extpar forest fraction in bc_phys file
${cdo} replace ${file_bc_phys}_z0_albedo.nc ${file_extpar}_forest.nc ${file_bc_phys}_z0_albedo_forest.nc
# clean forest
${rm} ${file_extpar}_FOR_D.nc ${file_extpar}_FOR_E.nc ${file_extpar}_FOR.nc

# NDVI
${cdo} selvar,NDVI ${file_extpar}.nc ${file_extpar}_NDVI.nc
# NDVI ==> veg_fract
${cdo} gtc,${lower_bound_NDVI_veg_fract} ${file_extpar}_NDVI.nc ${file_extpar}_NDVI_l_mask.nc
${cdo} setmisstoc,${lower_bound_NDVI_veg_fract} -setctomiss,0. -mul ${file_extpar}_NDVI.nc ${file_extpar}_NDVI_l_mask.nc ${file_extpar}_NDVI_l_bound.nc
${cdo} ltc,${upper_bound_NDVI_veg_fract} ${file_extpar}_NDVI_l_bound.nc ${file_extpar}_NDVI_u_mask.nc
${cdo} setmisstoc,${upper_bound_NDVI_veg_fract} -setctomiss,0. -mul ${file_extpar}_NDVI_l_bound.nc ${file_extpar}_NDVI_u_mask.nc ${file_extpar}_NDVI_bound.nc
${cdo} subc,${lower_bound_NDVI_veg_fract} ${file_extpar}_NDVI_bound.nc ${file_extpar}_NDVI_zero.nc
${cdo} setctomiss,0. -gtc,0. ${file_extpar}_NDVI_zero.nc ${file_extpar}_NDVI_zero_mask.nc
${cdo} setmisstoc,0. -mul ${file_extpar}_NDVI_zero.nc ${file_extpar}_NDVI_zero_mask.nc ${file_extpar}_NDVI_no_neg.nc
${cdo} mulc,${fact_NDVI_veg_fract} ${file_extpar}_NDVI_no_neg.nc ${file_extpar}_NDVI_veg_fract.nc
${cdo} chname,NDVI,veg_fract ${file_extpar}_NDVI_veg_fract.nc ${file_extpar}_veg_fract.nc
# replace old vegetation fraction by extpar vegetation fraction in bc_phys file
${cdo} replace ${file_bc_phys}_z0_albedo_forest.nc ${file_extpar}_veg_fract.nc ${file_bc_phys}_z0_albedo_forest_veg.nc
# clean NDVI ==> veg_fract
${rm} ${file_extpar}_NDVI_l_mask.nc ${file_extpar}_NDVI_l_bound.nc ${file_extpar}_NDVI_u_mask.nc ${file_extpar}_NDVI_bound.nc ${file_extpar}_NDVI_veg_fract.nc ${file_extpar}_NDVI_zero.nc ${file_extpar}_NDVI_zero_mask.nc ${file_extpar}_NDVI_no_neg.nc

# NDVI ==> lai_clim
${cdo} gtc,${lower_bound_NDVI_lai_clim} ${file_extpar}_NDVI.nc ${file_extpar}_NDVI_l_mask.nc
${cdo} setmisstoc,${lower_bound_NDVI_lai_clim} -setctomiss,0. -mul ${file_extpar}_NDVI.nc ${file_extpar}_NDVI_l_mask.nc ${file_extpar}_NDVI_l_bound.nc
${cdo} ltc,${upper_bound_NDVI_lai_clim} ${file_extpar}_NDVI_l_bound.nc ${file_extpar}_NDVI_u_mask.nc
${cdo} setmisstoc,${upper_bound_NDVI_lai_clim} -setctomiss,0. -mul ${file_extpar}_NDVI_l_bound.nc ${file_extpar}_NDVI_u_mask.nc ${file_extpar}_NDVI_bound.nc
${cdo} subc,${lower_bound_NDVI_lai_clim} ${file_extpar}_NDVI_bound.nc ${file_extpar}_NDVI_zero.nc
${cdo} setctomiss,0. -gtc,0. ${file_extpar}_NDVI_zero.nc ${file_extpar}_NDVI_zero_mask.nc
${cdo} setmisstoc,0. -mul ${file_extpar}_NDVI_zero.nc ${file_extpar}_NDVI_zero_mask.nc ${file_extpar}_NDVI_no_neg.nc
${cdo} mulc,${fact_NDVI_lai_clim} ${file_extpar}_NDVI_no_neg.nc ${file_extpar}_NDVI_lai_clim.nc
${cdo} chname,NDVI,lai_clim ${file_extpar}_NDVI_lai_clim.nc ${file_extpar}_lai_clim.nc
# replace old vegetation fraction by extpar vegetation fraction in bc_phys file
${cdo} replace ${file_bc_phys}_z0_albedo_forest_veg.nc ${file_extpar}_lai_clim.nc ${file_bc_phys}_extpar.nc
# clean NDVI ==> lai_clim
${rm} ${file_extpar}_NDVI_l_mask.nc ${file_extpar}_NDVI_l_bound.nc ${file_extpar}_NDVI_u_mask.nc ${file_extpar}_NDVI_bound.nc ${file_extpar}_NDVI_lai_clim.nc ${file_extpar}_NDVI_zero.nc ${file_extpar}_NDVI_zero_mask.nc ${file_extpar}_NDVI_no_neg.nc
# clean NDVI
${rm} ${file_extpar}_NDVI.nc ${file_extpar}_veg_fract.nc ${file_extpar}_lai_clim.nc

# clean bc_phys-files
${rm} ${file_bc_phys}_z0.nc ${file_bc_phys}_z0_albedo.nc ${file_extpar}_forest.nc ${file_bc_phys}_z0_albedo_forest.nc ${file_bc_phys}_z0_albedo_forest_veg.nc

exit
