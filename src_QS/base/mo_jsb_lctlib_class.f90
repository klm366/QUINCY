!> Contains function to read the landcover library file
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!>
!! This module provides the data from the landcover library file, which contains additional information on the
!! landcover data used in a particular run of JSBACH.
!! The name of the landcover library file is specified in the JSBACH configuration file (run.def) under the
!! keyword "LCT_FILE".
!!
MODULE mo_jsb_lctlib_class
  !
  ! This module provides the data from the landcover library file, which contains additional information on the
  ! landcover data used in a particular run of JSBACH.
  ! The name of the landcover library file is specified in the JSBACH configuration file (run.def) under the
  ! keyword "LCT_FILE".
  !
  ! Authors: Reiner Schnur, 7/2003
  !          Christian H. Reick, 8/2003
  !
!  USE mo_jsb_parallel, ONLY: my_process_is_stdio, my_process_is_mpi_parallel, p_io, p_bcast
  USE mo_kind,         ONLY: wp
  USE mo_exception,    ONLY: message_text, finish
  USE mo_util,         ONLY: int2string

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_lctlib_element, Read_lctlib

  TYPE t_lctlib_element
     ! --- general parameters ------------------------------------------------------------------------------------------------------
    INTEGER  :: LctNumber            !! Unique number for each landcover type (1, 2, 3, ..)
    CHARACTER(LEN=16)  :: LctName    !! Names of landcover types
    INTEGER  :: LandcoverClass       !! Landcover classes (not to be confused with landcover types!):
                                                 !! .. 0: Bare soil; 1:Glacier; 2: Lake;
                                                 !! .. 3: Natural forest; 4: Natural Grassland;
                                                 !! .. 5: Other natural vegetation; 6: crops; 7: pastures
    LOGICAL  :: NaturalVegFlag       !! "true" for natural vegetation landcover types, i.e. non-agricultural vegetation
    LOGICAL  :: ForestFlag           !! "true" for forest landcover types
    LOGICAL  :: GrassFlag            !! "True" for natural graslands
    LOGICAL  :: CropFlag             !! "true" for croplands
    LOGICAL  :: PastureFlag          !! "true" for landcover of type "pasture"
    LOGICAL  :: LakeFlag             !! "true" for a land surface class that signifies lake or ocean
    LOGICAL  :: GlacierFlag          !! "true" for a land surface class that signifies glaciers
    LOGICAL  :: BareSoilFlag         !! "true" if none of the foregoing landcover classes
    LOGICAL  :: C4flag               !! derived from ps_pathway, "true" if: ps_pathway == ic4phot

    INTEGER  :: growthform 
    INTEGER  :: ps_pathway
    INTEGER  :: phenology_type
    REAL(wp) :: lai_max
    REAL(wp) :: vegetation_height
    REAL(wp) :: sla
    REAL(wp) :: sigma_vis 
    REAL(wp) :: sigma_nir 
    REAL(wp) :: omega_clumping 
    REAL(wp) :: crown_shape_factor 
    REAL(wp) :: cn_leaf             
    REAL(wp) :: cn_leaf_min
    REAL(wp) :: cn_leaf_max
    REAL(wp) :: np_leaf  
    REAL(wp) :: np_leaf_min
    REAL(wp) :: np_leaf_max
    REAL(wp) :: k0_fn_struc
    REAL(wp) :: fn_oth_min
    REAL(wp) :: t_jmax_opt
    REAL(wp) :: t_jmax_omega
    REAL(wp) :: g0
    REAL(wp) :: g1_medlyn
    REAL(wp) :: g1_bberry
    REAL(wp) :: gmin

    ! turnover times
    REAL(wp) :: tau_leaf          , &
                tau_fine_root     , &
                tau_coarse_root   , &
                tau_branch        , &
                tau_sap_wood      , &
                tau_fruit         , &
                tau_seed_litter   , &
                tau_seed_est      , &
                tau_mycorrhiza

    ! N uptake parameters
    REAL(wp) :: vmax_uptake_n     , &
                vmax_uptake_p     , &
                bnf_base          , &
                vmax_symb_bnf     , &
                km_symb_bnf

    ! Vegetation dynamics parameters
    REAL(wp) :: lambda_est_light  , &
                k_est_light       , &
                seed_size         , &
                k1_mort_greff

    ! Phenology parameters
    REAL(wp) :: beta_soil_flush       , &
                beta_soil_senescence  , &
                gdd_req_max           , &
                k_gdd_dormance        , &
                t_air_senescence      , &
                min_leaf_age

    ! Allocation paramters
    REAL(wp) :: frac_sapwood_branch   , &
                wood_density          , &
                k_latosa              , &
                k_crtos               , &
                k_rtos                , &
                k2_fruit_alloc        , &
                allom_k1              , &
                allom_k2              , &
                phi_leaf_min          , &
                k_root                , &
                k_sapwood             , &
                c0_allom              , &
                fstore_target

    ! Soil
    REAL(wp) :: &
      & k_root_dist                   , &
      & k_som_fast_init               , &
      & k_som_slow_init

    ! copy of jsbach param
    REAL(wp) :: &
      & AlbedoLitterVIS               , &
      & AlbedoLitterNIR               , &
      & VegRoughness                  , &
      & MinVegRoughness               , &
      & MaxVegRoughness
    ! --- END QUINCY model parameters (MPI-BGC Jena) -----------------------------------------------------------------------------

  END TYPE t_lctlib_element

CONTAINS

  FUNCTION Read_lctlib(lctlib_file_name) RESULT(return_value)
    !
    ! Get landcover library, i.e. lookup table for landcover types
    !
    ! The structure of the lookup table in the input file is as follows. It can contain only three types of lines in arbitrary
    ! order:
    !
    !     Comment lines: These contain before a '#' only blanks
    !       Blank lines: These contain only blanks
    !        Data lines: The first nonblank characters form a keyword which is follwed by the data (separated by blanks)
    !
    ! Keywords are the same as the names of the lctlib components, with one exception:
    !
    !        NLCT : This is the keyword after which the number of landcover types, i.e. the number of data columns
    !               in the file is expected. This number has to be identical with the number of landcover types
    !               from the landcover data file (keyword: LCTLIB_FILE in the JSBACH configuration file).
    !
    ! After all other keywords NLCT columns of data are expected.
    !
    ! Example:
    !            # ---- LANDCOVER LIBRARY FILE -------
    !            NLCT 3
    !            LctNumber               3 5 9
    !            # Phenology types: 0= none 1=summergreen, 2=evergreen, 3=raingreen, 4=grasses
    !            PhenologyType          2 2 4
    !            # C4flag: 0=C3, 1=C4
    !            C4flag                 0 1 1
    !
    ! The first string on each line (except first line and comment lines) must correspond to the name of the lctlib component
    ! (case doesn't matter)
    ! For lctlib components of type LOGICAL use 0/1 in the file to indicate .FALSE./.TRUE.
    ! For lctlib components of type INTEGER the numbers on the line must be integer values
    ! For lctlib components of type REAL the numbers on the line can be either integer or floating point
    ! Comments can appear on any line, everything to the right of and including "#" is disregarded
    !
    ! The file can contain more keywords than needed --- therefore the same file can be used by several model components.
    !
#ifdef __QUINCY_STANDALONE__
    USE mo_util,        ONLY: tolower
#else
    USE mo_util_string, ONLY: tolower
#endif
    ! USE statements necessary for modifying QUINCY model input values of lctlib parameters
    USE mo_jsb_math_constants,      ONLY: one_day, one_year
    USE mo_jsb_physical_constants,  ONLY: molar_mass_C, molar_mass_N, molar_mass_P, Dwv, Dco2, Tzero
    USE mo_veg_constants,           ONLY: carbon_per_dryweight_leaf, sm2lm_grass, igrass, itree, icrop
    USE mo_q_assimi_constants,      ONLY: ic4phot

    CHARACTER(len=*), INTENT(in) :: lctlib_file_name
    TYPE(t_lctlib_element), POINTER    :: return_value(:)

    TYPE(t_lctlib_element), ALLOCATABLE :: lctlib(:)

    INTEGER, PARAMETER            :: lctlib_file_unit = 66

    CHARACTER(len=30)  :: key
    CHARACTER(len=256) :: line
    INTEGER            :: pos_comment, read_status
    INTEGER            :: pos,length
    CHARACTER(len=2)   :: blank_set = " "//achar(9) !< the blank characters: BLANK and TAB
    INTEGER,ALLOCATABLE:: itmp(:)                   !< temporary array used for input of logicals
    INTEGER            :: apu

    LOGICAL            :: exist_LctNumber              = .FALSE.
    LOGICAL            :: exist_LctName                = .FALSE.
    LOGICAL            :: exist_LandcoverClass         = .FALSE.
    LOGICAL            :: exist_AlbedoSnowVisMin       = .FALSE.
    LOGICAL            :: exist_AlbedoSnowVisMax       = .FALSE.
    LOGICAL            :: exist_AlbedoSnowNirMin       = .FALSE.
    LOGICAL            :: exist_AlbedoSnowNirMax       = .FALSE.
    LOGICAL            :: exist_AlbedoSnowMin          = .FALSE.
    LOGICAL            :: exist_AlbedoSnowMax          = .FALSE.
    LOGICAL            :: exist_AlbedoCanopyVIS        = .FALSE.
    LOGICAL            :: exist_AlbedoCanopyNIR        = .FALSE.
    LOGICAL            :: exist_AlbedoLitterVis        = .FALSE.
    LOGICAL            :: exist_AlbedoLitterNir        = .FALSE.
    LOGICAL            :: exist_nitrogenScalingFlag    = .FALSE.
    LOGICAL            :: exist_C4flag                 = .FALSE.
    LOGICAL            :: exist_CarboxRate             = .FALSE.
    LOGICAL            :: exist_ETransport             = .FALSE.
    LOGICAL            :: exist_VegHeight              = .FALSE.
    LOGICAL            :: exist_VegRoughness           = .FALSE.
    LOGICAL            :: exist_MinVegRoughness        = .FALSE.
    LOGICAL            :: exist_MaxVegRoughness        = .FALSE.
    LOGICAL            :: exist_PhenologyType          = .FALSE.
    LOGICAL            :: exist_CanopyResistanceMin    = .FALSE.
    LOGICAL            :: exist_MaxLAI                 = .FALSE.
    LOGICAL            :: exist_StemArea               = .FALSE.
    LOGICAL            :: exist_specificLeafArea_C     = .FALSE.
    LOGICAL            :: exist_alpha_nr_ind           = .FALSE.
    LOGICAL            :: exist_beta_nr_ind            = .FALSE.
    LOGICAL            :: exist_alpha_leaf             = .FALSE.
    LOGICAL            :: exist_beta_leaf              = .FALSE.
    LOGICAL            :: exist_knorr_Tau_w            = .FALSE.
    LOGICAL            :: exist_knorr_T_phi            = .FALSE.
    LOGICAL            :: exist_knorr_T_r              = .FALSE.
    LOGICAL            :: exist_knorr_Day_c            = .FALSE.
    LOGICAL            :: exist_knorr_Day_r            = .FALSE.
    LOGICAL            :: exist_knorr_k_l              = .FALSE.
    LOGICAL            :: exist_knorr_leaf_growth_rate = .FALSE.
    LOGICAL            :: exist_knorr_max_lai          = .FALSE.
    LOGICAL            :: exist_reserveC2leafC         = .FALSE.
    LOGICAL            :: exist_fract_npp_2_woodPool    = .FALSE.
    LOGICAL            :: exist_fract_npp_2_reservePool = .FALSE.
    LOGICAL            :: exist_fract_npp_2_exudates    = .FALSE.
    LOGICAL            :: exist_fract_green_2_herbivory = .FALSE.
    LOGICAL            :: exist_tau_c_woods        = .FALSE.
    LOGICAL            :: exist_LAI_shed_constant      = .FALSE.
    LOGICAL            :: exist_Max_C_content_woods       = .FALSE.
    LOGICAL            :: exist_ClumpinessFactor          = .FALSE.
    LOGICAL            :: exists_dynamic_PFT              = .FALSE.
    LOGICAL            :: exists_woody_PFT                = .FALSE.
    LOGICAL            :: exists_pasture_PFT              = .FALSE.
    LOGICAL            :: exists_bclimit_min_cold_mmtemp  = .FALSE.
    LOGICAL            :: exists_bclimit_max_cold_mmtemp  = .FALSE.
    LOGICAL            :: exists_bclimit_max_warm_mmtemp  = .FALSE.
    LOGICAL            :: exists_bclimit_min_temprange    = .FALSE.
    LOGICAL            :: exists_bclimit_min_gdd          = .FALSE.
    LOGICAL            :: exists_gdd_base                 = .FALSE.
    LOGICAL            :: exists_upper_tlim               = .FALSE.
    LOGICAL            :: exists_fract_wood_2_onSite       = .FALSE.
    LOGICAL            :: exists_fract_wood_2_paper        = .FALSE.
    LOGICAL            :: exists_fract_wood_2_construction = .FALSE.
    LOGICAL            :: exists_moist_extinction      = .FALSE.
    LOGICAL            :: exists_fuel_dens             = .FALSE.
    LOGICAL            :: exists_flame_length_f        = .FALSE.
    LOGICAL            :: exists_crown_length          = .FALSE.
    LOGICAL            :: exists_bark_par1             = .FALSE.
    LOGICAL            :: exists_bark_par2             = .FALSE.
    LOGICAL            :: exists_RCK                   = .FALSE.
    LOGICAL            :: exists_mort_prob             = .FALSE.
    LOGICAL            :: exists_woodlittersize_coef      = .FALSE.
    LOGICAL            :: exists_woodlit_coef             = .FALSE.
    LOGICAL            :: exists_leaflit_coef             = .FALSE.
    LOGICAL            :: exists_litveg_coef              = .FALSE.
    ! quincy model - start
    LOGICAL            :: exists_growthform = .FALSE.
    LOGICAL            :: exists_ps_pathway = .FALSE.
    LOGICAL            :: exists_phenology_type = .FALSE.
    LOGICAL            :: exists_lai_max = .FALSE.
    LOGICAL            :: exists_vegetation_height = .FALSE.
    LOGICAL            :: exists_sla = .FALSE.
    LOGICAL            :: exists_sigma_vis = .FALSE.
    LOGICAL            :: exists_sigma_nir = .FALSE.
    LOGICAL            :: exists_omega_clumping = .FALSE.
    LOGICAL            :: exists_crown_shape_factor = .FALSE.
    LOGICAL            :: exists_cn_leaf = .FALSE.
    LOGICAL            :: exists_cn_leaf_min = .FALSE.
    LOGICAL            :: exists_cn_leaf_max = .FALSE.
    LOGICAL            :: exists_np_leaf = .FALSE.
    LOGICAL            :: exists_np_leaf_min = .FALSE.
    LOGICAL            :: exists_np_leaf_max = .FALSE.
    LOGICAL            :: exists_k0_fn_struc = .FALSE.
    LOGICAL            :: exists_fn_oth_min = .FALSE.
    LOGICAL            :: exists_t_jmax_opt = .FALSE.
    LOGICAL            :: exists_t_jmax_omega = .FALSE.
    LOGICAL            :: exists_g0 = .FALSE.
    LOGICAL            :: exists_g1_medlyn = .FALSE.
    LOGICAL            :: exists_g1_bberry = .FALSE.
    LOGICAL            :: exists_gmin = .FALSE.
    LOGICAL            :: exists_tau_leaf = .FALSE.
    LOGICAL            :: exists_tau_fine_root = .FALSE.
    LOGICAL            :: exists_tau_coarse_root = .FALSE.
    LOGICAL            :: exists_tau_branch = .FALSE.
    LOGICAL            :: exists_tau_sap_wood = .FALSE.
    LOGICAL            :: exists_tau_fruit = .FALSE.
    LOGICAL            :: exists_tau_seed_litter = .FALSE.
    LOGICAL            :: exists_tau_seed_est = .FALSE.
    LOGICAL            :: exists_tau_mycorrhiza = .FALSE.
    LOGICAL            :: exists_vmax_uptake_n = .FALSE.
    LOGICAL            :: exists_vmax_uptake_p = .FALSE.
    LOGICAL            :: exists_bnf_base = .FALSE.
    LOGICAL            :: exists_vmax_symb_bnf = .FALSE.
    LOGICAL            :: exists_km_symb_bnf = .FALSE.
    LOGICAL            :: exists_lambda_est_light = .FALSE.
    LOGICAL            :: exists_k_est_light = .FALSE.
    LOGICAL            :: exists_seed_size = .FALSE.
    LOGICAL            :: exists_k1_mort_greff = .FALSE.
    LOGICAL            :: exists_beta_soil_flush = .FALSE.
    LOGICAL            :: exists_beta_soil_senescence = .FALSE.
    LOGICAL            :: exists_gdd_req_max = .FALSE.
    LOGICAL            :: exists_k_gdd_dormance = .FALSE.
    LOGICAL            :: exists_t_air_senescence = .FALSE.
    LOGICAL            :: exists_min_leaf_age = .FALSE.
    LOGICAL            :: exists_frac_sapwood_branch = .FALSE.
    LOGICAL            :: exists_wood_density = .FALSE.
    LOGICAL            :: exists_k_latosa = .FALSE.
    LOGICAL            :: exists_k_crtos = .FALSE.
    LOGICAL            :: exists_k_rtos = .FALSE.
    LOGICAL            :: exists_k2_fruit_alloc = .FALSE.
    LOGICAL            :: exists_allom_k1 = .FALSE.
    LOGICAL            :: exists_allom_k2 = .FALSE.
    LOGICAL            :: exists_phi_leaf_min = .FALSE.
    LOGICAL            :: exists_k_root = .FALSE.
    LOGICAL            :: exists_k_sapwood = .FALSE.
    LOGICAL            :: exists_c0_allom = .FALSE.
    LOGICAL            :: exists_fstore_target = .FALSE.
    LOGICAL            :: exists_k_root_dist = .FALSE.
    LOGICAL            :: exists_k_som_fast_init = .FALSE.
    LOGICAL            :: exists_k_som_slow_init = .FALSE.
    LOGICAL            :: exists_AlbedoLitterVIS = .FALSE.
    LOGICAL            :: exists_AlbedoLitterNIR = .FALSE.
    LOGICAL            :: exists_VegRoughness    = .FALSE.
    LOGICAL            :: exists_MinVegRoughness = .FALSE.
    LOGICAL            :: exists_MaxVegRoughness = .FALSE.
    ! quincy model - end
    INTEGER            :: i, j
    INTEGER            :: nlct, npft

    !>Read lctlib data and eventually allocate also memory for lctlib data on io-processor
    !>
!    IF (my_process_is_stdio()) THEN

       OPEN(unit=lctlib_file_unit, file=lctlib_file_name, form='FORMATTED', status='OLD', iostat=read_status)
       IF (read_status /= 0) CALL finish('init_lctlib','Error opening landcover library file')

       ! --- find keyword NLCT in landcover library file

       DO
          READ(lctlib_file_unit,'(A256)', IOSTAT=read_status) line
          IF (read_status /= 0) THEN
             CALL finish('init_lctlib','No keyword NLCT found in land cover library file '//TRIM(lctlib_file_name))
          END IF
          ! Look for comment
          pos_comment = SCAN(line,"#")                 ! Start position of comment, 0 if none present
          IF (pos_comment > 1) THEN
             line = TRIM(ADJUSTL(line(1:pos_comment))) ! Disregard everything to the right of "#", ..
             ! .. adjust to the left and trim to the right
          ELSE IF (pos_comment == 1) THEN
             CYCLE                                     ! Disregard whole line
          ELSE
             line = TRIM(ADJUSTL(line))                ! Only disregard blanks to the left and right
          ENDIF
          length=LEN_TRIM(line)
          IF(length== 0) CYCLE                 ! Line is empty

          pos = SCAN(line,blank_set)                   ! Position of first blank character
          IF(pos == 0)   CALL finish('init_lctlib',"Wrong syntax in lctlib definitions")

          READ(line(1:pos-1),'(A)') key

          IF(tolower(TRIM(key)) == 'nlct') THEN
             READ(line(pos:length),*,IOSTAT=read_status) nlct
             IF (read_status /= 0) THEN
                CALL finish('init_lctlib','Could not read number of landcover types (keyword: NLCT) from '//TRIM(lctlib_file_name))
             END IF
             EXIT ! found number of landcover types in landcover library file --- continue after loop
          END IF
       END DO

!    END IF

!    IF (my_process_is_mpi_parallel()) CALL p_bcast(nlct, p_io)

    ALLOCATE(lctlib(nlct))

!    IF (my_process_is_stdio()) THEN

       REWIND(unit=lctlib_file_unit) ! Go back to beginning of landcover library file

       ALLOCATE(itmp(nlct))

       DO
          READ(lctlib_file_unit,'(A256)', IOSTAT=read_status) line
          IF (read_status /= 0) EXIT                   ! Finished reading

          ! Look for comment
          pos_comment = SCAN(line,"#")                 ! Start position of comment, 0 if none present
          IF (pos_comment > 1) THEN
             line = TRIM(ADJUSTL(line(1:pos_comment))) ! Disregard everything to the right of "#", ..
             ! .. adjust to the left and trim to the right
          ELSE IF (pos_comment == 1) THEN
             CYCLE                                     ! Disregard whole line
          ELSE
             line = TRIM(ADJUSTL(line))                ! Only disregard blanks to the left and right
          ENDIF
          length=LEN_TRIM(line)
          IF(length== 0) CYCLE                 ! Line is empty

          pos = SCAN(line,blank_set)                   ! Position of first blank character
          IF(pos == 0)   CALL finish('init_lctlib',"Wrong syntax in lctlib definitions")

          READ(line(1:pos-1),'(A)') key
          key = tolower(key)

          IF(TRIM(key) == "nlct") CYCLE ! nlct already read above

          ! read the jsbach4 lctLib parameter
          SELECT CASE (TRIM(key))
          CASE ('lctnumber')
            READ(line(pos:length),*) lctlib(1:nlct)%LctNumber
            exist_LctNumber = .TRUE.
          CASE ('lctname')                                           ! Name of landcover type
            READ(line(pos:length),*) lctlib(1:nlct)%LctName
            exist_LctName = .TRUE.
          CASE ('landcoverclass')                                    ! Landcover class
            READ(line(pos:length),*) lctlib(1:nlct)%LandcoverClass
            exist_LandcoverClass = .TRUE.

          ! ------------------------------------------------------------------
          ! read the quincy lctLib parameter
          !IF (model_scheme_char == "jsbach_quincy" .OR. model_scheme_char == "quincy") THEN
          !  SELECT CASE (TRIM(key))
            CASE ('growthform')
              READ(line(pos:length),*) lctlib(1:nlct)%growthform
              exists_growthform = .TRUE.
            CASE ('ps_pathway')
              READ(line(pos:length),*) lctlib(1:nlct)%ps_pathway
              exists_ps_pathway = .TRUE.
            CASE ('phenology_type')
              READ(line(pos:length),*) lctlib(1:nlct)%phenology_type
              exists_phenology_type = .TRUE.
            CASE ('lai_max')
              READ(line(pos:length),*) lctlib(1:nlct)%lai_max
              exists_lai_max = .TRUE.
            CASE ('vegetation_height')
              READ(line(pos:length),*) lctlib(1:nlct)%vegetation_height
              exists_vegetation_height = .TRUE.
            CASE ('sla')
              READ(line(pos:length),*) lctlib(1:nlct)%sla
              exists_sla = .TRUE.
            CASE ('sigma_vis')
              READ(line(pos:length),*) lctlib(1:nlct)%sigma_vis
              exists_sigma_vis = .TRUE.
            CASE ('sigma_nir')
              READ(line(pos:length),*) lctlib(1:nlct)%sigma_nir
              exists_sigma_nir = .TRUE.
            CASE ('omega_clumping')
              READ(line(pos:length),*) lctlib(1:nlct)%omega_clumping
              exists_omega_clumping = .TRUE.
            CASE ('crown_shape_factor')
              READ(line(pos:length),*) lctlib(1:nlct)%crown_shape_factor
              exists_crown_shape_factor = .TRUE.
            CASE ('cn_leaf')
              READ(line(pos:length),*) lctlib(1:nlct)%cn_leaf
              exists_cn_leaf = .TRUE.
            CASE ('cn_leaf_min')
              READ(line(pos:length),*) lctlib(1:nlct)%cn_leaf_min
              exists_cn_leaf_min = .TRUE.
            CASE ('cn_leaf_max')
              READ(line(pos:length),*) lctlib(1:nlct)%cn_leaf_max
              exists_cn_leaf_max = .TRUE.
            CASE ('np_leaf')
              READ(line(pos:length),*) lctlib(1:nlct)%np_leaf
              exists_np_leaf = .TRUE.
            CASE ('np_leaf_min')
              READ(line(pos:length),*) lctlib(1:nlct)%np_leaf_min
              exists_np_leaf_min = .TRUE.
            CASE ('np_leaf_max')
              READ(line(pos:length),*) lctlib(1:nlct)%np_leaf_max
              exists_np_leaf_max = .TRUE.
            CASE ('k0_fn_struc')
              READ(line(pos:length),*) lctlib(1:nlct)%k0_fn_struc
              exists_k0_fn_struc = .TRUE.
            CASE ('fn_oth_min')
              READ(line(pos:length),*) lctlib(1:nlct)%fn_oth_min
              exists_fn_oth_min = .TRUE.
            CASE ('t_jmax_opt')
              READ(line(pos:length),*) lctlib(1:nlct)%t_jmax_opt
              exists_t_jmax_opt = .TRUE.
            CASE ('t_jmax_omega')
              READ(line(pos:length),*) lctlib(1:nlct)%t_jmax_omega
              exists_t_jmax_omega = .TRUE.
            CASE ('g0')
              READ(line(pos:length),*) lctlib(1:nlct)%g0
              exists_g0 = .TRUE.
            CASE ('g1_medlyn')
              READ(line(pos:length),*) lctlib(1:nlct)%g1_medlyn
              exists_g1_medlyn = .TRUE.
            CASE ('g1_bberry')
              READ(line(pos:length),*) lctlib(1:nlct)%g1_bberry
              exists_g1_bberry = .TRUE.
            CASE ('gmin')
              READ(line(pos:length),*) lctlib(1:nlct)%gmin
              exists_gmin = .TRUE.

            CASE ('tau_leaf')
              READ(line(pos:length),*) lctlib(1:nlct)%tau_leaf
              exists_tau_leaf = .TRUE.
            CASE ('tau_fine_root')
              READ(line(pos:length),*) lctlib(1:nlct)%tau_fine_root
              exists_tau_fine_root = .TRUE.
            CASE ('tau_coarse_root')
              READ(line(pos:length),*) lctlib(1:nlct)%tau_coarse_root
              exists_tau_coarse_root = .TRUE.
            CASE ('tau_branch')
              READ(line(pos:length),*) lctlib(1:nlct)%tau_branch
              exists_tau_branch = .TRUE.
            CASE ('tau_sap_wood')
              READ(line(pos:length),*) lctlib(1:nlct)%tau_sap_wood
              exists_tau_sap_wood = .TRUE.
            CASE ('tau_fruit')
              READ(line(pos:length),*) lctlib(1:nlct)%tau_fruit
              exists_tau_fruit = .TRUE.
            CASE ('tau_seed_litter')
              READ(line(pos:length),*) lctlib(1:nlct)%tau_seed_litter
              exists_tau_seed_litter = .TRUE.
            CASE ('tau_seed_est')
              READ(line(pos:length),*) lctlib(1:nlct)%tau_seed_est
              exists_tau_seed_est = .TRUE.
            CASE ('tau_mycorrhiza')
              READ(line(pos:length),*) lctlib(1:nlct)%tau_mycorrhiza
              exists_tau_mycorrhiza = .TRUE.

            CASE ('vmax_uptake_n')
              READ(line(pos:length),*) lctlib(1:nlct)%vmax_uptake_n
              exists_vmax_uptake_n = .TRUE.
            CASE ('vmax_uptake_p')
              READ(line(pos:length),*) lctlib(1:nlct)%vmax_uptake_p
              exists_vmax_uptake_p = .TRUE.
            CASE ('bnf_base')
              READ(line(pos:length),*) lctlib(1:nlct)%bnf_base
              exists_bnf_base = .TRUE.
            CASE ('vmax_symb_bnf')
              READ(line(pos:length),*) lctlib(1:nlct)%vmax_symb_bnf
              exists_vmax_symb_bnf = .TRUE.
            CASE ('km_symb_bnf')
              READ(line(pos:length),*) lctlib(1:nlct)%km_symb_bnf
              exists_km_symb_bnf = .TRUE.

            CASE ('lambda_est_light')
              READ(line(pos:length),*) lctlib(1:nlct)%lambda_est_light
              exists_lambda_est_light = .TRUE.
            CASE ('k_est_light')
              READ(line(pos:length),*) lctlib(1:nlct)%k_est_light
              exists_k_est_light = .TRUE.
            CASE ('seed_size')
              READ(line(pos:length),*) lctlib(1:nlct)%seed_size
              exists_seed_size = .TRUE.
            CASE ('k1_mort_greff')
              READ(line(pos:length),*) lctlib(1:nlct)%k1_mort_greff
              exists_k1_mort_greff = .TRUE.

            CASE ('beta_soil_flush')
              READ(line(pos:length),*) lctlib(1:nlct)%beta_soil_flush
              exists_beta_soil_flush = .TRUE.
            CASE ('beta_soil_senescence')
              READ(line(pos:length),*) lctlib(1:nlct)%beta_soil_senescence
              exists_beta_soil_senescence = .TRUE.
            CASE ('gdd_req_max')
              READ(line(pos:length),*) lctlib(1:nlct)%gdd_req_max
              exists_gdd_req_max = .TRUE.
            CASE ('k_gdd_dormance')
              READ(line(pos:length),*) lctlib(1:nlct)%k_gdd_dormance
              exists_k_gdd_dormance = .TRUE.
            CASE ('t_air_senescence')
              READ(line(pos:length),*) lctlib(1:nlct)%t_air_senescence
              exists_t_air_senescence = .TRUE.
            CASE ('min_leaf_age')
              READ(line(pos:length),*) lctlib(1:nlct)%min_leaf_age
              exists_min_leaf_age = .TRUE.

            CASE ('frac_sapwood_branch')
              READ(line(pos:length),*) lctlib(1:nlct)%frac_sapwood_branch
              exists_frac_sapwood_branch = .TRUE.
            CASE ('wood_density')
              READ(line(pos:length),*) lctlib(1:nlct)%wood_density
              exists_wood_density = .TRUE.
            CASE ('k_latosa')
              READ(line(pos:length),*) lctlib(1:nlct)%k_latosa
              exists_k_latosa = .TRUE.
            CASE ('k_crtos')
              READ(line(pos:length),*) lctlib(1:nlct)%k_crtos
              exists_k_crtos = .TRUE.
            CASE ('k_rtos')
              READ(line(pos:length),*) lctlib(1:nlct)%k_rtos
              exists_k_rtos = .TRUE.
            CASE ('k2_fruit_alloc')
              READ(line(pos:length),*) lctlib(1:nlct)%k2_fruit_alloc
              exists_k2_fruit_alloc = .TRUE.
            CASE ('allom_k1')
              READ(line(pos:length),*) lctlib(1:nlct)%allom_k1
              exists_allom_k1 = .TRUE.
            CASE ('allom_k2')
              READ(line(pos:length),*) lctlib(1:nlct)%allom_k2
              exists_allom_k2 = .TRUE.
            CASE ('phi_leaf_min')
              READ(line(pos:length),*) lctlib(1:nlct)%phi_leaf_min
              exists_phi_leaf_min = .TRUE.
            CASE ('k_root')
              READ(line(pos:length),*) lctlib(1:nlct)%k_root
              exists_k_root = .TRUE.
            CASE ('k_sapwood')
              READ(line(pos:length),*) lctlib(1:nlct)%k_sapwood
              exists_k_sapwood = .TRUE.
            CASE ('c0_allom')
              READ(line(pos:length),*) lctlib(1:nlct)%c0_allom
              exists_c0_allom = .TRUE.
            CASE ('fstore_target')
              READ(line(pos:length),*) lctlib(1:nlct)%fstore_target
              exists_fstore_target = .TRUE.

            CASE ('k_root_dist')
              READ(line(pos:length),*) lctlib(1:nlct)%k_root_dist
              exists_k_root_dist = .TRUE.
            CASE ('k_som_fast_init')
              READ(line(pos:length),*) lctlib(1:nlct)%k_som_fast_init
              exists_k_som_fast_init = .TRUE.
            CASE ('k_som_slow_init')
              READ(line(pos:length),*) lctlib(1:nlct)%k_som_slow_init
              exists_k_som_slow_init = .TRUE.

            CASE ('albedolittervis')
              READ(line(pos:length),*) lctlib(1:nlct)%AlbedoLitterVIS
              exists_AlbedoLitterVIS = .TRUE.
            CASE ('albedolitternir')
              READ(line(pos:length),*) lctlib(1:nlct)%AlbedoLitterNIR
              exists_AlbedoLitterNIR = .TRUE.
            CASE ('vegroughness')
              READ(line(pos:length),*) lctlib(1:nlct)%VegRoughness
              exists_VegRoughness = .TRUE.
            CASE ('minvegroughness')
              READ(line(pos:length),*) lctlib(1:nlct)%MinVegRoughness
              exists_MinVegRoughness = .TRUE.
            CASE ('maxvegroughness')
              READ(line(pos:length),*) lctlib(1:nlct)%MaxVegRoughness
              exists_MaxVegRoughness = .TRUE.
            CASE default
              ! nothing to do
            END SELECT
          !END IF
          ! quincy end
          ! ------------------------------------------------------------------
       END DO
       DEALLOCATE(itmp)

       !> test if the lctLib parameter could be read, if not finish()
       !>

       ! three main parameters
       IF(.NOT. exist_LctNumber) &
            CALL finish('init_lctlib','No data for LctNumber found in '//TRIM(lctlib_file_name))
       IF(.NOT. exist_LctName) &
!            CALL finish('init_lctlib','No data for LctName found in '//TRIM(lctlib_file_name))
            ! Set default
            lctlib%LctName = "undefined"
       IF(.NOT. exist_LandcoverClass) &
            CALL finish('init_lctlib','No data for LandcoverClass found in '//TRIM(lctlib_file_name))

       ! model scheme: jsbach_quincy / quincy
       !IF (model_scheme_char == "jsbach_quincy" .OR. model_scheme_char == "quincy") THEN
        IF(.NOT. exists_growthform) &
          CALL finish('init_lctlib','No data for growthform found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_ps_pathway) &
          CALL finish('init_lctlib','No data for ps_pathway found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_phenology_type) &
          CALL finish('init_lctlib','No data for phenology_type found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_lai_max) &
          CALL finish('init_lctlib','No data for lai_max found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_vegetation_height) &
          CALL finish('init_lctlib','No data for vegetation_height found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_sla) &
          CALL finish('init_lctlib','No data for sla found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_sigma_vis) &
          CALL finish('init_lctlib','No data for sigma_vis found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_sigma_nir) &
          CALL finish('init_lctlib','No data for sigma_nir found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_omega_clumping) &
          CALL finish('init_lctlib','No data for omega_clumping found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_crown_shape_factor) &
          CALL finish('init_lctlib','No data for crown_shape_factor found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_cn_leaf) &
          CALL finish('init_lctlib','No data for cn_leaf found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_cn_leaf_min) &
          CALL finish('init_lctlib','No data for cn_leaf_min found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_cn_leaf_max) &
          CALL finish('init_lctlib','No data for cn_leaf_max found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_np_leaf) &
          CALL finish('init_lctlib','No data for np_leaf found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_np_leaf_min) &
          CALL finish('init_lctlib','No data for np_leaf_min found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_np_leaf_max) &
          CALL finish('init_lctlib','No data for np_leaf_max found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_k0_fn_struc) &
          CALL finish('init_lctlib','No data for k0_fn_struc found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_fn_oth_min) &
          CALL finish('init_lctlib','No data for fn_oth_min found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_t_jmax_opt) &
          CALL finish('init_lctlib','No data for t_jmax_opt found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_t_jmax_omega) &
          CALL finish('init_lctlib','No data for t_jmax_omega found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_g0) &
          CALL finish('init_lctlib','No data for g0 found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_g1_medlyn) &
          CALL finish('init_lctlib','No data for g1_medlyn found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_g1_bberry) &
          CALL finish('init_lctlib','No data for g1_bberry found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_gmin) &
          CALL finish('init_lctlib','No data for gmin found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_tau_leaf) &
          CALL finish('init_lctlib','No data for tau_leaf found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_tau_fine_root) &
          CALL finish('init_lctlib','No data for tau_fine_root found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_tau_coarse_root) &
          CALL finish('init_lctlib','No data for tau_coarse_root found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_tau_branch) &
          CALL finish('init_lctlib','No data for tau_branch found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_tau_sap_wood) &
          CALL finish('init_lctlib','No data for tau_sap_wood found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_tau_fruit) &
          CALL finish('init_lctlib','No data for tau_fruit found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_tau_seed_litter) &
          CALL finish('init_lctlib','No data for tau_seed_litter found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_tau_seed_est) &
          CALL finish('init_lctlib','No data for tau_seed_est found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_tau_mycorrhiza) &
          CALL finish('init_lctlib','No data for tau_mycorrhiza found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_vmax_uptake_n) &
          CALL finish('init_lctlib','No data for vmax_uptake_n found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_vmax_uptake_p) &
          CALL finish('init_lctlib','No data for vmax_uptake_p found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_bnf_base) &
          CALL finish('init_lctlib','No data for bnf_base found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_vmax_symb_bnf) &
          CALL finish('init_lctlib','No data for vmax_symb_bnf found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_km_symb_bnf) &
          CALL finish('init_lctlib','No data for km_symb_bnf found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_lambda_est_light) &
          CALL finish('init_lctlib','No data for lambda_est_light found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_k_est_light) &
          CALL finish('init_lctlib','No data for k_est_light found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_seed_size) &
          CALL finish('init_lctlib','No data for seed_size found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_k1_mort_greff) &
          CALL finish('init_lctlib','No data for k1_mort_greff found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_beta_soil_flush) &
          CALL finish('init_lctlib','No data for beta_soil_flush found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_beta_soil_senescence) &
          CALL finish('init_lctlib','No data for beta_soil_senescence found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_gdd_req_max) &
          CALL finish('init_lctlib','No data for gdd_req_max found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_k_gdd_dormance) &
          CALL finish('init_lctlib','No data for k_gdd_dormance found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_t_air_senescence) &
          CALL finish('init_lctlib','No data for t_air_senescence found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_min_leaf_age) &
          CALL finish('init_lctlib','No data for min_leaf_age found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_frac_sapwood_branch) &
          CALL finish('init_lctlib','No data for frac_sapwood_branch found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_wood_density) &
          CALL finish('init_lctlib','No data for wood_density found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_k_latosa) &
          CALL finish('init_lctlib','No data for k_latosa found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_k_crtos) &
          CALL finish('init_lctlib','No data for k_crtos found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_k_rtos) &
          CALL finish('init_lctlib','No data for k_rtos found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_k2_fruit_alloc) &
          CALL finish('init_lctlib','No data for k2_fruit_alloc found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_allom_k1) &
          CALL finish('init_lctlib','No data for allom_k1 found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_allom_k2) &
          CALL finish('init_lctlib','No data for allom_k2 found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_phi_leaf_min) &
          CALL finish('init_lctlib','No data for phi_leaf_min found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_k_root) &
          CALL finish('init_lctlib','No data for k_root found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_k_sapwood) &
          CALL finish('init_lctlib','No data for k_sapwood found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_c0_allom) &
          CALL finish('init_lctlib','No data for c0_allom found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_fstore_target) &
          CALL finish('init_lctlib','No data for fstore_target found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_k_root_dist) &
          CALL finish('init_lctlib','No data for k_root_dist found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_k_som_fast_init) &
          CALL finish('init_lctlib','No data for k_som_fast_init found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_k_som_slow_init) &
          CALL finish('init_lctlib','No data for k_som_slow_init found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_AlbedoLitterVIS) &
          CALL finish('init_lctlib','No data for AlbedoLitterVIS found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_AlbedoLitterNIR) &
          CALL finish('init_lctlib','No data for AlbedoLitterNIR found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_VegRoughness) &
          CALL finish('init_lctlib','No data for VegRoughness found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_MinVegRoughness) &
          CALL finish('init_lctlib','No data for MinVegRoughness found in '//TRIM(lctlib_file_name))
        IF(.NOT. exists_MaxVegRoughness) &
          CALL finish('init_lctlib','No data for MaxVegRoughness found in '//TRIM(lctlib_file_name))
       !END IF
       ! quincy end

       CLOSE(unit=lctlib_file_unit) ! Closing access to lctlib file

       !------------------------------------------------------------------------------------------
       !>Modify the input values of some quincy lctlib parameters
       !>
       !>  loop over vegetation lct:
       !>    'jsbach_lite' & 'jsbach_pfts': calculations are applied to vegetation only; the lct = 1 is glacier (and has no values, i.e. is 0, in the quincy lctlib part)!
       !>    'quincy_eight_pfts':           calculations are applied to all 8 PFT
       !IF (model_scheme_char == "jsbach_quincy" .OR. model_scheme_char == "quincy") THEN
       !  IF (usecase == 'quincy_eight_pfts') THEN
           j = 1
       !  ELSE
       !    j = 2
       !  END IF

       ! sla
         DO i = j,nlct
             lctlib(i)%sla = 1._wp / carbon_per_dryweight_leaf / 1000.0_wp * molar_mass_C * lctlib(i)%sla
         END DO
       ! cn_leaf
         DO i = j,nlct
           lctlib(i)%cn_leaf = carbon_per_dryweight_leaf * 1000._wp * molar_mass_N / molar_mass_C / lctlib(i)%cn_leaf
         END DO
       ! cn_leaf_min
         DO i = j,nlct
           lctlib(i)%cn_leaf_min = carbon_per_dryweight_leaf * 1000._wp * molar_mass_N / molar_mass_C / lctlib(i)%cn_leaf_min
         END DO
       ! cn_leaf_max
         DO i = j,nlct
           lctlib(i)%cn_leaf_max = carbon_per_dryweight_leaf * 1000._wp * molar_mass_N / molar_mass_C / lctlib(i)%cn_leaf_max
         END DO
       ! np_leaf
         DO i = j,nlct
           lctlib(i)%np_leaf = carbon_per_dryweight_leaf * 1000._wp * molar_mass_P / molar_mass_C / lctlib(i)%cn_leaf / &
                                lctlib(i)%np_leaf
         END DO
       ! np_leaf_min
         DO i = j,nlct
           lctlib(i)%np_leaf_min = carbon_per_dryweight_leaf * 1000._wp * molar_mass_P / molar_mass_C / lctlib(i)%cn_leaf / &
                                   lctlib(i)%np_leaf_min
         END DO
       ! np_leaf_max
         DO i = j,nlct
           lctlib(i)%np_leaf_max = carbon_per_dryweight_leaf * 1000._wp * molar_mass_P / molar_mass_C / lctlib(i)%cn_leaf / &
                                   lctlib(i)%np_leaf_max
         END DO
       ! fn_oth_min
         DO i = j,nlct
           lctlib(i)%fn_oth_min = lctlib(i)%k0_fn_struc - 5.0_wp/molar_mass_n * 0.5_wp * 1.e3_wp * 71.4_wp * molar_mass_N / 1.e6_wp
         END DO
       ! tau_leaf
         DO i = j,nlct
           lctlib(i)%tau_leaf = lctlib(i)%tau_leaf / 12.0_wp
         END DO
       ! bnf_base
         DO i = j,nlct
          lctlib(i)%bnf_base = lctlib(i)%bnf_base / molar_mass_N / one_day / one_year * 1.e6_wp
         END DO
       ! vmax_symb_bnf
         DO i = j,nlct
          lctlib(i)%vmax_symb_bnf = lctlib(i)%vmax_symb_bnf / (molar_mass_N / molar_mass_C) / one_day / one_year * 1.e6_wp
         END DO
       ! km_symb_bnf
         DO i = j,nlct
          lctlib(i)%km_symb_bnf = lctlib(i)%km_symb_bnf / (molar_mass_C / molar_mass_N)
         END DO
       ! t_air_senescence
         DO i = j,nlct
           lctlib(i)%t_air_senescence = lctlib(i)%t_air_senescence + Tzero
         END DO
       ! wood_density
         DO i = j,nlct
           lctlib(i)%wood_density = lctlib(i)%wood_density * 1.e6_wp/ molar_mass_C
         END DO
       ! k_root
         DO i = j,nlct
           lctlib(i)%k_root = lctlib(i)%k_root * molar_mass_C * 1.e-10_wp
         END DO
       ! k_sapwood
         DO i = j,nlct
           lctlib(i)%k_sapwood = 1.e-3_wp * lctlib(i)%k_sapwood
         END DO
        ! k_rtos
         DO i = j,nlct
           SELECT CASE (lctlib(i)%growthform)
             CASE (itree)
               lctlib(i)%k_rtos = 1.75_wp * (lctlib(4)%sla * lctlib(i)%wood_density) / lctlib(i)%k_latosa
               !lctlib(i)%k_rtos = SQRT(lctlib(i)%k_root / lctlib(i)%k_sapwood * &
               !     lctlib(i)%tau_sap_wood / lctlib(i)%tau_fine_root * lctlib(i)%wood_density)
             CASE (igrass,icrop)
               lctlib(i)%k_rtos = 0.5_wp / sm2lm_grass
           END SELECT
         END DO
         ! c0_allom
         DO i = j,nlct
           IF(lctlib(i)%growthform == itree) THEN
             lctlib(i)%c0_allom = SQRT(lctlib(i)%k_root / lctlib(i)%k_sapwood * &
                                       lctlib(i)%tau_sap_wood / lctlib(i)%tau_fine_root * lctlib(i)%wood_density) * &
                                  lctlib(i)%tau_fine_root / lctlib(i)%tau_sap_wood
           ELSE
             lctlib(i)%c0_allom = 0.0_wp
           END IF
         END DO
       ! derive C4flag from ps_pathway
         DO i = 1,nlct
            IF ( lctlib(i)%ps_pathway == ic4phot ) THEN
              lctlib(i)%C4flag = .TRUE.
            ELSE
              lctlib(i)%C4flag = .FALSE.
            END IF
         END DO

       !END IF
       ! quincy

       !------------------------------------------------------------------------------------------

       ! --- translate landcover class information into landcover class flags

       lctlib(1:nlct)%BareSoilFlag   = .FALSE.
       lctlib(1:nlct)%GlacierFlag    = .FALSE.
       lctlib(1:nlct)%LakeFlag       = .FALSE.
       lctlib(1:nlct)%ForestFlag     = .FALSE.
       lctlib(1:nlct)%GrassFlag      = .FALSE.
       lctlib(1:nlct)%CropFlag       = .FALSE.
       lctlib(1:nlct)%PastureFlag    = .FALSE.
       lctlib(1:nlct)%NaturalVegFlag = .FALSE.

       DO i=1,nlct
          SELECT CASE (lctlib(i)%LandcoverClass)
          CASE (0)
             lctlib(i)%BareSoilFlag = .TRUE.
          CASE (1)
             lctlib(i)%GlacierFlag = .TRUE.
          CASE (2)
             lctlib(i)%LakeFlag = .TRUE.
          CASE (3)
             lctlib(i)%ForestFlag = .TRUE.
             lctlib(i)%NaturalVegFlag = .TRUE.
          CASE (4)
             lctlib(i)%GrassFlag = .TRUE.
             lctlib(i)%NaturalVegFlag = .TRUE.
          CASE (5)
             lctlib(i)%NaturalVegFlag = .TRUE.
          CASE (6)
             lctlib(i)%CropFlag = .TRUE.
          CASE (7)
             lctlib(i)%PastureFlag = .TRUE.
          CASE default
             WRITE(message_text,*) 'LandcoverClass ', lctlib(i)%LandcoverClass,' not allowed in ', TRIM(lctlib_file_name)
             CALL finish('init_lctlib', message_text)
          END SELECT
       END DO

       npft = 0
       DO i=1,nlct
         IF(lctlib(i)%naturalVegFlag .OR. lctlib(i)%cropFlag .OR. lctlib(i)%PastureFlag) THEN

           npft = npft + 1 ! count PFTs


!             ! Check consistency with PhenologyTypes
!             IF(lctlib(i)%naturalVegFlag .AND. &
!                  (lctlib(i)%PhenologyType .LE.0 .OR. lctlib(i)%PhenologyType > 4) ) THEN
!                CALL finish('init_lctlib',&
!                     'Inconsistent entries for one LandcoverClass and PhenologyType in '//TRIM(lctlib_file_name))
!             END IF
!             IF(lctlib(i)%CropFlag .AND. .NOT. lctlib(i)%PhenologyType==5 ) THEN
!                CALL finish('init_lctlib',&
!                     'Inconsistent entries for LandcoverClass "crops" and PhenologyType in '//TRIM(lctlib_file_name))
!             END IF

         END IF

!          !! --- check consistency between landcover classes and specifications for the dynamic vegetation

!          IF (lctlib(i)%dynamic_pft  .AND. .NOT. lctlib(i)%NaturalVegFlag) THEN
!             CALL finish('init_lctlib',&
!                     'Inconsistent entries: Every DYNAMIC_PFT must belong to LandcoverClass "natural" in '//TRIM(lctlib_file_name))
!          END IF

!          IF(lctlib(i)%woody_pft .AND. .NOT. lctlib(i)%NaturalVegFlag) THEN
!             CALL finish('init_lctlib',&
!                     'Inconsistent entries: Every WOODY_PFT must belong to LandcoverClass "natural" in '//TRIM(lctlib_file_name))
!          END IF

!          IF( lctlib(i)%GrassFlag .NEQV. (lctlib(i)%dynamic_pft .AND. .NOT. lctlib(i)%woody_PFT )) THEN
!             CALL finish('init_lctlib',&
!                  'Combination of DYNAMIC_PFT and WOODY_PFT is not consistent with LandcoverClass "grasses" in '//&
!                                                                                                     TRIM(lctlib_file_name))
!          END IF
       END DO

!    END IF !p_parallel_io


!    IF(my_process_is_mpi_parallel()) THEN
!       CALL p_bcast(npft, p_io)
!       CALL p_bcast(lctlib(:)%LctNumber, p_io)
!       DO i=1,nlct
!          CALL p_bcast(lctlib(i)%LctName, p_io)
!       ENDDO
!       CALL p_bcast(lctlib(:)%LandcoverClass, p_io)
!       CALL p_bcast(lctlib(:)%NaturalVegFlag,p_io)
!       CALL p_bcast(lctlib(:)%ForestFlag,p_io)
!       CALL p_bcast(lctlib(:)%GrassFlag,p_io)
!       CALL p_bcast(lctlib(:)%CropFlag,p_io)
!       CALL p_bcast(lctlib(:)%PastureFlag,p_io)
!       CALL p_bcast(lctlib(:)%LakeFlag, p_io)
!       CALL p_bcast(lctlib(:)%GlacierFlag, p_io)
!       CALL p_bcast(lctlib(:)%BareSoilFlag, p_io)
!       CALL p_bcast(lctlib(:)%C4flag, p_io)
!       CALL p_bcast(lctlib(:)%AlbedoSnowVisMin, p_io)
!       CALL p_bcast(lctlib(:)%AlbedoSnowVisMax, p_io)
!       CALL p_bcast(lctlib(:)%AlbedoSnowNirMin, p_io)
!       CALL p_bcast(lctlib(:)%AlbedoSnowNirMax, p_io)
!       CALL p_bcast(lctlib(:)%AlbedoSnowMin, p_io)
!       CALL p_bcast(lctlib(:)%AlbedoSnowMax, p_io)
!       CALL p_bcast(lctlib(:)%AlbedoCanopyVIS, p_io)
!       CALL p_bcast(lctlib(:)%AlbedoCanopyNIR, p_io)
!       CALL p_bcast(lctlib(:)%AlbedoLitterVis, p_io)
!       CALL p_bcast(lctlib(:)%AlbedoLitterNir, p_io)
!       CALL p_bcast(lctlib(:)%NitrogenScalingFlag,p_io)
!       CALL p_bcast(lctlib(:)%C4flag, p_io)
!       CALL p_bcast(lctlib(:)%CarboxRate, p_io)
!       CALL p_bcast(lctlib(:)%ETransport, p_io)
!       CALL p_bcast(lctlib(:)%VegHeight, p_io)
!       CALL p_bcast(lctlib(:)%VegRoughness, p_io)
!       CALL p_bcast(lctlib(:)%MinVegRoughness, p_io)
!       CALL p_bcast(lctlib(:)%MaxVegRoughness, p_io)
!       CALL p_bcast(lctlib(:)%PhenologyType, p_io)
!       CALL p_bcast(lctlib(:)%CanopyResistanceMin, p_io)
!       CALL p_bcast(lctlib(:)%MaxLAI, p_io)
!       CALL p_bcast(lctlib(:)%StemArea, p_io)
!       CALL p_bcast(lctlib(:)%specificLeafArea_C, p_io)
!       CALL p_bcast(lctlib(:)%knorr_Tau_w, p_io)
!       CALL p_bcast(lctlib(:)%knorr_T_phi, p_io)
!       CALL p_bcast(lctlib(:)%knorr_T_r, p_io)
!       CALL p_bcast(lctlib(:)%knorr_Day_c, p_io)
!       CALL p_bcast(lctlib(:)%knorr_Day_r, p_io)
!       CALL p_bcast(lctlib(:)%knorr_k_l, p_io)
!       CALL p_bcast(lctlib(:)%knorr_leaf_growth_rate, p_io)
!       CALL p_bcast(lctlib(:)%knorr_max_lai, p_io)
!       CALL p_bcast(lctlib(:)%reserveC2leafC, p_io)
!       CALL p_bcast(lctlib(:)%frac_npp_2_woodPool, p_io)
!       CALL p_bcast(lctlib(:)%frac_npp_2_reservePool, p_io)
!       CALL p_bcast(lctlib(:)%frac_npp_2_exudates, p_io)
!       CALL p_bcast(lctlib(:)%frac_green_2_herbivory, p_io)
!       CALL p_bcast(lctlib(:)%tau_Cpool_litter_leaf, p_io)
!       CALL p_bcast(lctlib(:)%tau_Cpool_litter_wood, p_io)
!       CALL p_bcast(lctlib(:)%tau_Cpool_woods, p_io)
!       CALL p_bcast(lctlib(:)%LAI_shed_constant, p_io)
!       CALL p_bcast(lctlib(:)%frac_C_litter_green2atmos, p_io)
!       CALL p_bcast(lctlib(:)%Max_C_content_woods, p_io)
!       CALL p_bcast(lctlib(:)%ClumpinessFactor, p_io)
!       CALL p_bcast(lctlib(:)%dynamic_pft, p_io)
!       CALL p_bcast(lctlib(:)%woody_pft, p_io)
!       CALL p_bcast(lctlib(:)%pasture_pft, p_io)
!       CALL p_bcast(lctlib(:)%bclimit_min_cold_mmtemp, p_io)
!       CALL p_bcast(lctlib(:)%bclimit_max_cold_mmtemp, p_io)
!       CALL p_bcast(lctlib(:)%bclimit_max_warm_mmtemp, p_io)
!       CALL p_bcast(lctlib(:)%bclimit_min_temprange, p_io)
!       CALL p_bcast(lctlib(:)%bclimit_min_gdd, p_io)
!       CALL p_bcast(lctlib(:)%gdd_base, p_io)
!       CALL p_bcast(lctlib(:)%upper_tlim, p_io)
!       CALL p_bcast(lctlib(:)%frac_wood_2_onSite, p_io)
!       CALL p_bcast(lctlib(:)%frac_wood_2_paper, p_io)
!       CALL p_bcast(lctlib(:)%frac_wood_2_construction, p_io)
!       CALL p_bcast(lctlib(:)%moist_extinction, p_io)
!       CALL p_bcast(lctlib(:)%fuel_dens, p_io)
!       CALL p_bcast(lctlib(:)%flame_length_f, p_io)
!       CALL p_bcast(lctlib(:)%crown_length, p_io)
!       CALL p_bcast(lctlib(:)%bark_par1, p_io)
!       CALL p_bcast(lctlib(:)%bark_par2, p_io)
!       CALL p_bcast(lctlib(:)%RCK, p_io)
!       CALL p_bcast(lctlib(:)%mort_prob, p_io)
!       DO i=1,5
!         CALL p_bcast(lctlib(:)%LitVeg_coef(i), p_io)
!         CALL p_bcast(lctlib(:)%LeafLit_coef(i), p_io)
!         CALL p_bcast(lctlib(:)%WoodLit_coef(i), p_io)
!       END DO
!       CALL p_bcast(lctlib(:)%WoodLitterSize, p_io)
!    END IF

    IF(npft>nlct) CALL finish('init_lctlib','npft larger than nlct: '//int2string(npft)//' > '//int2string(nlct))

!    ! Throw away values for glacier (first column)
!    nlct = nlct - 1
    ALLOCATE(return_value(nlct))
    return_value(1:nlct) = lctlib(1:nlct)

  END FUNCTION Read_lctlib

END MODULE mo_jsb_lctlib_class
