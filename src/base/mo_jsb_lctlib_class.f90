!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!! 
!! @par Revision History
!! First implementation into SLM by Jan Engel MPI-BGC Jena (20170518)
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

    INTEGER  :: growthform 
    INTEGER  :: ps_pathway
    INTEGER  :: phenology_type
    REAL(wp) :: sla                 
    REAL(wp) :: sigma_vis 
    REAL(wp) :: sigma_nir 
    REAL(wp) :: omega_clumping 
    REAL(wp) :: crown_shape_factor 
    REAL(wp) :: cn_leaf             
    REAL(wp) :: cn_leaf_min
    REAL(wp) :: cn_leaf_max
    REAL(wp) :: np_leaf             
    REAL(wp) :: k0_fn_struc
    REAL(wp) :: fn_oth_min
    REAL(wp) :: t_jmax_opt
    REAL(wp) :: t_jmax_omega
    REAL(wp) :: g0
    REAL(wp) :: g1
    REAL(wp) :: gmin
    REAL(wp) :: lambda_soil_ps
    REAL(wp) :: weibull_k_soil_ps
    REAL(wp) :: lambda_soil_gs
    REAL(wp) :: weibull_k_soil_gs

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
                bnf_base

! Vegetation dynamics parameters
    REAL(wp) :: lambda_est_light  , &
                k_est_light       , &
                seed_size

! Phenology parameters
    REAL(wp) :: beta_soil_flush       , &
                beta_soil_senescence  , &
                gdd_t_air_threshold   , &
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
    REAL(wp) :: k_root_dist


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
    USE mo_util,    ONLY: tolower

    ! USE statements necessary for modifying input values of lctlib parameters
    USE mo_base_constants,          ONLY: one_day, one_year
    USE mo_shared_constants,        ONLY: molar_mass_C, molar_mass_N, molar_mass_P, Dwv, Dco2, Tzero 
    USE mo_veg_constants,           ONLY: carbon_per_dryweight_leaf, sm2lm_grass, igrass, itree

    CHARACTER(len=*), INTENT(in) :: lctlib_file_name
    TYPE(t_lctlib_element), POINTER    :: return_value(:)

    TYPE(t_lctlib_element), ALLOCATABLE :: lctlib(:)

    INTEGER, PARAMETER            :: lctlib_file_unit = 66

    CHARACTER(len=30)  :: key
    CHARACTER(len=256) :: line
    INTEGER            :: pos_comment, read_status
    INTEGER            :: pos,length
    CHARACTER(len=2)   :: blank_set = " "//achar(9) !! the blank characters: BLANK and TAB
    INTEGER,ALLOCATABLE:: itmp(:)                   !! temporary array used for input of logicals
    INTEGER            :: i
    INTEGER :: nlct, npft

    ! Determine name of landcover library file

    ! Read lctlib data and eventually allocate also memory for lctlib data on io-processor

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

          SELECT CASE (TRIM(key))
          CASE ('lctnumber')
            READ(line(pos:length),*) lctlib(1:nlct)%LctNumber
!            exist_LctNumber = .TRUE.
          CASE ('lctname')                                           ! Name of landcover type
            READ(line(pos:length),*) lctlib(1:nlct)%LctName
!            exist_LctName = .TRUE.
          CASE ('landcoverclass')                                    ! Landcover class
            READ(line(pos:length),*) lctlib(1:nlct)%LandcoverClass
!            exist_LandcoverClass = .TRUE.

          CASE ('growthform')
            READ(line(pos:length),*) lctlib(1:nlct)%growthform
          CASE ('ps_pathway')
            READ(line(pos:length),*) lctlib(1:nlct)%ps_pathway
          CASE ('phenology_type')
            READ(line(pos:length),*) lctlib(1:nlct)%phenology_type
          CASE ('sla')
            READ(line(pos:length),*) lctlib(1:nlct)%sla
          CASE ('sigma_vis')
            READ(line(pos:length),*) lctlib(1:nlct)%sigma_vis
          CASE ('sigma_nir')
            READ(line(pos:length),*) lctlib(1:nlct)%sigma_nir
          CASE ('omega_clumping')
            READ(line(pos:length),*) lctlib(1:nlct)%omega_clumping
          CASE ('crown_shape_factor')
            READ(line(pos:length),*) lctlib(1:nlct)%crown_shape_factor
          CASE ('cn_leaf')
            READ(line(pos:length),*) lctlib(1:nlct)%cn_leaf
          CASE ('cn_leaf_min')
            READ(line(pos:length),*) lctlib(1:nlct)%cn_leaf_min
          CASE ('cn_leaf_max')
            READ(line(pos:length),*) lctlib(1:nlct)%cn_leaf_max
          CASE ('np_leaf')
            READ(line(pos:length),*) lctlib(1:nlct)%np_leaf
          CASE ('k0_fn_struc')
            READ(line(pos:length),*) lctlib(1:nlct)%k0_fn_struc
          CASE ('fn_oth_min')
            READ(line(pos:length),*) lctlib(1:nlct)%fn_oth_min
          CASE ('t_jmax_opt')
            READ(line(pos:length),*) lctlib(1:nlct)%t_jmax_opt
          CASE ('t_jmax_omega')
            READ(line(pos:length),*) lctlib(1:nlct)%t_jmax_omega
          CASE ('g0')
            READ(line(pos:length),*) lctlib(1:nlct)%g0
          CASE ('g1')
            READ(line(pos:length),*) lctlib(1:nlct)%g1
          CASE ('gmin')
            READ(line(pos:length),*) lctlib(1:nlct)%gmin
          CASE ('lambda_soil_ps')
            READ(line(pos:length),*) lctlib(1:nlct)%lambda_soil_ps
          CASE ('weibull_k_soil_ps')
            READ(line(pos:length),*) lctlib(1:nlct)%weibull_k_soil_ps
          CASE ('lambda_soil_gs')
            READ(line(pos:length),*) lctlib(1:nlct)%lambda_soil_gs
          CASE ('weibull_k_soil_gs')
            READ(line(pos:length),*) lctlib(1:nlct)%weibull_k_soil_gs

          CASE ('tau_leaf')
            READ(line(pos:length),*) lctlib(1:nlct)%tau_leaf
          CASE ('tau_fine_root')
            READ(line(pos:length),*) lctlib(1:nlct)%tau_fine_root
          CASE ('tau_coarse_root')
            READ(line(pos:length),*) lctlib(1:nlct)%tau_coarse_root
          CASE ('tau_branch')
            READ(line(pos:length),*) lctlib(1:nlct)%tau_branch
          CASE ('tau_sap_wood')
            READ(line(pos:length),*) lctlib(1:nlct)%tau_sap_wood
          CASE ('tau_fruit')
            READ(line(pos:length),*) lctlib(1:nlct)%tau_fruit
          CASE ('tau_seed_litter')
            READ(line(pos:length),*) lctlib(1:nlct)%tau_seed_litter
          CASE ('tau_seed_est')
            READ(line(pos:length),*) lctlib(1:nlct)%tau_seed_est
          CASE ('tau_mycorrhiza')
            READ(line(pos:length),*) lctlib(1:nlct)%tau_mycorrhiza

          CASE ('vmax_uptake_n')
            READ(line(pos:length),*) lctlib(1:nlct)%vmax_uptake_n
          CASE ('vmax_uptake_p')
            READ(line(pos:length),*) lctlib(1:nlct)%vmax_uptake_p
          CASE ('bnf_base')
            READ(line(pos:length),*) lctlib(1:nlct)%bnf_base

          CASE ('lambda_est_light')
            READ(line(pos:length),*) lctlib(1:nlct)%lambda_est_light
          CASE ('k_est_light')
            READ(line(pos:length),*) lctlib(1:nlct)%k_est_light
          CASE ('seed_size')
            READ(line(pos:length),*) lctlib(1:nlct)%seed_size

          CASE ('beta_soil_flush')
            READ(line(pos:length),*) lctlib(1:nlct)%beta_soil_flush
          CASE ('beta_soil_senescence')
            READ(line(pos:length),*) lctlib(1:nlct)%beta_soil_senescence
          CASE ('gdd_t_air_threshold')
            READ(line(pos:length),*) lctlib(1:nlct)%gdd_t_air_threshold
          CASE ('gdd_req_max')
            READ(line(pos:length),*) lctlib(1:nlct)%gdd_req_max
          CASE ('k_gdd_dormance')
            READ(line(pos:length),*) lctlib(1:nlct)%k_gdd_dormance
          CASE ('t_air_senescence')
            READ(line(pos:length),*) lctlib(1:nlct)%t_air_senescence
          CASE ('min_leaf_age')
            READ(line(pos:length),*) lctlib(1:nlct)%min_leaf_age

          CASE ('frac_sapwood_branch')
            READ(line(pos:length),*) lctlib(1:nlct)%frac_sapwood_branch
          CASE ('wood_density')
            READ(line(pos:length),*) lctlib(1:nlct)%wood_density
          CASE ('k_latosa')
            READ(line(pos:length),*) lctlib(1:nlct)%k_latosa
          CASE ('k_crtos')
            READ(line(pos:length),*) lctlib(1:nlct)%k_crtos
          CASE ('k_rtos')
            READ(line(pos:length),*) lctlib(1:nlct)%k_rtos
          CASE ('k2_fruit_alloc')
            READ(line(pos:length),*) lctlib(1:nlct)%k2_fruit_alloc
          CASE ('allom_k1')
            READ(line(pos:length),*) lctlib(1:nlct)%allom_k1
          CASE ('allom_k2')
            READ(line(pos:length),*) lctlib(1:nlct)%allom_k2
          CASE ('phi_leaf_min')
            READ(line(pos:length),*) lctlib(1:nlct)%phi_leaf_min
          CASE ('k_root')
            READ(line(pos:length),*) lctlib(1:nlct)%k_root
          CASE ('k_sapwood')
            READ(line(pos:length),*) lctlib(1:nlct)%k_sapwood
          CASE ('c0_allom')
            READ(line(pos:length),*) lctlib(1:nlct)%c0_allom
          CASE ('fstore_target')
            READ(line(pos:length),*) lctlib(1:nlct)%fstore_target

          CASE ('k_root_dist')
            READ(line(pos:length),*) lctlib(1:nlct)%k_root_dist



          CASE default
             ! nothing to do
          END SELECT
       END DO
       DEALLOCATE(itmp)


       CLOSE(unit=lctlib_file_unit) !! Closing access to lctlib file

        !------------------------------------------------------------------------------------------
        !> Modify the input values of some lctlib parameters
        !! SLM specific code !!

        !! sla
          DO i=1,nlct
              lctlib(i)%sla = 1._wp / carbon_per_dryweight_leaf / 1000.0_wp * molar_mass_C * lctlib(i)%sla
          END DO
        !! cn_leaf
          DO i=1,nlct
            lctlib(i)%cn_leaf = carbon_per_dryweight_leaf * 1000._wp * molar_mass_N / molar_mass_C / lctlib(i)%cn_leaf
          END DO
        !! cn_leaf_min
          DO i=1,nlct
            lctlib(i)%cn_leaf_min = carbon_per_dryweight_leaf * 1000._wp * molar_mass_N / molar_mass_C / lctlib(i)%cn_leaf_min
          END DO
        !! cn_leaf_max
          DO i=1,nlct
            lctlib(i)%cn_leaf_max = carbon_per_dryweight_leaf * 1000._wp * molar_mass_N / molar_mass_C / lctlib(i)%cn_leaf_max
          END DO
        !! np_leaf
          DO i=1,nlct
            lctlib(i)%np_leaf = carbon_per_dryweight_leaf * 1000._wp * molar_mass_P / molar_mass_C / lctlib(i)%cn_leaf / &
                                 lctlib(i)%np_leaf
          END DO
        !! fn_oth_min
          DO i=1,nlct
            lctlib(i)%fn_oth_min = lctlib(i)%k0_fn_struc - 5.0_wp/molar_mass_n * 0.5_wp * 1.e3_wp * 71.4_wp * molar_mass_N / 1.e6_wp
          END DO
        !! tau_leaf
          DO i=1,nlct
            lctlib(i)%tau_leaf = lctlib(i)%tau_leaf / 12.0_wp
          END DO
        !! bnf_base
          DO i=1,nlct
           lctlib(i)%bnf_base = lctlib(i)%bnf_base / molar_mass_N / one_day / one_year * 1.e6_wp
          END DO
        !! gdd_t_air_threshold
          DO i=1,nlct
            lctlib(i)%gdd_t_air_threshold = lctlib(i)%gdd_t_air_threshold + Tzero
          END DO
        !! t_air_senescence
          DO i=1,nlct
            lctlib(i)%t_air_senescence = lctlib(i)%t_air_senescence + Tzero
          END DO
        !! wood_density
          DO i=1,nlct
            lctlib(i)%wood_density = lctlib(i)%wood_density * 1.e6_wp/ molar_mass_C
          END DO
        !! k_root
          DO i=1,nlct
            lctlib(i)%k_root = lctlib(i)%k_root * molar_mass_C * 1.e-10_wp
          END DO
        !! k_sapwood
          DO i=1,nlct
            lctlib(i)%k_sapwood = 1.e-3_wp * lctlib(i)%k_sapwood
         END DO
         !! k_rtos
         DO i=1,nlct
            IF(lctlib(i)%growthform /= igrass) THEN
               lctlib(i)%k_rtos = 1.75_wp * (lctlib(4)%sla * lctlib(i)%wood_density) / lctlib(i)%k_latosa
               !lctlib(i)%k_rtos = SQRT(lctlib(i)%k_root / lctlib(i)%k_sapwood * & 
               !     lctlib(i)%tau_sap_wood / lctlib(i)%tau_fine_root * lctlib(i)%wood_density)
            ELSE
               lctlib(i)%k_rtos = 0.5_wp / sm2lm_grass
            END IF
         END DO
        !! c0_allom
          DO i=1,nlct
            IF(lctlib(i)%growthform == itree) THEN
              lctlib(i)%c0_allom = SQRT(lctlib(i)%k_root / lctlib(i)%k_sapwood * & 
                                        lctlib(i)%tau_sap_wood / lctlib(i)%tau_fine_root * lctlib(i)%wood_density) * & 
                                   lctlib(i)%tau_fine_root / lctlib(i)%tau_sap_wood
            ELSE
              lctlib(i)%c0_allom = 0.0_wp
            END IF
          END DO

        !------------------------------------------------------------------------------------------

       !! --- translate landcover class information into landcover class flags

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
