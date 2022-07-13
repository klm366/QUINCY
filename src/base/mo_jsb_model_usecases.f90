!>
!! This moduled defines various model usecases that can be enabled individually via namelist.
!! A model usecase defines:
!!    * the tiles that are enabled in the model and their nested structure
!!    * the processes that are active on each tile or its leaves
!!    * the fraction of each tile
!!    * the land-cover-type of each tile
!! 
!! @brief Contains methods for initialization of model usecases.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_jsb_model_usecases
#ifndef __NO_JSBACH__

  USE mo_exception,           ONLY: message, finish

  USE mo_jsb_model_class,     ONLY: t_jsb_model, MODEL_JSBACH, MODEL_QUINCY
  USE mo_jsb_process_class,   ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, RAD_, HD_, ASSIMI_, PHENO_, CARBON_, DISTURB_, &
    &                               FUEL_, ALCC_, FLCC_, WLCC_, NLCC_, PPLCC_, TLCC_, TCQ_, VEG_, SB_, SPQ_
  USE mo_jsb_process_class,   ONLY: ON_LEAFS_, ON_TILE_, AGGREGATE_ !, ON_SUBTREE_
  USE mo_jsb_tile,            ONLY: t_jsb_tile
  USE mo_jsb_lct_class,       ONLY: t_jsb_lct, &
    &                               LAND_TYPE, VEG_TYPE, GLACIER_TYPE, LAKE_TYPE, & !, BARE_TYPE   ! define the type of the tile
    &                               UPLAND_TYPE, LOWLAND_TYPE
!  USE mo_test_pools,          ONLY: test_pools                                                  ! to debug and test pools implementation

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_usecase

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_model_usecases'

CONTAINS

  ! ======================================================================================================= !
  !>initialize the usecase that determines the tiles and what processes are running where
  !>
  SUBROUTINE init_usecase(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_usecase'

    IF (model%config%use_lakes) THEN
      CALL message(TRIM(routine), 'Initializing model usecase "'//TRIM(model%config%usecase)//'" with lakes')
    ELSE
      CALL message(TRIM(routine), 'Initializing model usecase "'//TRIM(model%config%usecase)//'"')
    END IF

    ! select usecase based on namelist option
    SELECT CASE (TRIM(model%config%usecase))
    ! jsbach4 physics (with only SSE_ and PHENO_ at the veg_tile)
    CASE ('jsbach_lite')
      CALL init_usecase_lite(model)
    ! jsbach4 11 PFTs
    CASE ('jsbach_pfts')
      CALL init_usecase_pfts(model)
    CASE ('jsbach_hydrotiles')
      CALL init_usecase_hydrotiles(model)
    ! QUINCY 8 PFTs
    CASE ('quincy_eight_pfts')
      CALL init_usecase_quincy_eight_pfts(model)
    ! catch wrong usecase setting
    CASE DEFAULT
      CALL finish(TRIM(routine), 'No such usecase defined.')
    END SELECT

  END SUBROUTINE init_usecase

  ! ======================================================================================================= !
  !>init routine for usecase lite
  !>
  SUBROUTINE init_usecase_lite(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER :: box_tile, lake_tile, land_tile, glacier_tile, veg_tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_usecase_lite'

    ! overview of tile structure:
    !
    !  box
    !   |
    !   |-------
    !   |      |
    !  land  lakes
    !   |
    !   |-------------
    !   |            |      
    !  vegetation  glacier 

    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', &
      & processes      =(/A2L_,     L2A_,     RAD_,      HYDRO_,    TURB_    /), &
      & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_/), &
      & model_id=model%id, fract_varname='notsea')
    model%top => box_tile

#ifndef __NO_JSBACH_HD__
    IF (model%processes(HD_)%p%config%active) THEN
      CALL box_tile%Set_process_action(HD_, ON_TILE_)
    END IF
#endif

      ! Create sub-tiles
      IF (model%config%use_lakes) THEN
        lake_tile => t_jsb_tile('lake', 'Lake tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAKE_TYPE), parent=box_tile, fract_varname='fract_lake')
        land_tile => t_jsb_tile('land', 'Land tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAND_TYPE), parent=box_tile)
      ELSE
        land_tile => t_jsb_tile('land', 'Land tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAND_TYPE), &
          & parent=box_tile, &
          & fract_varname='equal')   ! Fraction = 1.
      END IF

        IF (model%config%use_glacier) THEN
          glacier_tile => t_jsb_tile('glac', 'Glacier tile', processes=(/SSE_/), &
            & lct=t_jsb_lct(GLACIER_TYPE), parent=land_tile, &
            & fract_varname='fract_glac')
        END IF

          veg_tile => t_jsb_tile('veg', 'Veg tile', processes=(/SSE_, PHENO_/), &
            & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
            & fract_varname='fract_veg')

    ! Next 3 lines only to avoid compiler warnings about pointers not being dereferenced
    IF (ASSOCIATED(glacier_tile)) CONTINUE
    IF (ASSOCIATED(lake_tile)) CONTINUE
    IF (ASSOCIATED(veg_Tile)) CONTINUE

    NULLIFY(box_tile, lake_tile, glacier_tile, land_tile, veg_tile)

  END SUBROUTINE init_usecase_lite

  ! ======================================================================================================= !
  !>init routine for usecase pfts
  !>
  SUBROUTINE init_usecase_pfts(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER :: box_tile, lake_tile, land_tile, glacier_tile, veg_tile, pft_tile

    INTEGER :: i

    INTEGER, PARAMETER :: npft = 11

    CHARACTER(len=5), PARAMETER :: pft_shortnames(npft) = [character(len=5) ::  &
      & 'pft01', 'pft02', 'pft03', 'pft04', 'pft05', 'pft06', &
      & 'pft07', 'pft08', 'pft09', 'pft10', 'pft11' ]

    CHARACTER(len=30), PARAMETER :: pft_longnames(npft) = [character(len=30) :: &
      & 'Tropical broadleaf evergreen', &
      & 'Tropical broadleaf deciduous', &
      & 'Extra-tropical evergreen',     &
      & 'Extra-tropical deciduous',     &
      & 'Raingreen shrubs',             &
      & 'Deciduous shrubs',             &
      & 'C3 grass',                     &
      & 'C4 grass',                     &
      & 'C3 pasture',                   &
      & 'C4 pasture',                   &
      & 'C3 crops'                    ]

    ! Index of PFT into lctlib columns
    ! id is one less than in lctlib file because the first column (glacier) has been thrown away
    INTEGER, PARAMETER :: pft_ids(npft) = (/ 1,2,3,4,9,10,11,12,14,15,19 /)

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_usecase_pfts'

    ! overview of tile structure:
    !
    !  box
    !   |
    !   |-------
    !   |      |
    !  land  lakes
    !   |
    !   |-------------
    !   |            |      
    !  vegetation  glacier 
    !   |
    !   |-----------------
    !   |      |      |
    !  pft01  pft02  pftXX

    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', &
      & processes      =(/A2L_,     L2A_,     RAD_  ,    TURB_    /), &
      & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_/), &
      & model_id=model%id, fract_varname='notsea')
    model%top => box_tile

#ifndef __NO_JSBACH_HD__
    IF (model%processes(HD_)%p%config%active) THEN
      CALL box_tile%Set_process_action(HD_, ON_TILE_)
    END IF
#endif

    IF (model%processes(PPLCC_)%p%config%active) THEN
      CALL box_tile%Set_process_action(PPLCC_, ON_LEAFS_)
    END IF

      IF (model%config%use_lakes) THEN
        lake_tile => t_jsb_tile('lake', 'Lake tile', &
          & processes=(/SEB_, HYDRO_/), process_actions=(/ON_TILE_, ON_TILE_/), &
          & lct=t_jsb_lct(LAKE_TYPE), parent=box_tile)
        land_tile => t_jsb_tile('land', 'Land tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAND_TYPE), parent=box_tile)
      ELSE
        land_tile => t_jsb_tile('land', 'Land tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAND_TYPE), parent=box_tile, &
          & fract_varname='equal')   ! Fraction = 1.
      END IF

        IF (model%config%use_glacier) THEN
          glacier_tile => t_jsb_tile('glac', 'Glacier tile', &
            & processes = (/HYDRO_, SSE_/), &
            & process_actions = (/ON_TILE_, ON_TILE_/), &
            & lct=t_jsb_lct(GLACIER_TYPE), parent=land_tile)
            ! & lct=t_jsb_lct(GLACIER_TYPE, lib_id=1), parent=land_tile)
        END IF

        !>veg_tile
        !!
        ! default jsb4: HYDRO_ & SSE_ are active at veg tile only to save runtime
        SELECT CASE (model%config%model_scheme)
        CASE (MODEL_JSBACH)
          veg_tile => t_jsb_tile('veg', 'Veg tile', &
            & processes       = (/HYDRO_,   SSE_,     ASSIMI_,   PHENO_/),    &
            & process_actions = (/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_/), &
            & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
            & fract_varname='fract_veg')

          IF (model%processes(ALCC_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(ALCC_, ON_TILE_, AGGREGATE_)
          END IF
          IF (model%processes(FLCC_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(FLCC_, ON_TILE_, AGGREGATE_)
          END IF
          IF (model%processes(WLCC_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(WLCC_, ON_TILE_, AGGREGATE_)
          END IF
          !
          IF (model%processes(NLCC_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(NLCC_, ON_TILE_, AGGREGATE_)
          END IF
        ! using quincy: 
        !  HYDRO_ is active at PFT tiles to make all the HYDRO_ variables available to the PFT tiles
        !    note: HYDRO_ results differ only marginally between ON_TILE_ & ON_LEAFS_ (tested by Reiner)
        !  SSE_ runs at PFT tiles instead of VEG tile, which is fine with the jsb4 tasks (info by Reiner)
        !    and needed by the quincy tasks because the lctlib is used
        !  VEG_ is active
        !  SB_  is active
        CASE (MODEL_QUINCY)
          veg_tile => t_jsb_tile('veg', 'Veg tile', &
            & processes       = (/HYDRO_,    SSE_,      ASSIMI_,   PHENO_,    VEG_,      SB_,       SPQ_/), &
            & process_actions = (/ON_LEAFS_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_/), &
            & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
            & fract_varname='fract_veg')
        END SELECT

        !   
        IF (model%processes(TLCC_)%p%config%active) THEN
          CALL veg_tile%Set_process_action(TLCC_, ON_TILE_, AGGREGATE_)
        END IF

          ! R: TBD: DISTURB_ ist kein Task von CARBON_,
          !         da es die cover fractions der pfts verändert und nur als Nebeneffekt auch die C-Pools
          !         verändert. Er würde im Prinzip auch ablaufen ohne den Prozess CARBON_ und würde auch
          !         auf alle anderen Stoffpools wirken. Außerdem soll CARBON_ später wie NITROGEN_
          !         oder andere Stoffpools gehandelt werden. Diese Prozesse haben bis
          !         auf ihre Fixierung kaum eigene Tasks sondern werden von anderen Prozessen verändert.
          !         => die Stoffpools sollten langfristig evtl. gar keine Prozesse sein
          !
          !         DISTURB_ braucht CARBON,
          !         da es C-Pools als Input für die Feuerberechnung braucht und sollte diese auch verändern
          !         wenn der CARBON Prozess existiert. Dies ist ein grundsätzliches Problem: Prozesse benötigen
          !         - oder sogar verändern - evtl. Variablen, die eigentlich zu einem anderen Prozess gehören.
          !         Das bedeutet der Speicher des anderen Prozesses muß angelegt werden auch wenn dieser
          !         Prozess gar nicht abläuft. Dazu müsste wohl am besten eine Prozedur die die Namelist-
          !         Einträge entwprechend umändert geschrieben werden oder wenn Usecases benutzt werden
          !         (wie hier) der entsprechende Prozess dennoch definiert werden, auch wenn in der Queue-
          !         Liste nachher kein Task dieses Prozesses auftaucht.
          !
          !         DISTURB_ braucht PHENO wg. veg_fract_correct s.o.
          !
          !         CARBON_ braucht DISTURB_ nicht!
          !
          !         Manche Prozeduren sollen nur ausgeführt werden wenn beide Prozesse aktiv sind.
          !         In DISTURB_ z.B. die dazugehörige C-Umverteilung bei Feuer und Winbruch.
          !         Dies habe ich mit "IF ( tile%Is_process_calculated(CARBON_) ) THEN" gemacht.
          !
          IF (model%processes(CARBON_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(CARBON_, ON_LEAFS_, AGGREGATE_)
          END IF
          IF (model%processes(TCQ_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(TCQ_, ON_LEAFS_, AGGREGATE_)
          END IF
          IF (model%processes(FUEL_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(FUEL_, ON_TILE_, AGGREGATE_)
          END IF
          IF (model%processes(DISTURB_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(DISTURB_, ON_LEAFS_, AGGREGATE_)
          END IF

            DO i=1,npft
              pft_tile => t_jsb_tile(pft_shortnames(i), TRIM(pft_longnames(i)), &
               & lct=t_jsb_lct(VEG_TYPE, lib_id=pft_ids(i)), parent=veg_tile)
            END DO

    ! ! Next 3 lines only to avoid compiler warnings about pointers not being dereferenced
    ! IF (ASSOCIATED(glacier_tile)) CONTINUE
    ! IF (ASSOCIATED(lake_tile)) CONTINUE
    ! IF (ASSOCIATED(pft_tile)) CONTINUE
        
    NULLIFY(box_tile, land_tile, veg_tile, pft_tile)
    IF (model%config%use_glacier) NULLIFY(glacier_tile)
    IF (model%config%use_lakes) NULLIFY(lake_tile)

    ! CALL test_pools()

  END SUBROUTINE init_usecase_pfts

    ! ======================================================================================================= !
  SUBROUTINE init_usecase_hydrotiles(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER :: box_tile, lake_tile, land_tile, glacier_tile, &
      & topo_tile, veg_tile, pft_tile

    INTEGER :: i

    INTEGER, PARAMETER :: npft = 11

    CHARACTER(len=5), PARAMETER :: pft_shortnames(npft) = [character(len=5) ::  &
      & 'pft01', 'pft02', 'pft03', 'pft04', 'pft05', 'pft06', &
      & 'pft07', 'pft08', 'pft09', 'pft10', 'pft11' ]

    CHARACTER(len=30), PARAMETER :: pft_longnames(npft) = [character(len=30) :: &
      & 'Tropical broadleaf evergreen', &
      & 'Tropical broadleaf deciduous', &
      & 'Extra-tropical evergreen',     &
      & 'Extra-tropical deciduous',     &
      & 'Raingreen shrubs',             &
      & 'Deciduous shrubs',             &
      & 'C3 grass',                     &
      & 'C4 grass',                     &
      & 'C3 pasture',                   &
      & 'C4 pasture',                   &
      & 'C3 crops'                    ]

    ! Index of PFT into lctlib columns
    ! id is one less than in lctlib file because the first column (glacier) has been thrown away
    INTEGER, PARAMETER :: pft_ids(npft) = (/ 1,2,3,4,9,10,11,12,14,15,19 /)

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_usecase_pfts'

    ! overview of tile structure:
    !
    !  box
    !   |
    !   |-------
    !   |      |
    !  land  lakes
    !   |
    !   |------------------------------
    !   |                |            |
    !  upland           lowland       glacier 
    !   |                |
    !  veg_upland       veg_lowland
    !   |                |
    !   |---------       |----------
    !   |        |       |         |
    !  pft01up  pftXXup pft01low  pftXXlow 
    !
    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', &
      & processes      =(/A2L_,     L2A_,     RAD_  ,    TURB_    /), &
      & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_/), &
      & model_id=model%id, fract_varname='notsea')
    model%top => box_tile

#ifndef __NO_JSBACH_HD__
    IF (model%processes(HD_)%p%config%active) THEN
      CALL box_tile%Set_process_action(HD_, ON_TILE_)
    END IF
#endif

    IF (model%processes(PPLCC_)%p%config%active) THEN
      CALL box_tile%Set_process_action(PPLCC_, ON_LEAFS_)
    END IF

      IF (model%config%use_lakes) THEN
        lake_tile => t_jsb_tile('lake', 'Lake tile', &
          & processes=(/SEB_, HYDRO_/), process_actions=(/ON_TILE_, ON_TILE_/), &
          & lct=t_jsb_lct(LAKE_TYPE), parent=box_tile)
        land_tile => t_jsb_tile('land', 'Land tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAND_TYPE), parent=box_tile)
      ELSE
        land_tile => t_jsb_tile('land', 'Land tile', &
          & processes=(/SEB_/), process_actions=(/ON_TILE_/), &
          & lct=t_jsb_lct(LAND_TYPE), parent=box_tile, &
          & fract_varname='equal')   ! Fraction = 1.
      END IF

        IF (model%config%use_glacier) THEN
          glacier_tile => t_jsb_tile('glac', 'Glacier tile', &
            & processes = (/HYDRO_, SSE_/), &
            & process_actions = (/ON_TILE_, ON_TILE_/), &
            & lct=t_jsb_lct(GLACIER_TYPE), parent=land_tile)
        END IF

        topo_tile => t_jsb_tile('upland', 'Upland tile', &
          ! Currently, topo_tiles have no own processes, later water lateral flow
          ! will be computed an these tiles
          ! & processes = (/HYDRO_, SSE_/), &
          ! & process_actions = (/ON_TILE_, ON_TILE_/), &
          & lct=t_jsb_lct(UPLAND_TYPE), parent=land_tile)

          ! HYDRO_ and SSE_ are active at veg tiles only to save runtime
          veg_tile => t_jsb_tile('veg_up', 'Veg tile (upland)', &
            & processes       = (/HYDRO_,   SSE_,     ASSIMI_,   PHENO_/),    &
            & process_actions = (/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_/), &
            & lct=t_jsb_lct(VEG_TYPE), parent=topo_tile,  &
            & fract_varname='fract_veg')

          IF (model%processes(ALCC_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(ALCC_, ON_TILE_, AGGREGATE_)
          END IF
          IF (model%processes(FLCC_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(FLCC_, ON_TILE_, AGGREGATE_)
          END IF
          IF (model%processes(WLCC_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(WLCC_, ON_TILE_, AGGREGATE_)
          END IF
          IF (model%processes(NLCC_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(NLCC_, ON_TILE_, AGGREGATE_)
          END IF
          IF (model%processes(TLCC_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(TLCC_, ON_TILE_, AGGREGATE_)
          END IF
          IF (model%processes(CARBON_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(CARBON_, ON_LEAFS_, AGGREGATE_)
          END IF
          IF (model%processes(TCQ_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(TCQ_, ON_LEAFS_, AGGREGATE_)
          END IF
          IF (model%processes(FUEL_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(FUEL_, ON_TILE_, AGGREGATE_)
          END IF
          IF (model%processes(DISTURB_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(DISTURB_, ON_LEAFS_, AGGREGATE_)
          END IF

            DO i=1,npft
              pft_tile => t_jsb_tile(pft_shortnames(i)//'_up', TRIM(pft_longnames(i)), &
               & lct=t_jsb_lct(VEG_TYPE, lib_id=pft_ids(i)), parent=veg_tile, &
               & fract_varname='fract_'//pft_shortnames(i))
            END DO

        topo_tile => t_jsb_tile('lowland', 'Lowland tile', &
          ! Currently, topo_tiles have no own processes, later water lateral flow
          ! will be computed an these tiles
          ! & processes = (/HYDRO_, SSE_/), &
          ! & process_actions = (/ON_TILE_, ON_TILE_/), &
          & lct=t_jsb_lct(LOWLAND_TYPE), parent=land_tile)

          veg_tile => t_jsb_tile('veg_low', 'Veg tile (Lowland)', &
            & processes       = (/HYDRO_,   SSE_,     ASSIMI_,   PHENO_/),    &
            & process_actions = (/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_/), &
            & lct=t_jsb_lct(VEG_TYPE), parent=topo_tile,  &
            & fract_varname='fract_veg')

          IF (model%processes(ALCC_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(ALCC_, ON_TILE_, AGGREGATE_)
          END IF
          IF (model%processes(FLCC_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(FLCC_, ON_TILE_, AGGREGATE_)
          END IF
          IF (model%processes(WLCC_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(WLCC_, ON_TILE_, AGGREGATE_)
          END IF
          IF (model%processes(NLCC_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(NLCC_, ON_TILE_, AGGREGATE_)
          END IF
          IF (model%processes(TLCC_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(TLCC_, ON_TILE_, AGGREGATE_)
          END IF
          IF (model%processes(CARBON_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(CARBON_, ON_LEAFS_, AGGREGATE_)
          END IF
          IF (model%processes(TCQ_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(TCQ_, ON_LEAFS_, AGGREGATE_)
          END IF
          IF (model%processes(FUEL_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(FUEL_, ON_TILE_, AGGREGATE_)
          END IF
          IF (model%processes(DISTURB_)%p%config%active) THEN
            CALL veg_tile%Set_process_action(DISTURB_, ON_LEAFS_, AGGREGATE_)
          END IF

            DO i=1,npft
              pft_tile => t_jsb_tile(pft_shortnames(i)//'_low', TRIM(pft_longnames(i)), &
               & lct=t_jsb_lct(VEG_TYPE, lib_id=pft_ids(i)), parent=veg_tile, &
               & fract_varname='fract_'//pft_shortnames(i))
            END DO
    ! ! Next 3 lines only to avoid compiler warnings about pointers not being dereferenced
    ! IF (ASSOCIATED(glacier_tile)) CONTINUE
    ! IF (ASSOCIATED(lake_tile)) CONTINUE
    ! IF (ASSOCIATED(pft_tile)) CONTINUE

    NULLIFY(box_tile, land_tile, topo_tile, veg_tile, pft_tile)
    IF (model%config%use_glacier) NULLIFY(glacier_tile)
    IF (model%config%use_lakes) NULLIFY(lake_tile)

  END SUBROUTINE init_usecase_hydrotiles

  ! ======================================================================================================= !
  !>init routine for usecase quincy_eight_pfts
  !>
  !>  works only in combination with 'model%config%model_scheme = MODEL_QUINCY'
  !>  code for other options not available in this usecase
  !>
  !>  the data/lctlib_quincy_nlct8.def is used for this usecase (hardcoded in: mo_jsb_base:jsbach_setup_tiles)
  !>
  SUBROUTINE init_usecase_quincy_eight_pfts(model)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: model

    CLASS(t_jsb_tile), POINTER    :: box_tile, lake_tile, land_tile, &                  ! tile used with this usecase
      &                              glacier_tile, veg_tile, pft_tile
    INTEGER                       :: i                                                  ! for looping
    INTEGER, PARAMETER            :: npft = 8                                           ! number of pft
    CHARACTER(len=5), PARAMETER   :: pft_shortnames(npft) = [character(len=5) :: &      ! pft short names
      &                              'pft01', 'pft02', 'pft03', 'pft04', &
      &                              'pft05', 'pft06', 'pft07', 'pft08' ]
    CHARACTER(len=30), PARAMETER  :: pft_longnames(npft) = [character(len=30) ::    &   ! pft long names
      &                              'TrBE: broad leaved evergreen (rainforest)',   &
      &                              'TeBE: broad leaved evergreen (xeric forest)', &
      &                              'TrBR: broad leaved deciduous (rain green)',   &
      &                              'TeBS: broad leaved deciduous (summer green)', &
      &                              'BNE:  needle-leaved evergreen',               &
      &                              'BNS:  needle-leaved deciduous',               &
      &                              'TeH:  C3 grass',                              &
      &                              'TrH:  C4 grass' ]
    ! Index of PFT (columns) in the lctlib file
    !   in the module mo_jsb_lctlib_class the 'SELECT CASE (model_usecase)' takes care that the lctlib columns 1:8 are used
    INTEGER, PARAMETER            :: pft_ids(npft) = (/ 1,2,3,4,5,6,7,8 /)
    CHARACTER(len=*), PARAMETER   :: routine = modname//':init_usecase_quincy_eight_pfts'

    ! overview of tile structure:
    !
    !  box
    !   |
    !   |
    !   |
    !  land
    !   |
    !   |
    !   |
    !  vegetation
    !   |
    !   |-----------------
    !   |      |      |
    !  pft01  pft02  pftXX

    ! QUINCY model must be enabled for this usecase
    SELECT CASE (model%config%model_scheme)
    CASE (MODEL_QUINCY)
      CALL message(TRIM(routine), 'Using the model_scheme MODEL_QUINCY')
    CASE DEFAULT
      CALL finish(TRIM(routine), 'The quincy_eight_pfts model usecase does work only with the model_scheme MODEL_QUINCY.')
    END SELECT

    !------------------------------------------------------------------------------------------------------ !
    ! box tile
    ! 
    ! Create root of the tile structure representing the whole grid box
    box_tile => t_jsb_tile('box', 'Top tile', &
      & processes      =(/A2L_,     L2A_,     RAD_     , TURB_    /), &
      & process_actions=(/ON_TILE_, ON_TILE_, ON_LEAFS_, ON_LEAFS_/), &
      & model_id=model%id, fract_varname='notsea')  ! fraction: all that "is not sea"
    model%top => box_tile

    !------------------------------------------------------------------------------------------------------ !
    ! land tile w/o process
    ! 
    land_tile => t_jsb_tile('land', 'Land tile', &
      ! add processes here only to make model init and jsb4_forcing working
      ! none of the tasks of these processes is running via task queue
      & processes       = (/SEB_    /), &
      & process_actions = (/ON_TILE_/), &
      & lct=t_jsb_lct(LAND_TYPE), parent=box_tile, &
      & fract_varname='equal')   ! Fraction of land tile = 1
                                 ! "equal" specifically means, no fraction input file is read and all tiles have the same fraction

    !------------------------------------------------------------------------------------------------------ !
    ! vegetation tile with PFT tiles as leafs
    ! 
    veg_tile => t_jsb_tile('veg', 'Veg tile', &
      & processes       = (/HYDRO_,    SSE_,     ASSIMI_,   PHENO_,    VEG_,      SB_,       SPQ_/), &
      & process_actions = (/ON_LEAFS_, ON_LEAFS_,ON_LEAFS_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_, ON_LEAFS_/), &
      & lct=t_jsb_lct(VEG_TYPE), parent=land_tile,  &
      & fract_varname='fract_veg')

    ! create PFT tiles 1:npft
    DO i = 1,npft
      ! via contructor of t_jsb_tile
      pft_tile => t_jsb_tile(pft_shortnames(i), &
       &                     TRIM(pft_longnames(i)), &
       &                     lct = t_jsb_lct(VEG_TYPE, lib_id = pft_ids(i)), &
       &                     parent = veg_tile)
    END DO

    ! catch compiler warnings and cleanup pointers
    !
    ! ! Next 3 lines only to avoid compiler warnings about pointers not being dereferenced
    ! IF (ASSOCIATED(glacier_tile)) CONTINUE
    ! IF (ASSOCIATED(lake_tile)) CONTINUE
    ! IF (ASSOCIATED(pft_tile)) CONTINUE
    NULLIFY(box_tile, land_tile, veg_tile, pft_tile)
    IF (model%config%use_glacier) NULLIFY(glacier_tile)
    IF (model%config%use_lakes) NULLIFY(lake_tile)

  END SUBROUTINE init_usecase_quincy_eight_pfts

#endif
END MODULE mo_jsb_model_usecases
