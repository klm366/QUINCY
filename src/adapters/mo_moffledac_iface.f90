MODULE mo_moffledac_iface
#ifndef __NO_JSBACH__

  USE mo_jsb_model_class,  ONLY: t_jsb_model
  USE mo_kind,             ONLY: wp
  USE mo_exception,        ONLY: finish

  IMPLICIT NONE
  PRIVATE

  INTEGER, SAVE :: jsbach_submodel_type = -1

  PUBLIC :: init_moffledac

  CHARACTER(len=*), PARAMETER :: modname = 'mo_moffledac_iface'

CONTAINS

  SUBROUTINE init_moffledac(model,model_id,patch,decomp)
    USE mo_jsb_domain_iface,   ONLY: pe_decomposed, t_patch
#ifdef HAVE_MOFFLEDAC
    USE mo_jsb_time_iface,     ONLY: get_date_components, get_time_dt, get_time_start, get_time_current, get_year_length
    USE mo_jsb_parallel_iface, ONLY: get_workgroup_parallel_setup
    USE mo_jsb_vertical_axes,  ONLY: jsbach_vgrids
    USE mo_exception,          ONLY: message_text, message, finish
    USE mo_input_calendar,     ONLY: INPUT_CALENDAR_DAY360, INPUT_CALENDAR_DAY365, INPUT_CALENDAR_DAY366, &
                                     INPUT_CALENDAR_JULIAN, INPUT_CALENDAR_GREGORIAN, INPUT_CALENDAR_UNDEFINED, &
                                     TICKS_PER_SEC
    USE mo_input,              ONLY: InputSubModelTypeGetNew, InputSubModelGetCurr, InputInit, InputCreate, InputOptImport, &
                                     InputDimLocalSet, InputDimEquivalence, InputDimAdd
#ifdef __ICON__
    USE mo_io_units,           ONLY: nnml, find_next_free_unit
    USE mo_parallel_config,    ONLY: nproma
    USE mo_input,              ONLY: InputFileAdd, InputDimInq, input_file_list
#else
    USE mo_input,              ONLY: InputDimLocalSetMulti, input_dim_list
    USE mo_filename,           ONLY: find_next_free_unit
    INTEGER :: nnml = 10
#endif
#endif
    TYPE (t_jsb_model),   POINTER              :: model
    INTEGER,              INTENT(in)           :: model_id
    TYPE (t_patch),       INTENT(in), OPTIONAL :: patch
    TYPE (pe_decomposed), INTENT(in), OPTIONAL :: decomp

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_moffledac'
#ifdef HAVE_MOFFLEDAC
    INTEGER  :: npe, mype, iope, comm               ! Parameters for parallel setup
    INTEGER  :: exp_start(4), run_start(4), dt, cal ! Parameters for time setup
    INTEGER  :: tmp(6), nblock
    INTEGER  :: nday, i, g, l
    REAL(wp) :: rtmp
#ifdef __ICON__
    TYPE (input_file_list), POINTER :: IFile
    REAL(wp), POINTER :: slm(:)
    INTEGER :: n
#else
    TYPE (input_dim_list),  POINTER :: dim_par
    LOGICAL, ALLOCATABLE :: lmask(:), lAtPE(:,:)
    INTEGER :: n_hemisphere_lat, lons, lone, lonsz, col_blocks, col_blocks_per_pe, baseLo, baseHi
#endif

    ! Get time parameters
    rtmp = get_time_dt()
    IF (rtmp*TICKS_PER_SEC /= FLOOR(rtmp*TICKS_PER_SEC)) CALL finish(routine,'Sorry - timesteps must be an integer number of ticks')
    dt = INT(rtmp)
    CALL get_date_components(get_time_start(),tmp(1),tmp(2),tmp(3),tmp(4),tmp(5),tmp(6)) ! Is get_time_start really start of exp?
    exp_start(1:3) = tmp(1:3)
    exp_start(4)   = (tmp(4)*3600+tmp(5)*60+tmp(6)) * TICKS_PER_SEC
    CALL get_date_components(get_time_current(),tmp(1),tmp(2),tmp(3),tmp(4),tmp(5),tmp(6))
    run_start(1:3) = tmp(1:3)
    run_start(4)   = (tmp(4)*3600+tmp(5)*60+tmp(6)) * TICKS_PER_SEC
    nday = get_year_length(1900) + get_year_length(1901) + get_year_length(2000)
    SELECT CASE (nday)
      CASE (1080)
        cal = INPUT_CALENDAR_DAY360
      CASE (1095)
        cal = INPUT_CALENDAR_DAY365
      CASE (1098)
        cal = INPUT_CALENDAR_DAY366
      CASE (1097)
        cal = INPUT_CALENDAR_JULIAN
      CASE (1096)
        cal = INPUT_CALENDAR_GREGORIAN
      CASE DEFAULT
        cal = INPUT_CALENDAR_UNDEFINED
    END SELECT

    ! Get parallel parameters
    CALL get_workgroup_parallel_setup(comm,npe,mype,iope)

    ! Global MofFlEDAc initialization
    i = nnml
    CALL InputCreate(exp_start(1), exp_start(2), exp_start(3), exp_start(4), calendar=cal,                     &
                     nml_setup=TRIM(model%config%namelist_filename), unt=i, funt=find_next_free_unit(10,1000), &
                     npe=npe, pe=mype, iope=iope, comm=comm)

    ! Submodel type
    IF (jsbach_submodel_type < 0) jsbach_submodel_type = InputSubModelTypeGetNew()

    ! Initialize new MofFlEDAc sub-model
    CALL InputInit(jsbach_submodel_type, dt, run_start(1), run_start(2), run_start(3), run_start(4), &
                   nml_file=TRIM(model%config%namelist_filename), unt = i)

    ! MofFlEDAc uses 0-based model numbers, JSBACH 1-based
    IF (InputSubModelGetCurr() /= model_id-1) CALL finish(routine,'Mismatch between general model id and moffledac''s submodel id')

    ! Get all input_ctl namelists regarding this model
    CALL InputOptImport(TRIM(model%config%namelist_filename))

    ! Add time dimension
    CALL InputDimAdd('time',lRecSep=.TRUE.)

    ! Setup the horizontal grid and spread it correctly across PEs
    IF (PRESENT(patch) .EQV. PRESENT(decomp)) CALL finish(routine,'Called with insufficient/inconsistent parameters')
    IF (PRESENT(patch)) THEN
      ! ICON setup
#ifdef __ICON__
      ! Get land sea mask
      CALL InputDimAdd('cell',alt_name='ncells')
      CALL InputFileAdd(TRIM(model%config%fract_filename),IFile,dt_file=0)
      ! Obtain world size (size of ncells dimension)
      CALL InputDimInq('cell',global=g,local=l)
      IF (g /= l) CALL finish(TRIM(routine),'cell dimension not parallelized. Should thus have same local and global size')
!      CALL InputVarAdd('cell',slm,'notsea',IFile,dt_update=0)
!      CALL InputGetData()
      n = SIZE(patch%cells%decomp_info%glb_index)
      DO WHILE (patch%cells%decomp_info%glb_index(n) == 0)
        n = n - 1
      ENDDO
      CALL InputDimLocalSet(dim_name='cell',idx=patch%cells%decomp_info%glb_index(1:n))
      nblock = (SIZE(patch%cells%decomp_info%glb_index)+nproma-1) / nproma
      CALL InputDimAdd('nproma',extent=nproma)
      CALL InputDimAdd('nblock',extent=nblock);
      CALL InputDimEquivalence('nproma,nblock','cell')
#else
      CALL finish(routine,'Inconsistent model setup')
#endif
    ELSE
#ifdef __ICON__
      CALL finish(routine,'Inconsistent model setup')
#else
      ! ECHAM setup
      CALL InputDimAdd('lon',alt_name='X',extent=decomp%nlon,lcyclic=.TRUE.)
      CALL InputDimAdd('lat',alt_name='Y',extent=decomp%nlat,ldec=.TRUE.)

      CALL InputDimAdd('nproma',extent=decomp%nproma)

      n_hemisphere_lat = INT(CEILING(REAL(decomp%nlat,wp)/REAL(decomp%nproca*2,wp)))
      lonsz            = INT(CEILING(REAL(decomp%nlon,wp)/REAL(decomp%nprocb  ,wp)))
      IF ((decomp%nlon/decomp%nprocb/decomp%nproma)*decomp%nprocb*decomp%nproma == decomp%nlon) THEN 
        ! Case is also covered by the general case below, but causes less overhead
        nblock = decomp%nlon*decomp%nlat/decomp%nproma
        CALL InputDimAdd('nblock',Extent=nblock,dim_ref=dim_par)
        CALL InputDimEquivalence('nproma,nblock','lon,lat')
        col_blocks         = decomp%nlon / decomp%nproma
        col_blocks_per_pe  = decomp%nlon / decomp%nproma / decomp%nprocb
        baseLo = col_blocks_per_pe*(decomp%set_b-1) + col_blocks * n_hemisphere_lat * (decomp%set_a-1)
        baseHi = col_blocks_per_pe*(decomp%set_b-1) - col_blocks * n_hemisphere_lat * (decomp%set_a-1) + nblock
        ALLOCATE(lmask(nblock))
        lmask(:) = .FALSE.
        DO i=1,n_hemisphere_lat
          baseHi = baseHi - col_blocks
          lmask(baseLo+1:baseLo+col_blocks_per_pe) = .TRUE.
          lmask(baseHi+1:baseHi+col_blocks_per_pe) = .TRUE.
          baseLo = baseLo + col_blocks
        ENDDO
        CALL InputDimLocalSet(dim_ref=dim_par,lmask=lmask)
        DEALLOCATE(lmask)
!        CALL InputDimPackedAdd('landpoints','nproma,nblock',lmask=slm>0.5_wp)
      ELSE
        baseLo =               n_hemisphere_lat *(decomp%set_a-1)
        baseHi = decomp%nlat - n_hemisphere_lat * decomp%set_a
        lons   = lonsz * (decomp%set_b-1) + 1
        lone   = MIN(lons+lonsz-1,decomp%nlon)
        ALLOCATE(lAtPE(decomp%nlon,decomp%nlat))
        lAtPE = .FALSE.
        lAtPE(lons:lone,baseLo+1:baseLo+n_hemisphere_lat) = .TRUE.
        lAtPE(lons:lone,baseHi+1:baseHi+n_hemisphere_lat) = .TRUE.
        CALL InputDimAdd('lon_iloc',Extent=lonsz) ! Just to get proper sizes of internal dimensions
        CALL InputDimLocalSetMulti('lon,lat',lmask=lAtPE)
        DEALLOCATE(lAtPE)
!        CALL InputDimPackedAdd('landpoints','lon_iloc,lat_iloc',lmask=slm>0.5_wp)
        ! Preparation of reading of ECHAM grid-point variables.
        CALL InputDimAdd('nblock')
        CALL InputDimEquivalence('nproma,nblock','lon_iloc,lat_iloc')
      ENDIF
#endif
    ENDIF

    ! Register all vertical axes known to JSBACH also for MofFlEDAc
    DO i=1,SIZE(jsbach_vgrids)
      CALL InputDimAdd(TRIM(jsbach_vgrids(i)%p%name),Extent=jsbach_vgrids(i)%p%n_levels)
    ENDDO
#endif

  END SUBROUTINE init_moffledac

#endif

END MODULE mo_moffledac_iface
