MODULE mo_moffledac_iface
#ifndef __NO_JSBACH__

  USE mo_jsb_model_class,  ONLY: t_jsb_model
  USE mo_kind,             ONLY: wp
  USE mo_exception,        ONLY: finish

  IMPLICIT NONE
  PRIVATE

  INTEGER, SAVE :: jsbach_submodel_type = -1
  INTEGER, SAVE, POINTER :: timer_moffledac(:) => NULL()
  INTEGER, SAVE :: nmoffledac_models = 1

  PUBLIC :: init_moffledac, init_moffledac_grids, get_jsb_var_name, get_var_name_add, timer_moffledac, var_has_fill

  CHARACTER(len=*), PARAMETER :: modname = 'mo_moffledac_iface'

CONTAINS

  SUBROUTINE init_moffledac(model,model_id)
#ifdef HAVE_MOFFLEDAC
    USE mo_jsb_time_iface,     ONLY: get_date_components, get_time_dt, get_time_start, get_time_experiment_start, get_year_length
!    USE mo_jsb_parallel_iface, ONLY: get_workgroup_parallel_setup
    USE mo_exception,          ONLY: finish
    USE mo_input_calendar,     ONLY: INPUT_CALENDAR_DAY360, INPUT_CALENDAR_DAY365, INPUT_CALENDAR_DAY366, &
                                     INPUT_CALENDAR_JULIAN, INPUT_CALENDAR_GREGORIAN, INPUT_CALENDAR_UNDEFINED, &
                                     TICKS_PER_SEC
    USE mo_input_arrays,       ONLY: ArrayAddData
    USE mo_input_strings,      ONLY: Val2Str
    USE mo_input,              ONLY: InputSubModelTypeGetNew, InputSubModelGetCurr, InputInit, InputCreate, InputOptImport
    USE mo_timer,              ONLY: new_timer
#ifdef __ICON__
    USE mo_io_units,           ONLY: nnml, find_next_free_unit
    USE mo_time_config,        ONLY: time_config
    USE mo_mpi,                ONLY: get_my_mpi_work_communicator, get_my_mpi_work_comm_size, get_my_mpi_work_id
#else
    USE mo_filename,           ONLY: find_next_free_unit
    USE mo_time_control,       ONLY: delta_time
    USE mo_io,                 ONLY: IO_timestep
    USE mo_mpi,                ONLY: p_pe, p_nprocs, p_all_comm
    INTEGER :: nnml = 10
#endif
#endif
    TYPE (t_jsb_model),   POINTER              :: model
    INTEGER,              INTENT(in)           :: model_id
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_moffledac'
#ifdef HAVE_MOFFLEDAC
    INTEGER  :: npe, mype, iope, comm               ! Parameters for parallel setup
    INTEGER  :: exp_start(4), run_start(4), dt, cal ! Parameters for time setup
    INTEGER  :: tmp(6), nday, i
    REAL(wp) :: rtmp

    ! Get time parameters
#ifdef __ICON__
    rtmp = time_config%dt_dyn_sec(model_id) !get_time_dt(model_id)
#else
   IF (delta_time /= 0._wp) THEN
      rtmp = delta_time
    ELSE
      rtmp = REAL(IO_timestep,wp)
    ENDIF
#endif
    IF (rtmp*TICKS_PER_SEC /= FLOOR(rtmp*TICKS_PER_SEC)) CALL finish(routine,'Sorry - timesteps must be an integer number of ticks')
    dt = INT(rtmp)
    CALL get_date_components(get_time_experiment_start(),tmp(1),tmp(2),tmp(3),tmp(4),tmp(5),tmp(6))
    exp_start(1:3) = tmp(1:3)
    exp_start(4)   = (tmp(4)*3600+tmp(5)*60+tmp(6)) * TICKS_PER_SEC
    CALL get_date_components(get_time_start(),tmp(1),tmp(2),tmp(3),tmp(4),tmp(5),tmp(6))
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
!    CALL get_workgroup_parallel_setup(comm,npe,mype,iope)
#ifdef __ICON__
    comm = get_my_mpi_work_communicator()
    npe  = get_my_mpi_work_comm_size()
    mype = get_my_mpi_work_id()
#else
    comm = p_all_comm
    npe  = p_nprocs
    mype = p_pe
#endif
    iope = 0 ! StW: TODO: This shouldn't be hardcoded?!
  !------------------------------------------------------------------------------

    ! Global MofFlEDAc initialization
    i = nnml
    CALL InputCreate(exp_start(1), exp_start(2), exp_start(3), exp_start(4), calendar=cal,                     &
                     nml_setup=TRIM(model%namelist_filename), unt=i, funt=find_next_free_unit(10,1000), &
                     npe=npe, pe=mype, iope=iope, comm=comm)

    ! Submodel type
    IF (jsbach_submodel_type < 0) jsbach_submodel_type = InputSubModelTypeGetNew()

    ! Initialize new MofFlEDAc sub-model
    CALL InputInit(jsbach_submodel_type, dt, run_start(1), run_start(2), run_start(3), run_start(4), &
                   nml_file=TRIM(model%namelist_filename), unt = i)

    ! MofFlEDAc uses 0-based model numbers, JSBACH 1-based
    IF (InputSubModelGetCurr() /= model_id-1) CALL finish(routine,'Mismatch between general model id and moffledac''s submodel id')

    ! Add new timer
    timer_moffledac => ArrayAddData(timer_moffledac,(/new_timer('moffledac_'//TRIM(Val2Str(model_id)))/),nmoffledac_models,3)

    ! Get all input_ctl namelists regarding this model
    CALL InputOptImport(TRIM(model%namelist_filename))

    ! Reading of some further namelists is demanded to start from a free, unopened file unit
    close(i)

#else
    ! Just preventing warnings about unused dummy arguments
#endif

  END SUBROUTINE init_moffledac

  SUBROUTINE init_moffledac_grids(model,patch,decomp)
#ifdef __ICON__
    USE mo_parallel_config,    ONLY: nproma
    USE mo_input,              ONLY: InputFileAdd, InputDimInq, input_file_list
#else
    USE mo_decomposition,      ONLY: pe_decomposed
    USE mo_input,              ONLY: InputDimLocalSetMulti, input_dim_list
#endif
    USE mo_jsb_domain_iface,   ONLY: t_patch
    USE mo_jsb_vertical_axes,  ONLY: jsbach_vgrids
    USE mo_exception,          ONLY: finish
    USE mo_input,              ONLY: InputDimLocalSet, InputDimEquivalence, InputDimAdd
#ifdef __ICON__
    INTEGER,              INTENT(in), OPTIONAL :: decomp ! Dummy placeholder for ECHAM compsition
#else
    TYPE (pe_decomposed), INTENT(in), OPTIONAL :: decomp
#endif
    TYPE (t_patch),       INTENT(in), OPTIONAL :: patch
    TYPE (t_jsb_model),   POINTER              :: model
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_moffledac_grids'

#ifdef HAVE_MOFFLEDAC

    INTEGER  :: nday, i, g, l , nblock
#ifdef __ICON__
    INTEGER :: n
#else
    TYPE (input_dim_list),  POINTER :: dim_par
    LOGICAL, ALLOCATABLE :: lmask(:), lAtPE(:,:)
    INTEGER :: n_hemisphere_lat, lons, lone, lonsz, col_blocks, col_blocks_per_pe, baseLo, baseHi
#endif

    ! Add time dimension
    CALL InputDimAdd('time',lRecSep=.TRUE.)

    ! Setup the horizontal grid and spread it correctly across PEs
    IF (PRESENT(patch) .EQV. PRESENT(decomp)) CALL finish(routine,'Called with insufficient/inconsistent parameters')
    IF (PRESENT(patch)) THEN
      ! ICON setup
#ifdef __ICON__
      ! Get land sea mask
      CALL InputDimAdd('cell',alt_name='ncells')
      CALL InputFileAdd(TRIM(model%config%fract_filename),dt_file=0)
      ! Obtain world size (size of ncells dimension)
      CALL InputDimInq('cell',global=g,local=l)
      IF (g /= l) CALL finish(TRIM(routine),'cell dimension not parallelized. Should thus have same local and global size')
      n = SIZE(patch%cells%decomp_info%glb_index)
      DO WHILE (patch%cells%decomp_info%glb_index(n) == 0)
        n = n - 1
      ENDDO
      CALL InputDimLocalSet(dim_name='cell',idx=patch%cells%decomp_info%glb_index(1:n))
      nblock = (SIZE(patch%cells%decomp_info%glb_index)+nproma-1) / nproma
      CALL InputDimAdd('nproma',extent=nproma)
      CALL InputDimAdd('nblock',extent=nblock,global=-1);
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
#else
    ! Just preventing warnings about unused dummy arguments
#endif

  END SUBROUTINE init_moffledac_grids

  FUNCTION get_jsb_var_name(varname,base,prefix,suffix)
  CHARACTER (len=64) :: get_jsb_var_name
  CHARACTER (len=*), INTENT(in) :: varname, base, prefix, suffix

  CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_jsb_var_name'
  CHARACTER (len=64) :: ret
  INTEGER :: st, en

    ret = TRIM(varname)
    st = 1
    en = LEN_TRIM(ret)
    IF (ret == '') THEN
      ret = TRIM(prefix)//'_'//TRIM(base)//'_'//TRIM(suffix)
      en = LEN_TRIM(ret)
      IF (prefix == '') st = 2
      IF (suffix == '') en = en - 1
    ENDIF

    get_jsb_var_name = ret(st:en)

  END FUNCTION get_jsb_var_name

  FUNCTION get_var_name_add(prefix)
  CHARACTER (len=4) :: get_var_name_add
  CHARACTER (len=*) :: prefix

  CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_var_name_add'

    IF ((prefix == '') .OR. (prefix(1:1) == '_')) THEN
      get_var_name_add = 'loc'
    ELSE
      get_var_name_add = 'loc_'
    ENDIF

  END FUNCTION get_var_name_add

  LOGICAL FUNCTION var_has_fill(var,val)
#ifdef HAVE_MOFFLEDAC
  USE mo_input, ONLY: input_var_list, InputDataGetSpecified, INPUT_DATA_FILL,INPUT_DATA_NO_STEP
  TYPE (input_var_list), POINTER  :: var
#else
  INTEGER, POINTER :: var
#endif
  REAL(wp), INTENT(out), OPTIONAL :: val

  REAL(wp), POINTER :: fill_val(:,:,:,:,:,:,:)

#ifdef HAVE_MOFFLEDAC
    fill_val => InputDataGetSpecified(var%aux_data,INPUT_DATA_FILL,INPUT_DATA_NO_STEP)
#else
    NULLIFY(fill_val)
#endif
    IF (ASSOCIATED(fill_val)) THEN
      IF (PRESENT(val)) val = fill_val(1,1,1,1,1,1,1)
      var_has_fill = .TRUE.
    ELSE
      var_has_fill = .FALSE.
      IF (PRESENT(val)) val = 0._wp
    ENDIF

  END FUNCTION var_has_fill

#endif

END MODULE mo_moffledac_iface
