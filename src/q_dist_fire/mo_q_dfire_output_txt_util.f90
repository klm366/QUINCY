!> Extended functions to complement QUINCY model output as ascii text for dfire
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
MODULE mo_q_dfire_output_txt_util

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message_text, message, finish

!  USE mo_qs_output_txt_util

  USE mo_jsb_process_class,   ONLY: SB_, Q_DFIRE_, Q_DFIRE_FUEL_

  USE mo_linked_list,         ONLY: t_q_linkedlist, t_q_linkedlist_node

  dsl4jsb_Use_config(SB_)

  IMPLICIT NONE

  PUBLIC  ! all is public

  INTEGER, PARAMETER, PRIVATE   :: ic = 1                             !< the point in ics:ice of which data will be written to file; per default = 1

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_q_dfire_output_util'

  ENUM, BIND(C) ! outputs operations
    ENUMERATOR :: EN_OP_AVERAGE = 1
  END ENUM

  ENUM, BIND(C) ! preprocessing array operations
    ENUMERATOR :: EN_PRE_OP_NULL = 1, &
                  EN_PRE_OP_SUM
  END ENUM

  !> -Use for scalar values
  !> -Use for array values, those are then SUM/... according to pre_operation
  TYPE :: t_qs_output_2d !< derived type to handle output operations, such as average
    CHARACTER(:), ALLOCATABLE :: name, units
    REAL(wp), POINTER :: data          !< Pointer to data
    REAL(wp), POINTER :: data_r1d(:)   !< pointer to array, preprocessing is required to deal with array -> to scalar value
    REAL(wp)          :: dataccum      !<
    REAL(wp)          :: counter       !< counter how many TS
    INTEGER           :: pre_operation !< options: SUM
    INTEGER           :: operation     !< average, inst, instantaneous(1st vaule of TS)
  CONTAINS
    PROCEDURE :: runOperation => o2d_runOperation
    PROCEDURE :: get => o2d_get
    PROCEDURE :: reset => o2d_reset
    FINAL :: finalize_qs_output_2d
  END TYPE t_qs_output_2d

  INTERFACE t_qs_output_2d
    MODULE PROCEDURE :: create_qs_output_2d, create_qs_output_2d_as_array
  END INTERFACE t_qs_output_2d


  TYPE :: t_qs_output_3d !< derived type to handle output operations, such as average
    CHARACTER(:), ALLOCATABLE :: name, units
    REAL(wp), POINTER     :: data(:)       !< Pointer to data
    REAL(wp), ALLOCATABLE :: dataccum(:)   !<
    REAL(wp)              :: counter       !< counter how many TS
    INTEGER               :: operation     !< average, inst, instantaneous(1st vaule of TS)
  CONTAINS
    PROCEDURE :: runOperation => o3d_runOperation
    PROCEDURE :: get => o3d_get
    PROCEDURE :: reset => o3d_reset
    FINAL :: finalize_qs_output_3d
  END TYPE t_qs_output_3d

  INTERFACE t_qs_output_3d
    MODULE PROCEDURE :: create_qs_output_3d
  END INTERFACE t_qs_output_3d

  REAL(wp), PARAMETER :: units_mols_to_micromols = 1000000._wp ! todo: move to dfire_constants

  INTERFACE add_to_list
    MODULE PROCEDURE :: add_to_list_o2d, add_to_list_o3d, add_to_list_of
  END INTERFACE add_to_list

  integer, save :: uniqueid = 1000 ! Generate unique file id

  TYPE :: t_qs_output_file
    integer :: fid
    character(:), allocatable :: filename
    type(t_q_linkedlist) :: vars
  CONTAINS
    procedure :: open
    procedure :: init
    generic, public :: add_var => add_var_by_object, add_var_by_name
    procedure :: add_var_by_object
    procedure :: add_var_by_name
    procedure :: write
    procedure :: update
    procedure :: reset
    procedure :: close
  END TYPE t_qs_output_file

  INTERFACE t_qs_output_file
    module procedure :: create_qs_output_file
  END INTERFACE t_qs_output_file


  !< Main interface to interact with Simple Outputs library
  TYPE :: t_qs_soutputs
    type(t_q_linkedlist) :: files
  CONTAINS
    !< user subroutines
    procedure :: add_file !< return file id
    procedure :: add_var => addvar_soutputs
!    !< Infraestructure subroutines
    procedure :: open => open_soutputs !< create file(s)
    procedure :: init => init_soutputs !< add headers into file
    procedure :: update => update_soutputs !< process with new data (if required, e.g mean)
    procedure :: write => write_soutputs !< write data into file
    procedure :: reset => reset_soutputs !< reset internal counters (if required, e.g mean counter)
    procedure :: close => close_soutputs !< close all files
  END TYPE t_qs_soutputs


  INTERFACE t_qs_soutputs
    module procedure :: create_qs_soutputs
  END INTERFACE t_qs_soutputs

CONTAINS

  !>
  !> Generete unique number (integer)
  !>
  INTEGER FUNCTION gen_unique_id()
    gen_unique_id = uniqueid
    uniqueid = uniqueid + 1
  END FUNCTION gen_unique_id


  ! ======================================================================================================= !
  !>
  !> Constructor for t_qs_soutputs
  !>
  FUNCTION create_qs_soutputs() RESULT (sio)
    TYPE(t_qs_soutputs) :: sio
    sio % files = t_q_linkedlist()
  END FUNCTION create_qs_soutputs


  ! ======================================================================================================= !
  !>
  !> Constructor for t_qs_output_file
  !>
  FUNCTION create_qs_output_file(filename) RESULT (ofile)
    CHARACTER(*), INTENT(in) :: filename
    TYPE(t_qs_output_file) :: ofile

    ofile % filename = filename
    ofile % fid = gen_unique_id()
    ofile % vars = t_q_linkedlist()
  END FUNCTION create_qs_output_file

  ! ======================================================================================================= !
  !>
  !> Add new variable into the file
  !>
  SUBROUTINE add_var_by_object(this, var2d, var3d)
    CLASS(t_qs_output_file), INTENT(inout) :: this
    TYPE(t_qs_output_2d), INTENT(in), OPTIONAL :: var2d
    TYPE(t_qs_output_3d), INTENT(in), OPTIONAL :: var3d

    IF (PRESENT(var2d)) THEN
      CALL add_to_list(this%vars, var2d)
    ELSEIF (PRESENT(var3d)) THEN
      CALL add_to_list(this%vars, var3d)
    ELSE
      CALL finish("add_var", "No 2d or 3d variable defined as an argument")
    END IF
  END SUBROUTINE add_var_by_object

  ! ======================================================================================================= !
  !>
  !> Directly add a new variable into the file
  !>
  SUBROUTINE add_var_by_name(this, name, units, pdata2d, pdata3d, preop, op)
    CLASS(t_qs_output_file), INTENT(inout) :: this
    CHARACTER(*), INTENT(in) :: name
    CHARACTER(*), INTENT(in) :: units
    INTEGER, INTENT(in) :: op
    INTEGER, INTENT(in), OPTIONAL :: preop
    REAL(wp), INTENT(in), POINTER, OPTIONAL :: pdata2d
    REAL(wp), INTENT(in), POINTER, OPTIONAL :: pdata3d(:)

    TYPE(t_qs_output_3d), POINTER :: newvar_o3d
    TYPE(t_qs_output_2d), POINTER :: newvar_o2d

    IF (PRESENT(pdata2d) .OR.  (PRESENT(pdata3d) .AND. PRESENT(pdata3d)) ) THEN
      ALLOCATE(newvar_o2d)
      IF (PRESENT(preop)) THEN
        ! When applying a preoperation, it must be given pdata3d instead of pdata2d
        newvar_o2d = t_qs_output_2d(name, units, pdata3d, preop, op)
      ELSE
        newvar_o2d = t_qs_output_2d(name, units, pdata2d, op)
      END IF
      CALL this % add_var_by_object(var2d=newvar_o2d)
    ELSE IF (PRESENT(pdata3d)) THEN
      ALLOCATE(newvar_o3d)
      newvar_o3d = t_qs_output_3d(name, units, pdata3d, op)
      CALL this % add_var_by_object(var3d=newvar_o3d)
    ELSE
      CALL finish("add_var_ny_name", "No pdata2d nor pdata3d found")
    END IF
  END SUBROUTINE add_var_by_name

  ! ======================================================================================================= !
  !>
  !> Initialize file (header)
  !>
  SUBROUTINE init(this)
    CLASS(t_qs_output_file), INTENT(inout) :: this
    CHARACTER(:), ALLOCATABLE :: header

    header = "Time"
    CALL this%vars%traverse(getHeaders)
  CONTAINS
    !> Extract var name
    SUBROUTINE getHeaders(node)
      TYPE(t_q_linkedlist_node), POINTER, INTENT(inout)  :: node

      SELECT TYPE(v=>node%value)
      CLASS IS (t_qs_output_2d)
        header = header // " " // v % name
      CLASS IS (t_qs_output_3d)
        header = header // " " // v % name
      CLASS DEFAULT
        ! unexpected type found, raise error
        CALL finish("getHeaders", "Unexpected class found")
      END SELECT
    END SUBROUTINE getHeaders
  END SUBROUTINE init


  ! ======================================================================================================= !
  !>
  !> Open file
  !>
  SUBROUTINE open(this)
    CLASS(t_qs_output_file), INTENT(in) :: this
    OPEN (UNIT=this%fid, FILE=this%filename, ACTION="write", STATUS="unknown")
  END SUBROUTINE open

  ! ======================================================================================================= !
  !>
  !> Write current data into txt file
  !>
  SUBROUTINE write(this, time_d)
    CLASS(t_qs_output_file), INTENT(inout) :: this
    REAL(wp), INTENT(in) :: time_d
    REAL(wp), ALLOCATABLE :: values(:)
    INTEGER :: nsz
    INTEGER :: ic

    ic = 1
    nsz = this % vars % getLength()
    ALLOCATE(values(nsz))
    CALL this % vars % traverse(getData)
  CONTAINS
    !> Extract var name
    SUBROUTINE getData(node)
      TYPE(t_q_linkedlist_node), POINTER, INTENT(inout)  :: node

      SELECT TYPE(v=>node%value)
      CLASS IS (t_qs_output_2d)
        values(ic) = v % get()
      CLASS IS (t_qs_output_3d)
        CALL finish("getData", "t_qs_output_3d not yet implemented")
      CLASS DEFAULT
        ! unexpected type found, raise error
        CALL finish("getData", "Unexpected class found")
      END SELECT
      ic = ic + 1
    END SUBROUTINE getData
  END SUBROUTINE write

  ! ======================================================================================================= !
  !>
  !> Apply an operation (average, ...) for each variable
  !>
  SUBROUTINE update(this)
    CLASS(t_qs_output_file), INTENT(inout) :: this
    CALL this % vars % traverse(runAllOperations)
  END SUBROUTINE update

  ! ======================================================================================================= !
  !>
  !> Reset all counters for Averaging variables
  !>
  SUBROUTINE reset(this)
    CLASS(t_qs_output_file), INTENT(inout) :: this
    CALL this % vars % traverse(resetAllVars)
  CONTAINS
    !> Extract var name
    SUBROUTINE resetAllVars(node)
      TYPE(t_q_linkedlist_node), POINTER, INTENT(inout)  :: node

      SELECT TYPE(v=>node%value)
      CLASS IS (t_qs_output_2d)
        CALL v % reset()
      CLASS IS (t_qs_output_3d)
        CALL v % reset()
      CLASS DEFAULT
        ! unexpected type found, raise error
        CALL finish("resetAllVars", "Unexpected class found")
      END SELECT
    END SUBROUTINE resetAllVars
  END SUBROUTINE reset


  ! ======================================================================================================= !
  !>
  !> Close file
  !>
  SUBROUTINE close(this)
    CLASS(t_qs_output_file), INTENT(in) :: this
    logical :: opened
    INQUIRE(this%fid, opened=opened)
    IF (opened) CLOSE(this%fid)
  END SUBROUTINE close

  ! ======================================================================================================= !
  !>
  !> Constructor for t_qs_output_3d
  !>
  FUNCTION create_qs_output_3d(name, units, pdata, op) RESULT (var3d)
    CHARACTER(*), INTENT(in) :: name
    CHARACTER(*), INTENT(in) :: units
    REAL(wp), INTENT(in), POINTER :: pdata(:)
    INTEGER, INTENT(in) :: op
    TYPE(t_qs_output_3d) :: var3d

    var3d % name = name
    var3d % units = units
    var3d % data => pdata
    var3d % operation = op
    ALLOCATE(var3d % dataccum(SIZE(var3d % data)))
    var3d % dataccum(:) = 0._wp
    var3d % counter = 0._wp
  END FUNCTION create_qs_output_3d


  ! ======================================================================================================= !
  !>
  !>
  !>
  SUBROUTINE finalize_qs_output_3d(this)
    TYPE(t_qs_output_3d), INTENT(inout) :: this

    IF (ALLOCATED(this % name)) DEALLOCATE(this % name)
    IF (ALLOCATED(this % units)) DEALLOCATE(this % units)
    NULLIFY(this % data)
  END SUBROUTINE finalize_qs_output_3d


  ! ======================================================================================================= !
  !>
  !> Constructor for t_qs_output_2d
  !>
  FUNCTION create_qs_output_2d(name, units, pdata, op) RESULT (var2d)
    CHARACTER(*), INTENT(in) :: name
    CHARACTER(*), INTENT(in) :: units
    REAL(wp), INTENT(in), POINTER :: pdata
    INTEGER, INTENT(in) :: op
    TYPE(t_qs_output_2d) :: var2d

    var2d % name = name
    var2d % units = units
    var2d % data => pdata
    var2d % data_r1d => NULL()
    var2d % pre_operation = EN_PRE_OP_NULL ! only enabled when data_r1d is defined
    var2d % operation = op
    var2d % dataccum = 0._wp
    var2d % counter = 0._wp
  END FUNCTION create_qs_output_2d


  ! ======================================================================================================= !
  !>
  !> Constructor for t_qs_output_2d
  !>
  FUNCTION create_qs_output_2d_as_array(name, units, parray, preop, op) RESULT (var2d)
    CHARACTER(*), INTENT(in) :: name
    CHARACTER(*), INTENT(in) :: units
    REAL(wp), INTENT(in), POINTER :: parray(:)
    INTEGER, INTENT(in) :: preop !< Apply operation into the array (preprocessing)
    INTEGER, INTENT(in) :: op
    TYPE(t_qs_output_2d) :: var2d

    var2d % name = name
    var2d % units = units
    ALLOCATE(var2d % data)
    var2d % data = 0._wp !< place sum data here
    var2d % data_r1d => parray
    var2d % operation = op
    var2d % pre_operation = EN_PRE_OP_SUM
    var2d % dataccum = 0._wp
    var2d % counter = 0._wp
  END FUNCTION create_qs_output_2d_as_array


  ! ======================================================================================================= !
  !>
  !>
  !>
  SUBROUTINE finalize_qs_output_2d(this)
    TYPE(t_qs_output_2d), INTENT(inout) :: this

    IF (ALLOCATED(this % name)) DEALLOCATE(this % name)
    IF (ALLOCATED(this % units)) DEALLOCATE(this % units)

    IF (ASSOCIATED(this % data_r1d)) THEN
      DEALLOCATE(this % data)
      NULLIFY(this % data_r1d)
    END IF

    NULLIFY(this % data)
  END SUBROUTINE finalize_qs_output_2d


  ! ======================================================================================================= !
  !>
  !>
  !>
  SUBROUTINE o2d_runOperation(this)
    CLASS(t_qs_output_2d), INTENT(inout) :: this

    ! preprocess array
    IF (ASSOCIATED(this % data_r1d)) THEN
      IF (this % pre_operation == EN_PRE_OP_SUM) THEN ! apply SUM operation
        this % data = SUM(this % data_r1d)
      ELSE
        CALL finish("o2d_runOperation", "Preprocessing operation not supported")
      END IF
    END IF

    ! compute
    IF (this % operation == EN_OP_AVERAGE) THEN ! accumulate data, at the end average
      this % dataccum = this % dataccum + this % data
      this % counter = this % counter + 1
    ELSE
      CALL finish("o2d_runOperation", "Operation not supported")
    END IF
  END SUBROUTINE o2d_runOperation


  ! ======================================================================================================= !
  !>
  !>
  !>
  FUNCTION o2d_get(this) RESULT (value)
    CLASS(t_qs_output_2d), INTENT(in) :: this
    REAL(wp) :: value
    IF (this % operation == EN_OP_AVERAGE) THEN
      value = this % dataccum / this % counter ! apply average
    ELSE
      CALL finish("o2d_runOperation", "Operation not supported")
    END IF
  END FUNCTION o2d_get


  ! ======================================================================================================= !
  !>
  !> Reinitialize counter to 0, useful for average
  !>
  SUBROUTINE o2d_reset(this)
    CLASS(t_qs_output_2d), INTENT(inout) :: this
    this % dataccum = 0._wp
    this % counter = 0._wp
  END SUBROUTINE o2d_reset


  ! ======================================================================================================= !
  !>
  !>
  !>
  SUBROUTINE o3d_runOperation(this)
    CLASS(t_qs_output_3d), INTENT(inout) :: this

    ! compute
    IF (this % operation == EN_OP_AVERAGE) THEN ! accumulate data, at the end average
      this % dataccum = this % dataccum + this % data
      this % counter = this % counter + 1
    ELSE
      CALL finish("o3d_runOperation", "Operation not supported")
    END IF
  END SUBROUTINE o3d_runOperation


  ! ======================================================================================================= !
  !>
  !>
  !>
  FUNCTION o3d_get(this) RESULT (value)
    CLASS(t_qs_output_3d), INTENT(in) :: this
    REAL(wp) :: value(SIZE(this%dataccum))
    IF (this % operation == EN_OP_AVERAGE) THEN
      value = this % dataccum / this % counter ! apply average
    ELSE
      CALL finish("o3d_get", "Operation not supported")
    END IF
  END FUNCTION o3d_get


  ! ======================================================================================================= !
  !>
  !> Reinitialize counter to 0, useful for average
  !>
  SUBROUTINE o3d_reset(this)
    CLASS(t_qs_output_3d), INTENT(inout) :: this
    this % dataccum = 0._wp
    this % counter = 0._wp
  END SUBROUTINE o3d_reset

  ! ======================================================================================================= !
  !>
  !> Auxiliary function to add t_qs_output_2d derive type into a list
  !>
  SUBROUTINE add_to_list_of(list, outputfile)
    TYPE(t_q_linkedlist), INTENT(inout) :: list
    TYPE(t_qs_output_file), INTENT(in), TARGET :: outputfile
    CLASS(*), POINTER :: ptr
    ptr => outputfile
    CALL list % append(ptr)
  END SUBROUTINE add_to_list_of


  ! ======================================================================================================= !
  !>
  !> Auxiliary function to add t_qs_output_2d derive type into a list
  !>
  SUBROUTINE add_to_list_o2d(list, outputvar)
    TYPE(t_q_linkedlist), INTENT(inout) :: list
    TYPE(t_qs_output_2d), INTENT(in), TARGET :: outputvar
    CLASS(*), POINTER :: ptr
    ptr => outputvar
    CALL list % append(ptr)
  END SUBROUTINE add_to_list_o2d


  ! ======================================================================================================= !
  !>
  !> Auxiliary function to add t_qs_output_2d derive type into a list
  !>
  SUBROUTINE add_to_list_o3d(list, outputvar)
    TYPE(t_q_linkedlist), INTENT(inout) :: list
    TYPE(t_qs_output_3d), INTENT(in), TARGET :: outputvar
    CLASS(*), POINTER :: ptr
    ptr => outputvar
    CALL list % append(ptr)
  END SUBROUTINE add_to_list_o3d

  ! ======================================================================================================= !
  !>
  !> Execute runOperation subroutine for each element
  !>
  SUBROUTINE runAllOperations(node)
    TYPE(t_q_linkedlist_node), POINTER, INTENT(inout)  :: node

    SELECT TYPE(v=>node%value)
    CLASS IS (t_qs_output_2d)
      CALL v % runOperation()
    CLASS IS (t_qs_output_3d)
      CALL v % runOperation()
    CLASS DEFAULT
      ! unexpected type found, raise error
      CALL finish("runAllOperations", "Unexpected class found")
    END SELECT
  END SUBROUTINE runAllOperations


  ! ======================================================================================================= !
  !>
  !> Execute reset subroutine for each element
  !>
  SUBROUTINE resetAll(node)
    TYPE(t_q_linkedlist_node), POINTER, INTENT(inout)  :: node

    SELECT TYPE(v=>node%value)
    CLASS IS (t_qs_output_2d)
      CALL v % reset()
    CLASS IS (t_qs_output_3d)
      CALL v % reset()
    CLASS DEFAULT
      ! unexpected type found, raise error
      CALL finish("resetAll", "Unexpected class found")
    END SELECT
  END SUBROUTINE resetAll


  ! ======================================================================================================= !
  !>
  !> Open files from soutputs
  !>
  SUBROUTINE open_soutputs(this)
    CLASS(t_qs_soutputs), INTENT(inout) :: this
    CALL this % files % traverse(openAllFiles) ! create files
  CONTAINS
      SUBROUTINE openAllFiles(node)
        TYPE(t_q_linkedlist_node), POINTER, INTENT(inout)  :: node

        SELECT TYPE(v=>node%value)
        CLASS IS (t_qs_output_file)
          CALL v % open()
        CLASS DEFAULT
          ! unexpected type found, raise error
          CALL finish("open_soutputs", "Unexpected class found", stacktrace=.TRUE.)
        END SELECT
      END SUBROUTINE openAllFiles
  END SUBROUTINE open_soutputs


  ! ======================================================================================================= !
  !>
  !> Initialize all files (e.g: add headers into files, ...)
  !>
  SUBROUTINE init_soutputs(this)
    CLASS(t_qs_soutputs), INTENT(inout) :: this
    CALL this % files % traverse(initAllFiles)
  CONTAINS
    SUBROUTINE initAllFiles(node)
      TYPE(t_q_linkedlist_node), POINTER, INTENT(inout)  :: node

      SELECT TYPE(v=>node%value)
      CLASS IS (t_qs_output_file)
        CALL v % init()
      CLASS DEFAULT
        ! unexpected type found, raise error
        CALL finish("init_soutputs", "Unexpected class found")
      END SELECT
    END SUBROUTINE initAllFiles
  END SUBROUTINE init_soutputs


  ! ======================================================================================================= !
  !>
  !> Update data from variables (if required, e.g: mean)
  !>
  SUBROUTINE update_soutputs(this)
    CLASS(t_qs_soutputs), INTENT(inout) :: this
    CALL this % files % traverse(runUpdates)
  CONTAINS
    SUBROUTINE runUpdates(node)
      TYPE(t_q_linkedlist_node), POINTER, INTENT(inout)  :: node

      SELECT TYPE(v=>node%value)
      CLASS IS (t_qs_output_file)
        CALL v % update()
      CLASS DEFAULT
        ! unexpected type found, raise error
        CALL finish("update_soutputs", "Unexpected class found")
      END SELECT
    END SUBROUTINE runUpdates
  END SUBROUTINE updatE_soutputs


  ! ======================================================================================================= !
  !>
  !> Write data into the files
  !>
  SUBROUTINE write_soutputs(this, time_d)
    CLASS(t_qs_soutputs), INTENT(inout) :: this
    REAL(wp), INTENT(in) :: time_d !< Time in days
    CALL this % files % traverse(writeAllFiles)
  CONTAINS
    SUBROUTINE writeAllFiles(node)
      TYPE(t_q_linkedlist_node), POINTER, INTENT(inout)  :: node

      SELECT TYPE(v=>node%value)
      CLASS IS (t_qs_output_file)
        CALL v % write(time_d)
      CLASS DEFAULT
        ! unexpected type found, raise error
        CALL finish("write_soutputs", "Unexpected class found")
      END SELECT
    END SUBROUTINE writeAllFiles
  END SUBROUTINE write_soutputs


  ! ======================================================================================================= !
  !>
  !> Reset internal counters at the end of TS
  !>
  SUBROUTINE reset_soutputs(this)
    CLASS(t_qs_soutputs), INTENT(inout) :: this
    CALL this % files % traverse(resetAllFiles)
  CONTAINS
    SUBROUTINE resetAllFiles(node)
      TYPE(t_q_linkedlist_node), POINTER, INTENT(inout)  :: node

      SELECT TYPE(v=>node%value)
      CLASS IS (t_qs_output_file)
        CALL v % reset()
      CLASS DEFAULT
        ! unexpected type found, raise error
        CALL finish("resetAllFiles", "Unexpected class found")
      END SELECT
    END SUBROUTINE resetAllFiles
  END SUBROUTINE reset_soutputs


  ! ======================================================================================================= !
  !>
  !> Reset internal counters at the end of TS
  !>
  SUBROUTINE close_soutputs(this)
    CLASS(t_qs_soutputs), INTENT(inout) :: this
    CALL this % files % traverse(closeAllFiles)
  CONTAINS
    SUBROUTINE closeAllFiles(node)
      TYPE(t_q_linkedlist_node), POINTER, INTENT(inout)  :: node

      SELECT TYPE(v=>node%value)
      CLASS IS (t_qs_output_file)
        CALL v % close()
      CLASS DEFAULT
        ! unexpected type found, raise error
        CALL finish("closeAllFiles", "Unexpected class found")
      END SELECT
    END SUBROUTINE closeAllFiles
  END SUBROUTINE close_soutputs


  ! ======================================================================================================= !
  !>
  !> Create new output file
  !>
  FUNCTION add_file(this, filename) RESULT (fid)
    CLASS(t_qs_soutputs), INTENT(inout) :: this
    CHARACTER(*), INTENT(in) :: filename
    INTEGER :: fid ! Actual position in the list.
    TYPE(t_qs_output_file), POINTER :: of
    ALLOCATE(of)
    of = t_qs_output_file(filename)
    CALL add_to_list(this % files, of)
    fid = this % files % getLength()
  END FUNCTION add_file


  ! ======================================================================================================= !
  !>
  !> Add variable into the file 'fid'
  !>
  SUBROUTINE addvar_soutputs(this, fid, name, units, pdata2d, pdata3d, preop, op)
    CLASS(t_qs_soutputs), INTENT(inout) :: this
    INTEGER, INTENT(in) :: fid
    CHARACTER(*), INTENT(in) :: name
    CHARACTER(*), INTENT(in) :: units
    INTEGER, INTENT(in) :: op
    INTEGER, INTENT(in), OPTIONAL :: preop
    REAL(wp), INTENT(in), POINTER, OPTIONAL :: pdata2d
    REAL(wp), INTENT(in), POINTER, OPTIONAL :: pdata3d(:)

    type(t_q_linkedlist_node), pointer :: inode
    TYPE(t_qs_output_file), pointer :: sfile

    ! find reference
    inode => this % files % getAtIndex(fid)
    SELECT TYPE(v => inode%value)
      TYPE IS (t_qs_output_file)
        sfile => v
    CLASS DEFAULT
      CALL finish("addvar_soutputs", "Unexpected class found")
    END SELECT

    CALL sfile % add_var(name, units, pdata2d, pdata3d, preop, op)
  END SUBROUTINE addvar_soutputs


END MODULE mo_q_dfire_output_txt_util

