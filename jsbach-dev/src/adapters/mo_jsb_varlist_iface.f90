!>
!! @brief Contains the interface to varlists/streams 
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2012-09-27)
!!
MODULE mo_jsb_varlist_iface

#ifdef __ICON__
  USE mo_var_list,    ONLY: new_var_list, &
                            default_var_list_settings, &
                            get_var_list, &
                            add_var
  USE mo_linked_list, ONLY: t_var_list
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_var_list, new_var_list, default_var_list_settings, get_var_list
  PUBLIC :: add_var
  
  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_varlist_iface'

END MODULE mo_jsb_varlist_iface
