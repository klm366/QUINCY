!>
!! @brief Contains adapter to different utilities
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2016-02-29)
!!
MODULE mo_jsb_utils_iface
#ifndef __NO_JSBACH__

  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: finish, message, message_text

#ifdef __ICON__
  USE mo_fortran_tools, ONLY: assign_if_present
#else
  USE mo_memory_base,   ONLY: assign_if_present
#endif

#endif
END MODULE mo_jsb_utils_iface
