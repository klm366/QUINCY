!>
!! @brief Contains the memory class for the soil processes.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_soil_memory_class

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_soil_memory, t_soil_phys_state, &
            t_soil_energy_state, t_soil_water_state, &
            t_soil_energy_flux, t_soil_water_flux, &
            copy_soil_mem, add_soil_mem

  TYPE t_soil_phys_state
    REAL(wp), POINTER :: &
      soil_depth   (:,  :),    & !< Soil depth derived from textures (bedrock) (for hydrology) [m]
      soil_depth_l (:,:,:),    & !< Depth of each soil layer that can be saturated with water (until bedrock) [m]
      vpor         (:,  :),    & !< Volumetric soil porosity []
      psdi         (:,  :),    & !< Soil pore size distribution index []
      bclapp       (:,  :),    & !< Exponent B in Clapp and Hornberger
      mpot         (:,  :),    & !< Matrix potential [m]
      shc          (:,  :),    & !< Saturated hydraulic conductivity [m/s]
      vhcap        (:,  :),    & !< Volumetric heat capacity of dry soil [J m-3 K-1]
      hdiff        (:,  :),    & !< Thermal diffusivity of soil [m2 s-1]
      alb_vis_soil (:,  :),    & !< Albedo of the soil in the visible range
      alb_nir_soil (:,  :)       !< Albedo of the soil in the NIR range
  END type t_soil_phys_state

  TYPE t_soil_energy_state
    REAL(wp), POINTER :: &
      t_soil        (:,:,:), &
      t_soil_acoef  (:,:,:), &
      t_soil_bcoef  (:,:,:), &
      hcap_grnd     (:,  :)     !< Heat capacity of the ground [J m-2 K-1]
  END type t_soil_energy_state

  TYPE t_soil_water_state
    REAL(wp), POINTER :: &
      wfract_soil   (:,  :),   & !< Wet (skin reservoir) fraction of soil []
      wsr_soil      (:,  :),   & !< Water content in skin reservoir on soil [m]
      sfract_soil   (:,  :),   & !< Snow fraction of soil []
      wsn_soil      (:,  :),   & !< Water content of snow reservoir of soil [m water equivalent]
      ws_l          (:,:,:),   & !< Water content in soil layers [m]
      ws_sat_l      (:,:,:),   & !< Water content of soil layers at saturation (derived from soil porosity) [m]
      ws_fc_l       (:,:,:),   & !< Water content of soil layers at field capacity (derived from volumetric field capacity)[m]
      ws_pwp_l      (:,:,:),   & !< Water content of soil layers at permanent wilting point (derived from vol. perm wilt point) [m]
      ws            (:  ,:)      !< Water content in the whole soil colums [m]
  END type t_soil_water_state

  TYPE t_soil_energy_flux
    REAL(wp), POINTER :: &
      grnd_hflx        (:,:)      !< ground heat flux
  END type t_soil_energy_flux

  TYPE t_soil_water_flux
    REAL(wp), POINTER :: &
      & runoff         (:,:), & !< Accumulated surface runoff
      & drainage       (:,:)     !< Accumulated drainage
  END type t_soil_water_flux

  TYPE t_soil_memory
     TYPE(t_soil_phys_state),   POINTER :: phys_state   => NULL()
     TYPE(t_soil_energy_state), POINTER :: energy_state => NULL()
     TYPE(t_soil_water_state),  POINTER :: water_state  => NULL()
     TYPE(t_soil_energy_flux),  POINTER :: energy_flux  => NULL()
     TYPE(t_soil_water_flux),   POINTER :: water_flux   => NULL()
  END type t_soil_memory

  ! Note: Would be nicer to be able to overload the "=" operator, but that would be more
  ! complicated because we only want to copy parts of the arrays (ics:ice,iblk)
  INTERFACE copy_soil_mem
    MODULE PROCEDURE copy_t_soil_energy_state
  END INTERFACE copy_soil_mem

  INTERFACE add_soil_mem
    MODULE PROCEDURE add_t_soil_energy_state
  END INTERFACE add_soil_mem

  CHARACTER(len=*), PARAMETER :: modname = 'mo_soil_memory_class'

CONTAINS

  SUBROUTINE copy_t_soil_energy_state(in, out, iblk, ics, ice)

    TYPE(t_soil_energy_state), INTENT(in)  :: in
    TYPE(t_soil_energy_state), INTENT(out) :: out
    INTEGER,                   INTENT(in)  :: iblk, ics, ice

    out%t_soil        (ics:ice, :, iblk) = in%t_soil        (ics:ice, :, iblk)
    out%t_soil_acoef  (ics:ice, :, iblk) = in%t_soil_acoef  (ics:ice, :, iblk)
    out%t_soil_bcoef  (ics:ice, :, iblk) = in%t_soil_bcoef  (ics:ice, :, iblk)

  END SUBROUTINE copy_t_soil_energy_state

  SUBROUTINE add_t_soil_energy_state(in1, in2, out, iblk, ics, ice)

    TYPE(t_soil_energy_state), INTENT(in)  :: in1, in2
    TYPE(t_soil_energy_state), INTENT(out) :: out
    INTEGER,                   INTENT(in)  :: iblk, ics, ice

    out%t_soil        (ics:ice, :, iblk) = in1%t_soil        (ics:ice, :, iblk) + in2%t_soil        (ics:ice, :, iblk)
    out%t_soil_acoef  (ics:ice, :, iblk) = in1%t_soil_acoef  (ics:ice, :, iblk) + in2%t_soil_acoef  (ics:ice, :, iblk)
    out%t_soil_bcoef  (ics:ice, :, iblk) = in1%t_soil_bcoef  (ics:ice, :, iblk) + in2%t_soil_bcoef  (ics:ice, :, iblk)

  END SUBROUTINE add_t_soil_energy_state

END MODULE mo_soil_memory_class
