!!
!! @brief Contains the routines for plant memory operations 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
MODULE mo_veg_memory

  IMPLICIT NONE

  PUBLIC :: VegPoolAdd,VegPoolSub,VegPoolMul,VegPoolDiv

CONTAINS

 FUNCTION VegPoolAdd(p1,p2) RESULT(p3)
     USE mo_base_memory
     IMPLICIT NONE
     TYPE(VegPool), INTENT(in) :: p1,p2
     TYPE(VegPool) :: p3
     p3%leaf = ElementAdd(p1%leaf,p2%leaf)
     p3%fine_root = ElementAdd(p1%fine_root,p2%fine_root)
     p3%coarse_root = ElementAdd(p1%coarse_root,p2%coarse_root)
     p3%sap_wood = ElementAdd(p1%sap_wood,p2%sap_wood)
     p3%heart_wood = ElementAdd(p1%heart_wood,p2%heart_wood)
     p3%labile = ElementAdd(p1%labile,p2%labile)
     p3%reserve = ElementAdd(p1%reserve,p2%reserve)
     p3%fruit = ElementAdd(p1%fruit,p2%fruit)
     p3%seed_bed = ElementAdd(p1%seed_bed,p2%seed_bed)
 END FUNCTION VegPoolAdd  

 FUNCTION VegPoolSub(p1,p2) RESULT(p3)
     USE mo_base_memory
     IMPLICIT NONE
     TYPE(VegPool), INTENT(in) :: p1,p2
     TYPE(VegPool) :: p3
     p3%leaf = ElementSub(p1%leaf,p2%leaf)
     p3%fine_root = ElementSub(p1%fine_root,p2%fine_root)
     p3%coarse_root = ElementSub(p1%coarse_root,p2%coarse_root)
     p3%sap_wood = ElementSub(p1%sap_wood,p2%sap_wood)
     p3%heart_wood = ElementSub(p1%heart_wood,p2%heart_wood)
     p3%labile = ElementSub(p1%labile,p2%labile)
     p3%reserve = ElementSub(p1%reserve,p2%reserve)
     p3%fruit = ElementSub(p1%fruit,p2%fruit)
     p3%seed_bed = ElementSub(p1%seed_bed,p2%seed_bed)
 END FUNCTION VegPoolSub  

 FUNCTION VegPoolMul(p1,p2) RESULT(p3)
     USE mo_base_memory
     IMPLICIT NONE
     TYPE(VegPool), INTENT(in) :: p1,p2
     TYPE(VegPool) :: p3
     p3%leaf = ElementMul(p1%leaf,p2%leaf)
     p3%fine_root = ElementMul(p1%fine_root,p2%fine_root)
     p3%coarse_root = ElementMul(p1%coarse_root,p2%coarse_root)
     p3%sap_wood = ElementMul(p1%sap_wood,p2%sap_wood)
     p3%heart_wood = ElementMul(p1%heart_wood,p2%heart_wood)
     p3%labile = ElementMul(p1%labile,p2%labile)
     p3%reserve = ElementMul(p1%reserve,p2%reserve)
     p3%fruit = ElementMul(p1%fruit,p2%fruit)
     p3%seed_bed = ElementMul(p1%seed_bed,p2%seed_bed)
 END FUNCTION VegPoolMul  

 FUNCTION VegPoolDiv(p1,p2) RESULT(p3)
     USE mo_base_memory
     IMPLICIT NONE
     TYPE(VegPool), INTENT(in) :: p1,p2
     TYPE(VegPool) :: p3
     p3%leaf = ElementDiv(p1%leaf,p2%leaf)
     p3%fine_root = ElementDiv(p1%fine_root,p2%fine_root)
     p3%coarse_root = ElementDiv(p1%coarse_root,p2%coarse_root)
     p3%sap_wood = ElementDiv(p1%sap_wood,p2%sap_wood)
     p3%heart_wood = ElementDiv(p1%heart_wood,p2%heart_wood)
     p3%labile = ElementDiv(p1%labile,p2%labile)
     p3%reserve = ElementDiv(p1%reserve,p2%reserve)
     p3%fruit = ElementDiv(p1%fruit,p2%fruit)
     p3%seed_bed = ElementDiv(p1%seed_bed,p2%seed_bed)
 END FUNCTION VegPoolDiv  

  SUBROUTINE init_plant_variable_accumulation(tile)
    ! set intial values for parameters that accumulate values (average) across a certain period
    ! called in plant_model_driver section 1.0
    ! the interval is defined in mo_veg_parameters
    USE mo_base_memory
    IMPLICIT NONE
    ! 0.1 Input
    ! 
    ! 0.2 Output
    TYPE(TileClass), DIMENSION(ntile), INTENT(out) :: tile 
    ! 0.3 Local
    ! 

    ! 
    ! 1.0 initialize default values
    ! 
    tile(:)%plant%accumulated%par_sunlit        = 1._wp
    tile(:)%plant%accumulated%par_shaded        = 1._wp
    tile(:)%plant%accumulated%fleaf_sunlit      = 1._wp
    tile(:)%plant%accumulated%beta_air          = 1._wp
    tile(:)%plant%accumulated%beta_soil_ps      = 1._wp
    tile(:)%plant%accumulated%beta_soil_gs      = 1._wp
  END SUBROUTINE init_plant_variable_accumulation

END MODULE mo_veg_memory
