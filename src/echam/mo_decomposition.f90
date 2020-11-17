!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!
MODULE mo_decomposition
  !
  ! This module defines the data types holding meta information
  !   for data decomposition in parallel mode
  !
  ! Authors:
  !
  ! A. Rhodin, MPI, August 1999, original source
  ! A. Rhodin, DWD/MPI, August 2001, ffsl decomposition 
  ! T. Diehl,  DKRZ, September 2001, spitfire decomposition
  ! L. Kornblueh, merging of ffsl, spitfire 
  ! A. Rhodin, DWD/MPI, May 2002, blocking (nproma)

  IMPLICIT NONE

  PRIVATE
  !
  ! Data type declaration
  !
  PUBLIC :: pe_decomposed        ! data type with decomposition of a single PE 
  !
  ! Information regarding the decomposition for 
  ! Flux-Form Semi-Lagrangian transport scheme

  TYPE pe_ffsl
    !
    SEQUENCE
    ! 
    INTEGER           :: nlat  ! number of latitudes on this pe                
    INTEGER           :: lats  ! southern  latitude  on this pe (S=1,N=nlatg)  
    INTEGER           :: latn  ! northern  latitudes on this pe                
    INTEGER           :: pe_s  ! processor id of neighbour in the south        
    INTEGER           :: pe_n  ! processor id of neighbour in the north        
    INTEGER           :: pe_x  ! processor to exchange hemispheres with        
    INTEGER           :: nlatx ! number of latitudes to exchange with pe_x     
  END TYPE pe_ffsl

  ! Data type for local information of decomposition per PE 

  TYPE pe_decomposed
     !
     SEQUENCE
     !
     ! Information regarding the whole model domain and discretization

     INTEGER          :: nlon     ! number of longitudes
     INTEGER          :: nlat     ! number of latitudes
     INTEGER          :: npts = -1! number of global compressed grid points (e.g. land)
     LOGICAL, POINTER :: mask(:,:) => NULL() ! global mask
     INTEGER          :: nlev     ! number of levels
     INTEGER          :: nm       ! max. wavenumber used(triangular truncation)
     INTEGER ,POINTER :: nnp(:)   ! number of points on each m-column
     INTEGER ,POINTER :: nmp(:)   ! displacement of the first point of
                                  ! m-columns with respect to the first point 
                                  ! of the first m-column

     ! General information on the decomposition

     INTEGER          :: d_nprocs ! # of PEs for debugging/non-debugging domain
     INTEGER          :: spe      ! index # of first PE of domain
     INTEGER          :: epe      ! index # of last PE of domain
     INTEGER          :: nprocb   !#of PEs for dimension that counts longitudes
     INTEGER          :: nproca   !#of PEs for dimension that counts latitudes

     INTEGER ,POINTER :: mapmesh(:,:) ! indirection array mapping from a
                                      ! logical 2-d mesh to the processor index
                                      ! numbers

     ! local information

     INTEGER          :: pe        ! PE id 
     INTEGER          :: set_b     ! PE id in direction of longitudes
     INTEGER          :: set_a     ! PE id in direction of latitudes 
     LOGICAL          :: col_1d    ! 1d column model(s) on this pe

     ! Grid space decomposition 

     INTEGER          :: nglat     ! number of latitudes  on PE
     INTEGER          :: nglatmax  ! max. number of latitudes per PE
     INTEGER          :: nglon     ! number of longitudes on PE
     INTEGER          :: nglonmax  ! max. number of longitudes per PE
     INTEGER          :: nglh(2)   ! number of latitudes on each hemisphere
     INTEGER          :: glats(2)  ! start values of latitudes
     INTEGER          :: glate(2)  ! end values of latitudes
     INTEGER          :: glons(2)  ! start values of longitudes
     INTEGER          :: glone(2)  ! end values of longitudes
     INTEGER ,POINTER :: glat(:)   ! global latitude index N->S
     INTEGER ,POINTER :: glon(:)   ! offset to global longitude

     INTEGER          :: ngpts = -1! number of compressed grid points on PE
     INTEGER          :: gptss = -1! start value of grid indices
     INTEGER          :: gptse = -1! end value of grid indices

     ! Irregular grid

     INTEGER          :: ngpblks   ! number of rows
     INTEGER          :: nproma    ! number of columns
     INTEGER          :: npromz    ! number of columns in last row
     LOGICAL          :: lreg      ! flag for regular grid 

     ! Fourier space decomposition

     LOGICAL          :: lfused    ! true if this PE used in Fourier space
     INTEGER          :: nflat     ! number of latitudes on PE
     INTEGER          :: nflev     ! number of levels on PE
     INTEGER          :: nflevp1   ! number of levels+1 on PE
     INTEGER          :: flats(2)  ! start values of latitudes (on row)
     INTEGER          :: flate(2)  ! end values of latitudes (on row)
     INTEGER          :: flevs     ! start values of levels (on column)
     INTEGER          :: fleve     ! end values of levels (on column)

     ! Legendre space decomposition

     ! Row of PEs with same set_a
     INTEGER          :: nlm       ! number of local wave numbers m handled
     INTEGER ,POINTER :: lm(:)     ! actual local wave numbers m handled
     INTEGER          :: lnsp      ! number of complex spectral coefficients  
     INTEGER ,POINTER :: nlmp(:)   ! displacement of the first point of columns
     INTEGER ,POINTER :: nlnp(:)   ! number of points on each column
     INTEGER          :: nlnm0     ! number of coeff. with m=0 on this pe
     INTEGER ,POINTER :: intr(:)   ! index array used by transpose routine

     ! Column of PEs with same set_b
     INTEGER          :: nllev     ! number of levels
     INTEGER          :: nllevp1   ! number of levels+1
     INTEGER          :: llevs     ! start values of levels
     INTEGER          :: lleve     ! end values of levels

     ! Spectral space decomposition

     ! local PE
     INTEGER          :: snsp      ! number of spectral coefficients
     INTEGER          :: snsp2     ! 2*number of spectral coefficients
     INTEGER          :: ssps      ! first spectral coefficient
     INTEGER          :: sspe      ! last  spectral coefficient

     LOGICAL          :: lfirstc   ! true, if first global coeff (m=0,n=0) on
     ! this PE (for nudging)
     INTEGER          :: ifirstc   ! location of first global coeff on this PE
     INTEGER ,POINTER :: np1(:)    ! value of (n+1) for all coeffs of this PE 
     INTEGER ,POINTER :: np1i(:)   ! index of value of (n+1) for all coeffs
                                   ! of this PE 
     INTEGER ,POINTER :: mymsp(:)  ! value of m for all coeffs of this PE
     INTEGER          :: nns       ! number of different n-values for this PE
     INTEGER ,POINTER :: nindex(:) ! the nns elements contain the
     ! values of (n+1)

     INTEGER          :: nsm       ! number of wavenumbers per PE
     INTEGER ,POINTER :: sm (:)    ! actual local wave numbers handled
     INTEGER ,POINTER :: snnp(:)   ! number of n coeff. per wave number m
     INTEGER ,POINTER :: snn0(:)   ! first coeff. n for a given m
     INTEGER          :: nsnm0     ! number of coeffs with m=0 on this pe

     ! Flux-Form Semi-Lagrangian transport scheme decomposition

     TYPE (pe_ffsl)  :: ffsl
    
    ! Mark last element of derived type to allow proper usage later 
    ! in derived MPI datatypes 
    INTEGER            :: last          ! just dummy

  END TYPE pe_decomposed

END MODULE mo_decomposition
