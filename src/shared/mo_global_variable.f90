!>
!! @par Copyright
!!  This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!!  Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!
!! @brief diverse subroutines and functions for tile & forcing variables
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2014-03-17)
!!

MODULE mo_global_variable

  USE mo_kind,      ONLY: wp

  IMPLICIT NONE

  PUBLIC :: update_pid_batch, load_global_batch, update_trainable_parameters, load_trained_parameters

  REAL(wp) ::                                                      &
    & layer1_weight_global(53, 64*12), layer1_bias_global(64*12),  &
    & layer2_weight_global(64, 64*12), layer2_bias_global(64*12),  &
    & layer3_weight_global(64, 64*12), layer3_bias_global(64*12),  &
    & layer4_weight_global(64,  1*12), layer4_bias_global( 1*12)

  REAL(wp) :: batch_global(1280, 57)   ! 160*8
  REAL(wp) :: batch_pid   (160,  57)   ! 

  CHARACTER(len=*), parameter :: path = '/work/mj0143/b301058/fire_online_2_koppen_lv2/central/'

  INTEGER :: file_idx = 0
  INTEGER :: first_iblk = 0
  INTEGER :: step = 0

  CHARACTER(len=*), PARAMETER :: modname = 'mo_global_variable'

CONTAINS
  
  SUBROUTINE get_first_iblk(x)
    INTEGER, INTENT(in) :: x

    step = step + 1
    IF (step .eq. 1) THEN
      first_iblk = x
      print*, first_iblk
    ENDIF
  END SUBROUTINE get_first_iblk


  SUBROUTINE update_pid_batch(                          &
    & V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10, &
    & V11, V12, V13, V14, V15, V16, V17, V18, V19, V20, &
    & V21, V22, V23, V24, V25, V26, V27, V28, V29, V30, &
    & V31, V32, V33, V34, V35, V36, V37, V38, V39, V40, &
    & V41, V42, V43, V44, V45, V46, V47, V48, V49, V50, &
    & V51, V52, V53, RHC,   Y,   R,   K,   n, iblk)

    ! Input Arguments
    REAL(wp), INTENT(in) ::                               &
      & V1,  V2,  V3,  V4,  V5,  V6,  V7,  V8,  V9,  V10, &
      & V11, V12, V13, V14, V15, V16, V17, V18, V19, V20, &
      & V21, V22, V23, V24, V25, V26, V27, V28, V29, V30, &
      & V31, V32, V33, V34, V35, V36, V37, V38, V39, V40, &
      & V41, V42, V43, V44, V45, V46, V47, V48, V49, V50, &
      & V51, V52, V53, RHC,   Y,   R,   K

    INTEGER, INTENT(in) ::  n, iblk

    CHARACTER(len=200) :: batch_name
    INTEGER :: pid, getpid, idx, i
    pid = getpid()
    pid = mod(pid, 8)

    idx = n + 8*(iblk-1)

    batch_pid(idx,1)  = V1;  batch_pid(idx,2)  = V2 ; batch_pid(idx,3)  = V3
    batch_pid(idx,4)  = V4;  batch_pid(idx,5)  = V5 ; batch_pid(idx,6)  = V6
    batch_pid(idx,7)  = V7;  batch_pid(idx,8)  = V8 ; batch_pid(idx,9)  = V9
    batch_pid(idx,10) = V10; batch_pid(idx,11) = V11; batch_pid(idx,12) = V12
    batch_pid(idx,13) = V13; batch_pid(idx,14) = V14; batch_pid(idx,15) = V15
    batch_pid(idx,16) = V16; batch_pid(idx,17) = V17; batch_pid(idx,18) = V18
    batch_pid(idx,19) = V19; batch_pid(idx,20) = V20; batch_pid(idx,21) = V21
    batch_pid(idx,22) = V22; batch_pid(idx,23) = V23; batch_pid(idx,24) = V24
    batch_pid(idx,25) = V25; batch_pid(idx,26) = V26; batch_pid(idx,27) = V27
    batch_pid(idx,28) = V28; batch_pid(idx,29) = V29; batch_pid(idx,30) = V30
    batch_pid(idx,31) = V31; batch_pid(idx,32) = V32; batch_pid(idx,33) = V33
    batch_pid(idx,34) = V34; batch_pid(idx,35) = V35; batch_pid(idx,36) = V36
    batch_pid(idx,37) = V37; batch_pid(idx,38) = V38; batch_pid(idx,39) = V39
    batch_pid(idx,40) = V40; batch_pid(idx,41) = V41; batch_pid(idx,42) = V42
    batch_pid(idx,43) = V43; batch_pid(idx,44) = V44; batch_pid(idx,45) = V45
    batch_pid(idx,46) = V46; batch_pid(idx,47) = V47; batch_pid(idx,48) = V48
    batch_pid(idx,49) = V49; batch_pid(idx,50) = V50; batch_pid(idx,51) = V51
    batch_pid(idx,52) = V52; batch_pid(idx,53) = V53

    batch_pid(idx,54) = RHC; batch_pid(idx,55) = Y  ; batch_pid(idx,56) = R
    batch_pid(idx,57) = K

    101 format('batch_p',I1,'.bin') 
    write(batch_name, 101) pid

    IF (iblk .eq. 1) THEN
      open(101,file=path//batch_name, access='direct', form='unformatted', recl=160*57*4)
      !DO i=1,16
      write(101,rec=1) batch_pid
      !ENDDO
    ENDIF

  END SUBROUTINE update_pid_batch


  SUBROUTINE load_global_batch()
    ! Input Arguments
    CHARACTER(len=200) :: batch_0, batch_1, batch_2, batch_3, batch_4, batch_5, batch_6, batch_7
    INTEGER :: pid, getpid, i
    REAL(wp) :: u

    !pid = getpid()
    !pid = mod(pid, 8)+1

    !CALL RANDOM_NUMBER(u)
    !i = INT(FLOOR((16)*u))


    111 format('batch_p0.bin'); write(batch_0, 111) 
    112 format('batch_p1.bin'); write(batch_1, 112)
    113 format('batch_p2.bin'); write(batch_2, 113)
    114 format('batch_p3.bin'); write(batch_3, 114) 
    115 format('batch_p4.bin'); write(batch_4, 115)
    116 format('batch_p5.bin'); write(batch_5, 116)
    117 format('batch_p6.bin'); write(batch_6, 117)
    118 format('batch_p7.bin'); write(batch_7, 118)

    open(111,file=path//batch_0, access='direct', form='unformatted', recl=160*57*4)
    read(111,rec=1) batch_global(1:160,:)
    open(112,file=path//batch_1, access='direct', form='unformatted', recl=160*57*4)
    read(112,rec=1) batch_global(161:320,:)
    open(113,file=path//batch_2, access='direct', form='unformatted', recl=160*57*4)
    read(113,rec=1) batch_global(321:480,:)
    open(114,file=path//batch_3, access='direct', form='unformatted', recl=160*57*4)
    read(114,rec=1) batch_global(481:640,:)
    open(115,file=path//batch_4, access='direct', form='unformatted', recl=160*57*4)
    read(115,rec=1) batch_global(641:800,:)
    open(116,file=path//batch_5, access='direct', form='unformatted', recl=160*57*4)
    read(116,rec=1) batch_global(801:960,:)
    open(117,file=path//batch_6, access='direct', form='unformatted', recl=160*57*4)
    read(117,rec=1) batch_global(961:1120,:)
    open(118,file=path//batch_7, access='direct', form='unformatted', recl=160*57*4)
    read(118,rec=1) batch_global(1121:1280,:)

  END SUBROUTINE load_global_batch


  SUBROUTINE update_trainable_parameters(     &
    & layer1_weight, layer1_bias,             & ! in
    & layer2_weight, layer2_bias,             & ! in
    & layer3_weight, layer3_bias,             & ! in
    & layer4_weight, layer4_bias              ) ! in

    ! Input Arguments
    REAL(wp), INTENT(in) ::                          &
      & layer1_weight(53,64*12), layer1_bias(64*12), &
      & layer2_weight(64,64*12), layer2_bias(64*12), &
      & layer3_weight(64,64*12), layer3_bias(64*12), &
      & layer4_weight(64, 1*12), layer4_bias( 1*12)

    INTEGER :: pid, getpid, i
    INTEGER :: file_idx_, file_idx_w
    CHARACTER(len=200) :: w1, w2, w3, w4
    CHARACTER(len=200) :: b1, b2, b3, b4

    !file_idx_ = MOD(file_idx, 512) + 1
    !file_idx = file_idx + 1
    !file_idx_w = MOD(file_idx, 512) + 1

    !pid = getpid()
    !pid = mod(pid, 8)

    !print*, pid, file_idx, file_idx_, file_idx_w

    201 format('W1.bin'); write(w1, 201)
    202 format('W2.bin'); write(w2, 202)
    203 format('W3.bin'); write(w3, 203) 
    204 format('W4.bin'); write(w4, 204) 

    211 format('b1.bin'); write(b1, 211)
    212 format('b2.bin'); write(b2, 212)
    213 format('b3.bin'); write(b3, 213)
    214 format('b4.bin'); write(b4, 214)

    open(201,file=path//w1, access='direct', form='unformatted', recl=53*64*12*4)
    open(202,file=path//w2, access='direct', form='unformatted', recl=64*64*12*4)
    open(203,file=path//w3, access='direct', form='unformatted', recl=64*64*12*4)
    open(204,file=path//w4, access='direct', form='unformatted', recl=64* 1*12*4)

    open(211,file=path//b1, access='direct', form='unformatted', recl=64*12*4)
    open(212,file=path//b2, access='direct', form='unformatted', recl=64*12*4)
    open(213,file=path//b3, access='direct', form='unformatted', recl=64*12*4)
    open(214,file=path//b4, access='direct', form='unformatted', recl= 1*12*4)

    DO i = 1,8
      write(201,rec=i) layer1_weight
      write(202,rec=i) layer2_weight
      write(203,rec=i) layer3_weight
      write(204,rec=i) layer4_weight

      write(211,rec=i) layer1_bias
      write(212,rec=i) layer2_bias
      write(213,rec=i) layer3_bias
      write(214,rec=i) layer4_bias
    ENDDO

    layer1_weight_global(:,:) = layer1_weight(:,:)
    layer2_weight_global(:,:) = layer2_weight(:,:)
    layer3_weight_global(:,:) = layer3_weight(:,:)
    layer4_weight_global(:,:) = layer4_weight(:,:)

    layer1_bias_global(:)     = layer1_bias(:)
    layer2_bias_global(:)     = layer2_bias(:)
    layer3_bias_global(:)     = layer3_bias(:)
    layer4_bias_global(:)     = layer4_bias(:)

  END SUBROUTINE update_trainable_parameters


  SUBROUTINE load_trained_parameters(         &
    & layer1_weight, layer1_bias,             & ! in
    & layer2_weight, layer2_bias,             & ! in
    & layer3_weight, layer3_bias,             & ! in
    & layer4_weight, layer4_bias              ) ! in

    ! Input Arguments
    REAL(wp), INTENT(inout) ::                       &
      & layer1_weight(53,64*12), layer1_bias(64*12), &
      & layer2_weight(64,64*12), layer2_bias(64*12), &
      & layer3_weight(64,64*12), layer3_bias(64*12), &
      & layer4_weight(64, 1*12), layer4_bias( 1*12)

    INTEGER :: pid, getpid, i
    CHARACTER(len=200) :: w1, w2, w3 ,w4
    CHARACTER(len=200) :: b1, b2, b3 ,b4

    pid = getpid()
    pid = mod(pid, 8)+1

    301 format('W1.bin'); write(w1, 301)
    302 format('W2.bin'); write(w2, 302)
    303 format('W3.bin'); write(w3, 303) 
    304 format('W4.bin'); write(w4, 304) 

    311 format('b1.bin'); write(b1, 311)
    312 format('b2.bin'); write(b2, 312)
    313 format('b3.bin'); write(b3, 313)
    314 format('b4.bin'); write(b4, 314)

    open(301,file=path//w1, access='direct', form='unformatted', recl=53*64*12*4)
    open(302,file=path//w2, access='direct', form='unformatted', recl=64*64*12*4)
    open(303,file=path//w3, access='direct', form='unformatted', recl=64*64*12*4)
    open(304,file=path//w4, access='direct', form='unformatted', recl=64* 1*12*4)

    open(311,file=path//b1, access='direct', form='unformatted', recl=64*12*4)
    open(312,file=path//b2, access='direct', form='unformatted', recl=64*12*4)
    open(313,file=path//b3, access='direct', form='unformatted', recl=64*12*4)
    open(314,file=path//b4, access='direct', form='unformatted', recl= 1*12*4)

    read(301,rec=pid) layer1_weight
    read(302,rec=pid) layer2_weight
    read(303,rec=pid) layer3_weight
    read(304,rec=pid) layer4_weight

    read(311,rec=pid) layer1_bias
    read(312,rec=pid) layer2_bias
    read(313,rec=pid) layer3_bias
    read(314,rec=pid) layer4_bias

    layer1_weight_global(:,:) = layer1_weight(:,:)
    layer2_weight_global(:,:) = layer2_weight(:,:)
    layer3_weight_global(:,:) = layer3_weight(:,:)
    layer4_weight_global(:,:) = layer4_weight(:,:)

    layer1_bias_global(:)     = layer1_bias(:)
    layer2_bias_global(:)     = layer2_bias(:)
    layer3_bias_global(:)     = layer3_bias(:)
    layer4_bias_global(:)     = layer4_bias(:)

  END SUBROUTINE load_trained_parameters


END MODULE mo_global_variable
