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

  PUBLIC :: initialize_trainable_parameters, update_trainable_parameters

  REAL(wp) ::                                                   &
    & layer1_weight_global(39,  128), layer1_bias_global(128),  &
    & layer2_weight_global(128, 128), layer2_bias_global(128),  &
    & layer3_weight_global(128, 128), layer3_bias_global(128),  &
    & layer4_weight_global(128, 128), layer4_bias_global(128),  &
    & layer5_weight_global(128, 128), layer5_bias_global(128),  &
    & layer6_weight_global(128,   1), layer6_bias_global(1)

  
  CHARACTER(len=*), parameter :: path = '/work/mj0143/b301058/fire_online_0_base_update/central/'

  !COMMON /trainable_parameters/ layer1_weight_global, layer1_bias_global, &
  !  layer2_weight_global, layer2_bias_global, layer3_weight_global, layer3_bias_global, &
  !  layer4_weight_global, layer4_bias_global, layer5_weight_global, layer5_bias_global, &
  !  layer6_weight_global, layer6_bias_global

  CHARACTER(len=*), PARAMETER :: modname = 'mo_global_variable'

CONTAINS

  SUBROUTINE initialize_trainable_parameters( &
    & layer1_weight, layer1_bias,             & ! in
    & layer2_weight, layer2_bias,             & ! in
    & layer3_weight, layer3_bias,             & ! in
    & layer4_weight, layer4_bias,             & ! in
    & layer5_weight, layer5_bias,             & ! in
    & layer6_weight, layer6_bias              ) ! in

    ! Input Arguments
    REAL(wp), INTENT(in) ::                       &
      & layer1_weight(39, 128), layer1_bias(128), &
      & layer2_weight(128,128), layer2_bias(128), &
      & layer3_weight(128,128), layer3_bias(128), &
      & layer4_weight(128,128), layer4_bias(128), &
      & layer5_weight(128,128), layer5_bias(128), &
      & layer6_weight(128,  1), layer6_bias(1)

    INTEGER :: pid, getpid
    CHARACTER(len=200) :: w1, w2, w3 ,w4, w5, w6
    CHARACTER(len=200) :: b1, b2, b3 ,b4, b5, b6

    pid = getpid()
    pid = mod(pid, 8)

    11 format('W1_',I1,'.bin') 
    12 format('W2_',I1,'.bin')
    13 format('W3_',I1,'.bin')
    14 format('W4_',I1,'.bin')
    15 format('W5_',I1,'.bin')
    16 format('W6_',I1,'.bin')
    21 format('b1_',I1,'.bin') 
    22 format('b2_',I1,'.bin')
    23 format('b3_',I1,'.bin')
    24 format('b4_',I1,'.bin')
    25 format('b5_',I1,'.bin')
    26 format('b6_',I1,'.bin')

    write(w1, 11) pid
    write(w2, 12) pid
    write(w3, 13) pid
    write(w4, 14) pid
    write(w5, 15) pid
    write(w6, 16) pid
    write(b1, 21) pid
    write(b2, 22) pid
    write(b3, 23) pid
    write(b4, 24) pid
    write(b5, 25) pid
    write(b6, 26) pid

    layer1_weight_global(:,:) = layer1_weight(:,:)
    layer2_weight_global(:,:) = layer2_weight(:,:)
    layer3_weight_global(:,:) = layer3_weight(:,:)
    layer4_weight_global(:,:) = layer4_weight(:,:)
    layer5_weight_global(:,:) = layer5_weight(:,:)
    layer6_weight_global(:,:) = layer6_weight(:,:)

    layer1_bias_global(:) = layer1_bias(:)
    layer2_bias_global(:) = layer2_bias(:)
    layer3_bias_global(:) = layer3_bias(:)
    layer4_bias_global(:) = layer4_bias(:)
    layer5_bias_global(:) = layer5_bias(:)
    layer6_bias_global(:) = layer6_bias(:)

    open(701,file=path//w1, access='direct', form='unformatted', recl=39*128*4)
    write(701,rec=1) layer1_weight_global
    open(702,file=path//w2, access='direct', form='unformatted', recl=128*128*4)
    write(702,rec=1) layer2_weight_global
    open(703,file=path//w3, access='direct', form='unformatted', recl=128*128*4)
    write(703,rec=1) layer3_weight_global
    open(704,file=path//w4, access='direct', form='unformatted', recl=128*128*4)
    write(704,rec=1) layer4_weight_global
    open(705,file=path//w5, access='direct', form='unformatted', recl=128*128*4)
    write(705,rec=1) layer5_weight_global
    open(706,file=path//w6, access='direct', form='unformatted', recl=128*1*4)
    write(706,rec=1) layer6_weight_global

    open(801,file=path//b1, access='direct', form='unformatted', recl=128*4)
    write(801,rec=1) layer1_bias_global
    open(802,file=path//b2, access='direct', form='unformatted', recl=128*4)
    write(802,rec=1) layer2_bias_global
    open(803,file=path//b3, access='direct', form='unformatted', recl=128*4)
    write(803,rec=1) layer3_bias_global
    open(804,file=path//b4, access='direct', form='unformatted', recl=128*4)
    write(804,rec=1) layer4_bias_global
    open(805,file=path//b5, access='direct', form='unformatted', recl=128*4)
    write(805,rec=1) layer5_bias_global
    open(806,file=path//b6, access='direct', form='unformatted', recl=1*4)
    write(806,rec=1) layer6_bias_global

  END SUBROUTINE initialize_trainable_parameters


  SUBROUTINE prior_trainable_parameters(      &
    & layer1_weight, layer1_bias,             & ! inout
    & layer2_weight, layer2_bias,             & ! inout
    & layer3_weight, layer3_bias,             & ! inout
    & layer4_weight, layer4_bias,             & ! inout
    & layer5_weight, layer5_bias,             & ! inout
    & layer6_weight, layer6_bias              ) ! inout

    ! Input Arguments
    REAL(wp), INTENT(inout) ::                    &
      & layer1_weight(39, 128), layer1_bias(128), &
      & layer2_weight(128,128), layer2_bias(128), &
      & layer3_weight(128,128), layer3_bias(128), &
      & layer4_weight(128,128), layer4_bias(128), &
      & layer5_weight(128,128), layer5_bias(128), &
      & layer6_weight(128,  1), layer6_bias(1)

    INTEGER :: pid, getpid
    CHARACTER(len=200) :: w1, w2, w3 ,w4, w5, w6
    CHARACTER(len=200) :: b1, b2, b3 ,b4, b5, b6

    pid = getpid()
    pid = mod(pid, 8)

    11 format('W1_',I1,'.bin') 
    12 format('W2_',I1,'.bin')
    13 format('W3_',I1,'.bin')
    14 format('W4_',I1,'.bin')
    15 format('W5_',I1,'.bin')
    16 format('W6_',I1,'.bin')
    21 format('b1_',I1,'.bin') 
    22 format('b2_',I1,'.bin')
    23 format('b3_',I1,'.bin')
    24 format('b4_',I1,'.bin')
    25 format('b5_',I1,'.bin')
    26 format('b6_',I1,'.bin')

    write(w1, 11) pid
    write(w2, 12) pid
    write(w3, 13) pid
    write(w4, 14) pid
    write(w5, 15) pid
    write(w6, 16) pid
    write(b1, 21) pid
    write(b2, 22) pid
    write(b3, 23) pid
    write(b4, 24) pid
    write(b5, 25) pid
    write(b6, 26) pid

    open(101,file=path//w1, access='direct', form='unformatted', recl=39*128*4)
    read(101,rec=1) layer1_weight_global
    open(102,file=path//w2, access='direct', form='unformatted', recl=128*128*4)
    read(102,rec=1) layer2_weight_global
    open(103,file=path//w3, access='direct', form='unformatted', recl=128*128*4)
    read(103,rec=1) layer3_weight_global
    open(104,file=path//w4, access='direct', form='unformatted', recl=128*128*4)
    read(104,rec=1) layer4_weight_global
    open(105,file=path//w5, access='direct', form='unformatted', recl=128*128*4)
    read(105,rec=1) layer5_weight_global
    open(106,file=path//w6, access='direct', form='unformatted', recl=128*1*4)
    read(106,rec=1) layer6_weight_global

    open(201,file=path//b1, access='direct', form='unformatted', recl=128*4)
    read(201,rec=1) layer1_bias_global
    open(202,file=path//b2, access='direct', form='unformatted', recl=128*4)
    read(202,rec=1) layer2_bias_global
    open(203,file=path//b3, access='direct', form='unformatted', recl=128*4)
    read(203,rec=1) layer3_bias_global
    open(204,file=path//b4, access='direct', form='unformatted', recl=128*4)
    read(204,rec=1) layer4_bias_global
    open(205,file=path//b5, access='direct', form='unformatted', recl=128*4)
    read(205,rec=1) layer5_bias_global
    open(206,file=path//b6, access='direct', form='unformatted', recl=1*4)
    read(206,rec=1) layer6_bias_global


    layer1_weight(:,:) = layer1_weight_global(:,:)
    layer2_weight(:,:) = layer2_weight_global(:,:)
    layer3_weight(:,:) = layer3_weight_global(:,:)
    layer4_weight(:,:) = layer4_weight_global(:,:)
    layer5_weight(:,:) = layer5_weight_global(:,:)
    layer6_weight(:,:) = layer6_weight_global(:,:)

    layer1_bias(:) = layer1_bias_global(:)
    layer2_bias(:) = layer2_bias_global(:)
    layer3_bias(:) = layer3_bias_global(:)
    layer4_bias(:) = layer4_bias_global(:)
    layer5_bias(:) = layer5_bias_global(:)
    layer6_bias(:) = layer6_bias_global(:)

  END SUBROUTINE prior_trainable_parameters


  SUBROUTINE update_trainable_parameters(     &
    & layer1_weight, layer1_bias,             & ! in
    & layer2_weight, layer2_bias,             & ! in
    & layer3_weight, layer3_bias,             & ! in
    & layer4_weight, layer4_bias,             & ! in
    & layer5_weight, layer5_bias,             & ! in
    & layer6_weight, layer6_bias              ) ! in

    ! Input Arguments
    REAL(wp), INTENT(in) ::                       &
      & layer1_weight(39, 128), layer1_bias(128), &
      & layer2_weight(128,128), layer2_bias(128), &
      & layer3_weight(128,128), layer3_bias(128), &
      & layer4_weight(128,128), layer4_bias(128), &
      & layer5_weight(128,128), layer5_bias(128), &
      & layer6_weight(128,  1), layer6_bias(1)

INTEGER :: pid, getpid
    CHARACTER(len=200) :: w1, w2, w3 ,w4, w5, w6
    CHARACTER(len=200) :: b1, b2, b3 ,b4, b5, b6

    pid = getpid()
    pid = mod(pid, 8)

    11 format('W1_',I1,'.bin') 
    12 format('W2_',I1,'.bin')
    13 format('W3_',I1,'.bin')
    14 format('W4_',I1,'.bin')
    15 format('W5_',I1,'.bin')
    16 format('W6_',I1,'.bin')
    21 format('b1_',I1,'.bin') 
    22 format('b2_',I1,'.bin')
    23 format('b3_',I1,'.bin')
    24 format('b4_',I1,'.bin')
    25 format('b5_',I1,'.bin')
    26 format('b6_',I1,'.bin')

    write(w1, 11) pid
    write(w2, 12) pid
    write(w3, 13) pid
    write(w4, 14) pid
    write(w5, 15) pid
    write(w6, 16) pid
    write(b1, 21) pid
    write(b2, 22) pid
    write(b3, 23) pid
    write(b4, 24) pid
    write(b5, 25) pid
    write(b6, 26) pid
    

    open(301,file=path//w1, access='direct', form='unformatted', recl=39*128*4)
    layer1_weight_global(:,:) = layer1_weight(:,:)
    write(301,rec=1) layer1_weight_global

    open(302,file=path//w2, access='direct', form='unformatted', recl=128*128*4)
    layer2_weight_global(:,:) = layer2_weight(:,:)
    write(302,rec=1) layer2_weight_global

    open(303,file=path//w2, access='direct', form='unformatted', recl=128*128*4)
    layer3_weight_global(:,:) = layer3_weight(:,:)
    write(303,rec=1) layer3_weight_global

    open(304,file=path//w2, access='direct', form='unformatted', recl=128*128*4)
    layer4_weight_global(:,:) = layer4_weight(:,:)
    write(304,rec=1) layer4_weight_global

    open(305,file=path//w2, access='direct', form='unformatted', recl=128*128*4)
    layer5_weight_global(:,:) = layer5_weight(:,:)
    write(305,rec=1) layer5_weight_global

    open(306,file=path//w2, access='direct', form='unformatted', recl=128*1*4)
    layer6_weight_global(:,:) = layer6_weight(:,:)
    write(306,rec=1) layer6_weight_global

    open(401,file=path//b1, access='direct', form='unformatted', recl=128*4)
    layer1_bias_global(:) = layer1_bias(:)
    write(401,rec=1) layer1_bias_global

    open(402,file=path//b2, access='direct', form='unformatted', recl=128*4)
    layer2_bias_global(:) = layer2_bias(:)
    write(402,rec=1) layer2_bias_global

    open(403,file=path//b3, access='direct', form='unformatted', recl=128*4)
    layer3_bias_global(:) = layer3_bias(:)
    write(403,rec=1) layer3_bias_global

    open(404,file=path//b4, access='direct', form='unformatted', recl=128*4)
    layer4_bias_global(:) = layer4_bias(:)
    write(404,rec=1) layer4_bias_global

    open(405,file=path//b5, access='direct', form='unformatted', recl=128*4)
    layer5_bias_global(:) = layer5_bias(:)
    write(405,rec=1) layer5_bias_global

    open(406,file=path//b6, access='direct', form='unformatted', recl=1*4)
    layer6_bias_global(:) = layer6_bias(:)
    write(406,rec=1) layer6_bias_global

  END SUBROUTINE update_trainable_parameters


  SUBROUTINE merge_trainable_parameters()

    REAL(wp) ::                                         &
      & layer1_weight_(39, 128,8), layer1_bias_(128,8), &
      & layer2_weight_(128,128,8), layer2_bias_(128,8), &
      & layer3_weight_(128,128,8), layer3_bias_(128,8), &
      & layer4_weight_(128,128,8), layer4_bias_(128,8), &
      & layer5_weight_(128,128,8), layer5_bias_(128,8), &
      & layer6_weight_(128,  1,8), layer6_bias_(1,8),   &

      & layer1_weight(39, 128),    layer1_bias(128),    &
      & layer2_weight(128,128),    layer2_bias(128),    &
      & layer3_weight(128,128),    layer3_bias(128),    &
      & layer4_weight(128,128),    layer4_bias(128),    &
      & layer5_weight(128,128),    layer5_bias(128),    &
      & layer6_weight(128,  1),    layer6_bias(1)

    INTEGER :: i, i_, read_idx

    open(501,file='/work/mj0143/b301058/fire_online_0_base_update/central/W1_0.bin', &
      access='direct', form='unformatted', recl=39*128*4)
    read(501,rec=1) layer1_weight_(:,:,1)
    open(502,file='/work/mj0143/b301058/fire_online_0_base_update/central/W1_1.bin', &
      access='direct', form='unformatted', recl=39*128*4)
    read(502,rec=1) layer1_weight_(:,:,2)
    open(503,file='/work/mj0143/b301058/fire_online_0_base_update/central/W1_2.bin', &
      access='direct', form='unformatted', recl=39*128*4)
    read(503,rec=1) layer1_weight_(:,:,3)
    open(504,file='/work/mj0143/b301058/fire_online_0_base_update/central/W1_3.bin', &
      access='direct', form='unformatted', recl=39*128*4)
    read(504,rec=1) layer1_weight_(:,:,4)
    open(505,file='/work/mj0143/b301058/fire_online_0_base_update/central/W1_4.bin', &
      access='direct', form='unformatted', recl=39*128*4)
    read(505,rec=1) layer1_weight_(:,:,5)
    open(506,file='/work/mj0143/b301058/fire_online_0_base_update/central/W1_5.bin', &
      access='direct', form='unformatted', recl=39*128*4)
    read(506,rec=1) layer1_weight_(:,:,6)
    open(507,file='/work/mj0143/b301058/fire_online_0_base_update/central/W1_6.bin', &
      access='direct', form='unformatted', recl=39*128*4)
    read(507,rec=1) layer1_weight_(:,:,7)
    open(508,file='/work/mj0143/b301058/fire_online_0_base_update/central/W1_7.bin', &
      access='direct', form='unformatted', recl=39*128*4)
    read(508,rec=1) layer1_weight_(:,:,8)
    layer1_weight(:,:) = SUM(layer1_weight_, dim=3)/8._wp

    write(501,rec=1) layer1_weight
    write(502,rec=1) layer1_weight
    write(503,rec=1) layer1_weight
    write(504,rec=1) layer1_weight
    write(505,rec=1) layer1_weight
    write(506,rec=1) layer1_weight
    write(507,rec=1) layer1_weight
    write(508,rec=1) layer1_weight


    open(511,file='/work/mj0143/b301058/fire_online_0_base_update/central/W2_0.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(511,rec=1) layer2_weight_(:,:,1)
    open(512,file='/work/mj0143/b301058/fire_online_0_base_update/central/W2_1.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(512,rec=1) layer2_weight_(:,:,2)
    open(513,file='/work/mj0143/b301058/fire_online_0_base_update/central/W2_2.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(513,rec=1) layer2_weight_(:,:,3)
    open(514,file='/work/mj0143/b301058/fire_online_0_base_update/central/W2_3.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(514,rec=1) layer2_weight_(:,:,4)
    open(515,file='/work/mj0143/b301058/fire_online_0_base_update/central/W2_4.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(515,rec=1) layer2_weight_(:,:,5)
    open(516,file='/work/mj0143/b301058/fire_online_0_base_update/central/W2_5.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(516,rec=1) layer2_weight_(:,:,6)
    open(517,file='/work/mj0143/b301058/fire_online_0_base_update/central/W2_6.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(517,rec=1) layer2_weight_(:,:,7)
    open(518,file='/work/mj0143/b301058/fire_online_0_base_update/central/W2_7.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(518,rec=1) layer2_weight_(:,:,8)
    layer2_weight(:,:) = SUM(layer2_weight_, dim=3)/8._wp

    write(511,rec=1) layer2_weight
    write(512,rec=1) layer2_weight
    write(513,rec=1) layer2_weight
    write(514,rec=1) layer2_weight
    write(515,rec=1) layer2_weight
    write(516,rec=1) layer2_weight
    write(517,rec=1) layer2_weight
    write(518,rec=1) layer2_weight


    open(521,file='/work/mj0143/b301058/fire_online_0_base_update/central/W3_0.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(521,rec=1) layer3_weight_(:,:,1)
    open(522,file='/work/mj0143/b301058/fire_online_0_base_update/central/W3_1.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(522,rec=1) layer3_weight_(:,:,2)
    open(523,file='/work/mj0143/b301058/fire_online_0_base_update/central/W3_2.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(523,rec=1) layer3_weight_(:,:,3)
    open(524,file='/work/mj0143/b301058/fire_online_0_base_update/central/W3_3.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(524,rec=1) layer3_weight_(:,:,4)
    open(525,file='/work/mj0143/b301058/fire_online_0_base_update/central/W3_4.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(525,rec=1) layer3_weight_(:,:,5)
    open(526,file='/work/mj0143/b301058/fire_online_0_base_update/central/W3_5.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(526,rec=1) layer3_weight_(:,:,6)
    open(527,file='/work/mj0143/b301058/fire_online_0_base_update/central/W3_6.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(527,rec=1) layer3_weight_(:,:,7)
    open(528,file='/work/mj0143/b301058/fire_online_0_base_update/central/W3_7.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(528,rec=1) layer3_weight_(:,:,8)
    layer3_weight(:,:) = SUM(layer3_weight_, dim=3)/8._wp

    write(521,rec=1) layer3_weight
    write(522,rec=1) layer3_weight
    write(523,rec=1) layer3_weight
    write(524,rec=1) layer3_weight
    write(525,rec=1) layer3_weight
    write(526,rec=1) layer3_weight
    write(527,rec=1) layer3_weight
    write(528,rec=1) layer3_weight


    open(531,file='/work/mj0143/b301058/fire_online_0_base_update/central/W4_0.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(531,rec=1) layer4_weight_(:,:,1)
    open(532,file='/work/mj0143/b301058/fire_online_0_base_update/central/W4_1.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(532,rec=1) layer4_weight_(:,:,2)
    open(533,file='/work/mj0143/b301058/fire_online_0_base_update/central/W4_2.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(533,rec=1) layer4_weight_(:,:,3)
    open(534,file='/work/mj0143/b301058/fire_online_0_base_update/central/W4_3.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(534,rec=1) layer4_weight_(:,:,4)
    open(535,file='/work/mj0143/b301058/fire_online_0_base_update/central/W4_4.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(535,rec=1) layer4_weight_(:,:,5)
    open(536,file='/work/mj0143/b301058/fire_online_0_base_update/central/W4_5.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(536,rec=1) layer4_weight_(:,:,6)
    open(537,file='/work/mj0143/b301058/fire_online_0_base_update/central/W4_6.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(537,rec=1) layer4_weight_(:,:,7)
    open(538,file='/work/mj0143/b301058/fire_online_0_base_update/central/W4_7.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(538,rec=1) layer4_weight_(:,:,8)
    layer4_weight(:,:) = SUM(layer4_weight_, dim=3)/8._wp

    write(531,rec=1) layer4_weight
    write(532,rec=1) layer4_weight
    write(533,rec=1) layer4_weight
    write(534,rec=1) layer4_weight
    write(535,rec=1) layer4_weight
    write(536,rec=1) layer4_weight
    write(537,rec=1) layer4_weight
    write(538,rec=1) layer4_weight


    open(541,file='/work/mj0143/b301058/fire_online_0_base_update/central/W5_0.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(541,rec=1) layer5_weight_(:,:,1)
    open(542,file='/work/mj0143/b301058/fire_online_0_base_update/central/W5_1.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(542,rec=1) layer5_weight_(:,:,2)
    open(543,file='/work/mj0143/b301058/fire_online_0_base_update/central/W5_2.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(543,rec=1) layer5_weight_(:,:,3)
    open(544,file='/work/mj0143/b301058/fire_online_0_base_update/central/W5_3.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(544,rec=1) layer5_weight_(:,:,4)
    open(545,file='/work/mj0143/b301058/fire_online_0_base_update/central/W5_4.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(545,rec=1) layer5_weight_(:,:,5)
    open(546,file='/work/mj0143/b301058/fire_online_0_base_update/central/W5_5.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(546,rec=1) layer5_weight_(:,:,6)
    open(547,file='/work/mj0143/b301058/fire_online_0_base_update/central/W5_6.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(547,rec=1) layer5_weight_(:,:,7)
    open(548,file='/work/mj0143/b301058/fire_online_0_base_update/central/W5_7.bin', &
      access='direct', form='unformatted', recl=128*128*4)
    read(548,rec=1) layer5_weight_(:,:,8)
    layer5_weight(:,:) = SUM(layer5_weight_, dim=3)/8._wp

    write(541,rec=1) layer5_weight
    write(542,rec=1) layer5_weight
    write(543,rec=1) layer5_weight
    write(544,rec=1) layer5_weight
    write(545,rec=1) layer5_weight
    write(546,rec=1) layer5_weight
    write(547,rec=1) layer5_weight
    write(548,rec=1) layer5_weight


    open(551,file='/work/mj0143/b301058/fire_online_0_base_update/central/W6_0.bin', &
      access='direct', form='unformatted', recl=128*1*4)
    read(551,rec=1) layer6_weight_(:,:,1)
    open(552,file='/work/mj0143/b301058/fire_online_0_base_update/central/W6_1.bin', &
      access='direct', form='unformatted', recl=128*1*4)
    read(552,rec=1) layer6_weight_(:,:,2)
    open(553,file='/work/mj0143/b301058/fire_online_0_base_update/central/W6_2.bin', &
      access='direct', form='unformatted', recl=128*1*4)
    read(553,rec=1) layer6_weight_(:,:,3)
    open(554,file='/work/mj0143/b301058/fire_online_0_base_update/central/W6_3.bin', &
      access='direct', form='unformatted', recl=128*1*4)
    read(554,rec=1) layer6_weight_(:,:,4)
    open(555,file='/work/mj0143/b301058/fire_online_0_base_update/central/W6_4.bin', &
      access='direct', form='unformatted', recl=128*1*4)
    read(555,rec=1) layer6_weight_(:,:,5)
    open(556,file='/work/mj0143/b301058/fire_online_0_base_update/central/W6_5.bin', &
      access='direct', form='unformatted', recl=128*1*4)
    read(556,rec=1) layer6_weight_(:,:,6)
    open(557,file='/work/mj0143/b301058/fire_online_0_base_update/central/W6_6.bin', &
      access='direct', form='unformatted', recl=128*1*4)
    read(557,rec=1) layer6_weight_(:,:,7)
    open(558,file='/work/mj0143/b301058/fire_online_0_base_update/central/W6_7.bin', &
      access='direct', form='unformatted', recl=128*1*4)
    read(558,rec=1) layer6_weight_(:,:,8)
    layer6_weight(:,:) = SUM(layer6_weight_, dim=3)/8._wp

    write(551,rec=1) layer6_weight
    write(552,rec=1) layer6_weight
    write(553,rec=1) layer6_weight
    write(554,rec=1) layer6_weight
    write(555,rec=1) layer6_weight
    write(556,rec=1) layer6_weight
    write(557,rec=1) layer6_weight
    write(558,rec=1) layer6_weight


    open(601,file='/work/mj0143/b301058/fire_online_0_base_update/central/b1_0.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(601,rec=1) layer1_bias_(:,1)
    open(602,file='/work/mj0143/b301058/fire_online_0_base_update/central/b1_1.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(602,rec=1) layer1_bias_(:,2)
    open(603,file='/work/mj0143/b301058/fire_online_0_base_update/central/b1_2.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(603,rec=1) layer1_bias_(:,3)
    open(604,file='/work/mj0143/b301058/fire_online_0_base_update/central/b1_3.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(604,rec=1) layer1_bias_(:,4)
    open(605,file='/work/mj0143/b301058/fire_online_0_base_update/central/b1_4.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(605,rec=1) layer1_bias_(:,5)
    open(606,file='/work/mj0143/b301058/fire_online_0_base_update/central/b1_5.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(606,rec=1) layer1_bias_(:,6)
    open(607,file='/work/mj0143/b301058/fire_online_0_base_update/central/b1_6.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(607,rec=1) layer1_bias_(:,7)
    open(608,file='/work/mj0143/b301058/fire_online_0_base_update/central/b1_7.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(608,rec=1) layer1_bias_(:,8)
    layer1_bias(:) = SUM(layer1_bias_, dim=2)/8._wp

    write(601,rec=1) layer1_bias
    write(602,rec=1) layer1_bias
    write(603,rec=1) layer1_bias
    write(604,rec=1) layer1_bias
    write(605,rec=1) layer1_bias
    write(606,rec=1) layer1_bias
    write(607,rec=1) layer1_bias
    write(608,rec=1) layer1_bias


    open(611,file='/work/mj0143/b301058/fire_online_0_base_update/central/b2_0.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(611,rec=1) layer2_bias_(:,1)
    open(612,file='/work/mj0143/b301058/fire_online_0_base_update/central/b2_1.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(612,rec=1) layer2_bias_(:,2)
    open(613,file='/work/mj0143/b301058/fire_online_0_base_update/central/b2_2.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(613,rec=1) layer2_bias_(:,3)
    open(614,file='/work/mj0143/b301058/fire_online_0_base_update/central/b2_3.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(614,rec=1) layer2_bias_(:,4)
    open(615,file='/work/mj0143/b301058/fire_online_0_base_update/central/b2_4.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(615,rec=1) layer2_bias_(:,5)
    open(616,file='/work/mj0143/b301058/fire_online_0_base_update/central/b2_5.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(616,rec=1) layer2_bias_(:,6)
    open(617,file='/work/mj0143/b301058/fire_online_0_base_update/central/b2_6.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(617,rec=1) layer2_bias_(:,7)
    open(618,file='/work/mj0143/b301058/fire_online_0_base_update/central/b2_7.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(618,rec=1) layer2_bias_(:,8)
    layer2_bias(:) = SUM(layer2_bias_, dim=2)/8._wp

    write(611,rec=1) layer2_bias
    write(612,rec=1) layer2_bias
    write(613,rec=1) layer2_bias
    write(614,rec=1) layer2_bias
    write(615,rec=1) layer2_bias
    write(616,rec=1) layer2_bias
    write(617,rec=1) layer2_bias
    write(618,rec=1) layer2_bias


    open(621,file='/work/mj0143/b301058/fire_online_0_base_update/central/b3_0.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(621,rec=1) layer3_bias_(:,1)
    open(622,file='/work/mj0143/b301058/fire_online_0_base_update/central/b3_1.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(622,rec=1) layer3_bias_(:,2)
    open(623,file='/work/mj0143/b301058/fire_online_0_base_update/central/b3_2.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(623,rec=1) layer3_bias_(:,3)
    open(624,file='/work/mj0143/b301058/fire_online_0_base_update/central/b3_3.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(624,rec=1) layer3_bias_(:,4)
    open(625,file='/work/mj0143/b301058/fire_online_0_base_update/central/b3_4.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(625,rec=1) layer3_bias_(:,5)
    open(626,file='/work/mj0143/b301058/fire_online_0_base_update/central/b3_5.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(626,rec=1) layer3_bias_(:,6)
    open(627,file='/work/mj0143/b301058/fire_online_0_base_update/central/b3_6.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(627,rec=1) layer3_bias_(:,7)
    open(628,file='/work/mj0143/b301058/fire_online_0_base_update/central/b3_7.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(628,rec=1) layer3_bias_(:,8)
    layer3_bias(:) = SUM(layer3_bias_, dim=2)/8._wp

    write(621,rec=1) layer3_bias
    write(622,rec=1) layer3_bias
    write(623,rec=1) layer3_bias
    write(624,rec=1) layer3_bias
    write(625,rec=1) layer3_bias
    write(626,rec=1) layer3_bias
    write(627,rec=1) layer3_bias
    write(628,rec=1) layer3_bias


    open(631,file='/work/mj0143/b301058/fire_online_0_base_update/central/b4_0.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(631,rec=1) layer4_bias_(:,1)
    open(632,file='/work/mj0143/b301058/fire_online_0_base_update/central/b4_1.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(632,rec=1) layer4_bias_(:,2)
    open(633,file='/work/mj0143/b301058/fire_online_0_base_update/central/b4_2.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(633,rec=1) layer4_bias_(:,3)
    open(634,file='/work/mj0143/b301058/fire_online_0_base_update/central/b4_3.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(634,rec=1) layer4_bias_(:,4)
    open(635,file='/work/mj0143/b301058/fire_online_0_base_update/central/b4_4.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(635,rec=1) layer4_bias_(:,5)
    open(636,file='/work/mj0143/b301058/fire_online_0_base_update/central/b4_5.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(636,rec=1) layer4_bias_(:,6)
    open(637,file='/work/mj0143/b301058/fire_online_0_base_update/central/b4_6.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(637,rec=1) layer4_bias_(:,7)
    open(638,file='/work/mj0143/b301058/fire_online_0_base_update/central/b4_7.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(638,rec=1) layer4_bias_(:,8)
    layer4_bias(:) = SUM(layer4_bias_, dim=2)/8._wp

    write(631,rec=1) layer4_bias
    write(632,rec=1) layer4_bias
    write(633,rec=1) layer4_bias
    write(634,rec=1) layer4_bias
    write(635,rec=1) layer4_bias
    write(636,rec=1) layer4_bias
    write(637,rec=1) layer4_bias
    write(638,rec=1) layer4_bias


    open(641,file='/work/mj0143/b301058/fire_online_0_base_update/central/b5_0.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(641,rec=1) layer5_bias_(:,1)
    open(642,file='/work/mj0143/b301058/fire_online_0_base_update/central/b5_1.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(642,rec=1) layer5_bias_(:,2)
    open(643,file='/work/mj0143/b301058/fire_online_0_base_update/central/b5_2.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(643,rec=1) layer5_bias_(:,3)
    open(644,file='/work/mj0143/b301058/fire_online_0_base_update/central/b5_3.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(644,rec=1) layer5_bias_(:,4)
    open(645,file='/work/mj0143/b301058/fire_online_0_base_update/central/b5_4.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(645,rec=1) layer5_bias_(:,5)
    open(646,file='/work/mj0143/b301058/fire_online_0_base_update/central/b5_5.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(646,rec=1) layer5_bias_(:,6)
    open(647,file='/work/mj0143/b301058/fire_online_0_base_update/central/b5_6.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(647,rec=1) layer5_bias_(:,7)
    open(648,file='/work/mj0143/b301058/fire_online_0_base_update/central/b5_7.bin', &
      access='direct', form='unformatted', recl=128*4)
    read(648,rec=1) layer5_bias_(:,8)
    layer5_bias(:) = SUM(layer5_bias_, dim=2)/8._wp

    write(641,rec=1) layer5_bias
    write(642,rec=1) layer5_bias
    write(643,rec=1) layer5_bias
    write(644,rec=1) layer5_bias
    write(645,rec=1) layer5_bias
    write(646,rec=1) layer5_bias
    write(647,rec=1) layer5_bias
    write(648,rec=1) layer5_bias


    open(651,file='/work/mj0143/b301058/fire_online_0_base_update/central/b6_0.bin', &
      access='direct', form='unformatted', recl=1*4)
    read(651,rec=1) layer6_bias_(:,1)
    open(652,file='/work/mj0143/b301058/fire_online_0_base_update/central/b6_1.bin', &
      access='direct', form='unformatted', recl=1*4)
    read(652,rec=1) layer6_bias_(:,2)
    open(653,file='/work/mj0143/b301058/fire_online_0_base_update/central/b6_2.bin', &
      access='direct', form='unformatted', recl=1*4)
    read(653,rec=1) layer6_bias_(:,3)
    open(654,file='/work/mj0143/b301058/fire_online_0_base_update/central/b6_3.bin', &
      access='direct', form='unformatted', recl=1*4)
    read(654,rec=1) layer6_bias_(:,4)
    open(655,file='/work/mj0143/b301058/fire_online_0_base_update/central/b6_4.bin', &
      access='direct', form='unformatted', recl=1*4)
    read(655,rec=1) layer6_bias_(:,5)
    open(656,file='/work/mj0143/b301058/fire_online_0_base_update/central/b6_5.bin', &
      access='direct', form='unformatted', recl=1*4)
    read(656,rec=1) layer6_bias_(:,6)
    open(657,file='/work/mj0143/b301058/fire_online_0_base_update/central/b6_6.bin', &
      access='direct', form='unformatted', recl=1*4)
    read(657,rec=1) layer6_bias_(:,7)
    open(658,file='/work/mj0143/b301058/fire_online_0_base_update/central/b6_7.bin', &
      access='direct', form='unformatted', recl=1*4)
    read(658,rec=1) layer6_bias_(:,8)
    layer6_bias(:) = SUM(layer6_bias_, dim=2)/8._wp

    write(651,rec=1) layer6_bias
    write(652,rec=1) layer6_bias
    write(653,rec=1) layer6_bias
    write(654,rec=1) layer6_bias
    write(655,rec=1) layer6_bias
    write(656,rec=1) layer6_bias
    write(657,rec=1) layer6_bias
    write(658,rec=1) layer6_bias

  END SUBROUTINE merge_trainable_parameters

END MODULE mo_global_variable
