PROGRAM TransitionMatrix
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: n3 = 3, n11 = 11
    REAL(8) :: transition_3x3(n3, n3)
    REAL(8) :: transition_11x11(n11, n11)
    REAL(8) :: vegetation_map(n11)
    INTEGER :: category_mapping(n11)
    INTEGER :: i, j
    REAL(8) :: transition_prob
    REAL(8) :: row_sum
  
    ! Step 1: Define the 3x3 transition matrix (example values)
    DATA transition_3x3 / &
       0.6D0, 0.2D0, 0.2D0, &
       0.1D0, 0.7D0, 0.2D0, &
       0.05D0, 0.15D0, 0.8D0 /
    
    ! Step 2: Define the potential vegetation map fractions for the 11 land cover types
    DATA vegetation_map / &
       0.05D0, 0.10D0, 0.10D0, 0.10D0, &
       0.05D0, 0.05D0, 0.10D0, 0.10D0, &
       0.15D0, 0.10D0, 0.10D0 /
    
    ! Step 3: Define the mapping of 11 land cover types to the 3 broad categories
    ! Mapping: 0 = Crops, 1 = Pasture, 2 = Natural
    DATA category_mapping / &
       3, 3, 3, 3, 3, 3, 3, 3, &  ! Types 1-8 (Natural)
       2, 2, 1 /                  ! Types 9-11 (Pasture and Crops)
    
    ! Step 4: Initialize the 11x11 transition matrix to zero
    transition_11x11 = 0.0D0
  
    ! Step 5: Populate the 11x11 transition matrix
    DO i = 1, n11
       DO j = 1, n11
          ! Get the broad categories for each land cover type
          transition_prob = transition_3x3(category_mapping(i), category_mapping(j))
          
          ! Distribute the transition using the vegetation map fractions
          transition_11x11(i, j) = transition_prob * vegetation_map(i) * vegetation_map(j)
       END DO
    END DO
  
    ! Step 6: Normalize the rows of the 11x11 matrix so that each row sums to 1
    DO i = 1, n11
       row_sum = 0.0D0
       DO j = 1, n11
          row_sum = row_sum + transition_11x11(i, j)
       END DO
       IF (row_sum /= 0.0D0) THEN
          DO j = 1, n11
             transition_11x11(i, j) = transition_11x11(i, j) / row_sum
          END DO
       END IF
    END DO
  
    ! Step 7: Print the resulting 11x11 transition matrix
    PRINT *, '11x11 Transition Matrix:'
    DO i = 1, n11
       WRITE(*,'(11F8.5)') (transition_11x11(i, j), j = 1, n11)
    END DO
  
  END PROGRAM TransitionMatrix
  

ic = 1
DO ic = ics, ice

  DO WHILE (ASSOCIATED(current_tile))
   i_tile = i_tile + 1

   ! Define a 1D vector `TransMtrxVec` for the 3x3 matrix, dimensioned as 9 elements
   ! Read and set transition matrix
   DO idx = 1, 9
      ! Calculate the original row (`i`) and column (`j`) positions from the 1D index `idx
      i = (idx - 1) / 3 + 1
      j = MOD(idx - 1, 3) + 1

      ! Set diagonal elements to 0.0; otherwise, read from input
      IF (i .EQ. j) THEN
         cf_transition(ic, idx, :) = 0.0
      ELSE
         cf_transition(ic, idx, :) = input_file_trans%Read_2d(variable_name = TRIM(transitions(i*j)))
      END IF
   END DO

   ! Loop over all pairs of PFTs and extend the transition matrix for PFTs
   DO idx = 1, npft * npft
      ! Calculate the row (`i`) and column (`j`) indices in the `npft x npft` matrix
      i = (idx - 1) / npft + 1
      j = MOD(idx - 1, npft) + 1

      ! Use the mappings to access the corresponding element in `TransMtrxVec`
      row = LC_type2LC_class_mapping(i)
      col = LC_type2LC_class_mapping(j)
      trans_idx = (row - 1) * 3 + col  ! Mapping to 3x3 transition matrix vector index

      transition_prob = cf_transition(:, trans_idx, :)

      ! Calculate transition probabilities for PFTs
      cf_transition_pfts(ic, idx, :) = cf_current_year(ic, i, :) * transition_prob * cf_current_year(ic, j, :)
   END DO

   ! Initialize transition gain and loss
   cf_transition_gain = 0.0
   cf_transition_loss = 0.0

   ! Calculate transition gain and loss of each PFT in one loop
   DO idx = 1, npft * npft
      ! Determine the (i, j) position for gain and loss calculations
      i = (idx - 1) / npft + 1
      j = MOD(idx - 1, npft) + 1

      ! Calculate transition gain for PFT i -> j, and loss for PFT j -> i
      cf_transition_gain(ic, i_tile, :) = cf_transition_gain(ic, i_tile, :) + &
                                         cf_transition_pfts(ic, idx, :) * cf_current_year(ic, j, :)

      cf_transition_loss(ic, i_tile, :) = cf_transition_loss(ic, i_tile, :) + &
                                         cf_transition_pfts(ic, idx, :) * cf_current_year(ic, i, :)
   END DO

   ! Update cf_current_year by subtracting transition losses and adding transition gains
   cf_current_year(ic, i_tile, :) = cf_current_year(ic, i_tile, :) + cf_transition_gain(ic, i_tile, :) - cf_transition_loss(ic, i_tile, :)

   ! Get the next sibling tile
   current_tile => current_tile%Get_next_sibling_tile()
  END DO
END DO



! If diagonal elements are set to represent the fraction remaining unchanged, apply them directly
! uncomment the part: ! Adjust the diagonal entries so each row sums to 1
! Initialize an array to store the updated cf_current_year
! TODO:
  DO i = 1, npft
   ! For each i-th PFT, multiply its row in cf_current_year by the corresponding diagonal element
   diag_idx = (i - 1) * npft + i  ! Index for the diagonal element in TransMtrx_pfts vector
   cf_current_year(:, ic, i) = cf_current_year(:, ic, i) *cf_transition_pfts(:,npft*npft,:)
  END DO




