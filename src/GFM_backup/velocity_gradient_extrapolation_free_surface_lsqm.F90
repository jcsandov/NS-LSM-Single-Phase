



subroutine velocity_gradient_extrapolation_free_surface_lsqm( a_coeff_vector    , alpha_vec      , & 
                                                              velocity_gradient , du_dx_fs_lsqm      )

   implicit none

   real (kind = rdf), dimension(1:27), intent(in) :: a_coeff_vector ! in matrix form for an 
                                                                       ! easier handling 

   real (kind = rdf), dimension(1:3), intent(in) :: alpha_vec ! position vector between xs and
                                                              ! xneasrest
   real (kind = rdf), dimension(1:3,1:3), intent(in) :: velocity_gradient ! velocity_gradient at
                                                                          ! inearest, jnearest, knearest

   real (kind = rdf), dimension(1:3,1:3), intent(inout) :: du_dx_fs_lsqm ! velocity gradient at the
                                                                         ! free-surface

   ! local variables
   integer :: igfs, jgfs
   real (kind = rdf), dimension(1:3,1:3) :: fmatrix ! f_ij matrix for linear function approximation 

   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !                            _   _
   !                           |     |  
   !                           | b11 |   (1)
   !                           | c11 |   (2)
   !                           | d11 |   (3)
   !                           | --- |  
   !                           | b21 |   (4)
   !                           | c21 |   (5)
   !                           | d21 |   (6)
   !                           | --- |  
   !                           | b31 |   (7)
   !                           | c31 |   (8)
   !                           | d31 |   (9)
   !                           |     |  
   !                           | --- |  
   !                           |     |  
   !                           | b12 |   (10)
   !                           | c12 |   (11)
   !                           | d12 |   (12)
   !                           | --- |  
   !                           | b22 |   (13)
   !   a_coeff_vector    =     | c22 |   (14)
   !                           | d22 |   (15)
   !                           | --- |  
   !                           | b32 |   (16)
   !                           | c32 |   (17)
   !                           | d32 |   (18)
   !                           |     |  
   !                           | --- |  
   !                           |     |  
   !                           | b13 |   (19)
   !                           | c13 |   (20)
   !                           | d13 |   (21)
   !                           | --- |  
   !                           | b23 |   (22)
   !                           | c23 |   (23)
   !                           | d23 |   (24)
   !                           | --- |  
   !                           | b33 |   (25)
   !                           | c33 |   (26)
   !                           | d33 |   (27)
   !                            -   -  
   !
   !
   ! aij is the coefficient for extrapolating ∂u_i/∂x_j
   !
   !
   ! du_dx_fs_lsqm : cartesian velocity gradient at the free-surface computed
   ! using the coefficients obtained by lsqm.
   ! 
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   ! f_ij = b_ij * α + c_ij * β + d_ij * γ

   fmatrix(1,1) =    a_coeff_vector(1)  * alpha_vec(1) &
                   + a_coeff_vector(2)  * alpha_vec(2) &
                   + a_coeff_vector(3)  * alpha_vec(3)

   fmatrix(2,1) =    a_coeff_vector(4)  * alpha_vec(1) &
                   + a_coeff_vector(5)  * alpha_vec(2) &
                   + a_coeff_vector(6)  * alpha_vec(3)

   fmatrix(3,1) =    a_coeff_vector(7)  * alpha_vec(1) &
                   + a_coeff_vector(8)  * alpha_vec(2) &
                   + a_coeff_vector(9)  * alpha_vec(3)

   fmatrix(1,2) =    a_coeff_vector(10) * alpha_vec(1) &
                   + a_coeff_vector(11) * alpha_vec(2) &
                   + a_coeff_vector(12) * alpha_vec(3)

   fmatrix(2,2) =    a_coeff_vector(13) * alpha_vec(1) &
                   + a_coeff_vector(14) * alpha_vec(2) &
                   + a_coeff_vector(15) * alpha_vec(3)

   fmatrix(3,2) =    a_coeff_vector(16) * alpha_vec(1) &
                   + a_coeff_vector(17) * alpha_vec(2) &
                   + a_coeff_vector(18) * alpha_vec(3)

   fmatrix(1,3) =    a_coeff_vector(19) * alpha_vec(1) &
                   + a_coeff_vector(20) * alpha_vec(2) &
                   + a_coeff_vector(21) * alpha_vec(3)

   fmatrix(2,3) =    a_coeff_vector(22) * alpha_vec(1) &
                   + a_coeff_vector(23) * alpha_vec(2) &
                   + a_coeff_vector(24) * alpha_vec(3)

   fmatrix(3,3) =    a_coeff_vector(25) * alpha_vec(1) &
                   + a_coeff_vector(26) * alpha_vec(2) &
                   + a_coeff_vector(27) * alpha_vec(3)


   ! xs : free-surface location
   !
   ! ∂u_i/∂x^j (xs) = ∂u_i/∂x_j (xs + α) - f_ij
   ! Where ∂u_i/∂x_j (xs + α) is the velocity gradient at the nearest water-phase node
   ! with a valid computed gradient (using the rules defined in the velocity_curv_gradient_tensor
   ! subroutine)

   do jgfs = 1,3 ! j-index for the velocity gradient at the free-surface loop
   do igfs = 1,3 ! i-index for for the velocity gradient the free-surface loop

      du_dx_fs_lsqm(igfs,jgfs) = velocity_gradient(igfs,jgfs) - fmatrix(igfs,jgfs)

   end do
   end do

end subroutine velocity_gradient_extrapolation_free_surface_lsqm