
subroutine pressure_gradient_extrapolation_free_surface_lsqm( a_coeff_vector_p, alpha_vec, & 
                                                              pressure_gradient, dp_dx_fs_lsqm)

   implicit none

   real (kind = rdf), dimension(1:12), intent(in) :: a_coeff_vector_p ! in matrix form for an 
                                                                         ! easier handling 

   real (kind = rdf), dimension(1:3), intent(in) :: alpha_vec ! position vector between xs and
                                                              ! xneasrest
   real (kind = rdf), dimension(1:3), intent(in) :: pressure_gradient ! velocity_gradient at
                                                                          ! inearest, jnearest, knearest

   real (kind = rdf), dimension(1:3), intent(inout) :: dp_dx_fs_lsqm ! velocity gradient at the
                                                                         ! free-surface

   ! local variables

   integer :: igfs, jgfs

   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !                       __        __
   !                      |  a1 a2 a3  |
   !                      |  b1 b2 b3  |
   !  a_coeff_vector_p =  |  c1 c2 c3  |
   !                      |  d1 d2 d3  |
   !                       --        --
   !
   !
   ! aij is the coefficient for extrapolating ∂p_i/∂x_j
   !
   !
   ! dp_dx_fs_lsqm : cartesian pressure gradient at the free-surface computed
   ! using the coefficients obtained by lsqm.
   ! 
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   ! xs : free-surface location
   !
   ! ∂p/∂x_i (xs) = ∂p/∂x_i (xs + α) - f_i
   ! Where ∂u_i/∂x_j (xs + α) is the velocity gradient at the nearest water-phase node
   ! with a valid computed gradient (using the rules defined in the velocity_curv_gradient_tensor
   ! subroutine)


   dp_dx_fs_lsqm(1) = pressure_gradient(1) - ( a_coeff_vector_p(1) + a_coeff_vector_p(2)  * alpha_vec(1) &
                                                                   + a_coeff_vector_p(3)  * alpha_vec(2) &
                                                                   + a_coeff_vector_p(4)  * alpha_vec(3) )


   dp_dx_fs_lsqm(2) = pressure_gradient(2) - ( a_coeff_vector_p(5) + a_coeff_vector_p(6)  * alpha_vec(1) &
                                                                   + a_coeff_vector_p(7)  * alpha_vec(2) &
                                                                   + a_coeff_vector_p(8)  * alpha_vec(3) )

   dp_dx_fs_lsqm(3) = pressure_gradient(3) - ( a_coeff_vector_p(9) + a_coeff_vector_p(10) * alpha_vec(1) &
                                                                   + a_coeff_vector_p(11) * alpha_vec(2) &
                                                                   + a_coeff_vector_p(12) * alpha_vec(3) )


end subroutine pressure_gradient_extrapolation_free_surface_lsqm