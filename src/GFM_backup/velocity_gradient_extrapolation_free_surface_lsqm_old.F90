

subroutine velocity_gradient_extrapolation_free_surface_lsqm_old( a_coeff_vector    , alpha_vec      , & 
                                                                  velocity_gradient , du_dx_fs_lsqm      )

   implicit none

   real (kind = rdf), dimension(0:3,1:9), intent(in) :: a_coeff_vector ! in matrix form for an 
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
   !                     __                                   __
   !                    |  a11 a21 a31 a12 a22 a32 a13 a23 a33  |
   !                    |  b11 b21 b31 b12 b22 b32 b13 b23 b33  |
   !  a_coeff_vector =  |  c11 c21 c31 c12 c22 c32 c13 c23 c33  |
   !                    |  d11 d21 d31 d12 d22 d32 d13 d23 d33  |
   !                     --                                   --
   !
   !
   ! aij is the coefficient for extrapolating ∂u_i/∂x_j
   !
   !
   ! du_dx_fs_lsqm : cartesian velocity gradient at the free-surface computed
   ! using the coefficients obtained by lsqm.
   ! 
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   ! f_ij = a_ij + b_ij * α + c_ij * β + d_ij * γ

   fmatrix(1,1) = a_coeff_vector(0,1) + a_coeff_vector(1,1)*alpha_vec(1) &
                                    & + a_coeff_vector(2,1)*alpha_vec(2) &
                                    & + a_coeff_vector(3,1)*alpha_vec(3)

   fmatrix(2,1) = a_coeff_vector(0,2) + a_coeff_vector(1,2)*alpha_vec(1) &
                                    & + a_coeff_vector(2,2)*alpha_vec(2) &
                                    & + a_coeff_vector(3,2)*alpha_vec(3)

   fmatrix(3,1) = a_coeff_vector(0,3) + a_coeff_vector(1,3)*alpha_vec(1) &
                                    & + a_coeff_vector(2,3)*alpha_vec(2) &
                                    & + a_coeff_vector(3,3)*alpha_vec(3)

   fmatrix(1,2) = a_coeff_vector(0,4) + a_coeff_vector(1,4)*alpha_vec(1) &
                                    & + a_coeff_vector(2,4)*alpha_vec(2) &
                                    & + a_coeff_vector(3,4)*alpha_vec(3)

   fmatrix(2,2) = a_coeff_vector(0,5) + a_coeff_vector(1,5)*alpha_vec(1) &
                                    & + a_coeff_vector(2,5)*alpha_vec(2) &
                                    & + a_coeff_vector(3,5)*alpha_vec(3)

   fmatrix(3,2) = a_coeff_vector(0,6) + a_coeff_vector(1,6)*alpha_vec(1) &
                                    & + a_coeff_vector(2,6)*alpha_vec(2) &
                                    & + a_coeff_vector(3,6)*alpha_vec(3)

   fmatrix(1,3) = a_coeff_vector(0,7) + a_coeff_vector(1,7)*alpha_vec(1) &
                                    & + a_coeff_vector(2,7)*alpha_vec(2) &
                                    & + a_coeff_vector(3,7)*alpha_vec(3)

   fmatrix(2,3) = a_coeff_vector(0,8) + a_coeff_vector(1,8)*alpha_vec(1) &
                                    & + a_coeff_vector(2,8)*alpha_vec(2) &
                                    & + a_coeff_vector(3,8)*alpha_vec(3)

   fmatrix(3,3) = a_coeff_vector(0,9) + a_coeff_vector(1,9)*alpha_vec(1) &
                                    & + a_coeff_vector(2,9)*alpha_vec(2) &
                                    & + a_coeff_vector(3,9)*alpha_vec(3)


   ! xs : free-surface location
   !
   ! ∂u_i/∂x^j (xs) = ∂u_i/∂x_j (xs + α) - f_ij
   ! Where ∂u_i/∂x_j (xs + α) is the velocity gradient at the nearest water-phase node
   ! with a valid computed gradient (using the rules defined in the velocity_curv_gradient_tensor
   ! subroutine)

   do jgfs = 1,3 ! j-index for the velocity gradient at the free-surface loop
      do igfs= 1,3 ! i-index for for the velocity gradient the free-surface loop

            du_dx_fs_lsqm(igfs,jgfs) = velocity_gradient(igfs,jgfs) - fmatrix(igfs,jgfs)

      end do
   end do

end subroutine velocity_gradient_extrapolation_free_surface_lsqm_old
