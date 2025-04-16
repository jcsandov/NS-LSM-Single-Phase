subroutine velocity_extrapolation_free_surface_lsqm ( alpha_vec, i, j, k, du_dx_fs_lsqm, u_fs_lsqm)

   implicit none

   real (kind = rdf), dimension(1:3), intent(in) :: alpha_vec ! position vector between xs and
                                                              ! x(i,j,k) (xs and x as vectors)
   integer, intent(in) :: i,j,k
   real (kind = rdf), dimension(1:3,1:3), intent(in) :: du_dx_fs_lsqm ! velocity gradient at the
                                                                         ! free-surface
   
   ! I declared this variable as inout to see it and initialise it in the main code before calling                                                                         
   real (kind = rdf), dimension(1:3), intent(inout) :: u_fs_lsqm ! velocity at the free-surface

   ! local variables

   integer :: iufs, tsum

   do iufs = 1,3 ! index for velocity free-surface computation (ufs). iufs = 1 -> u, iufs = 2 -> v... 
      
      ! u_fs_lsqm(iufs) = q(iufs+1) because q(2,i,j,k) = u(i,j,k), q(3,i,j,k) = v(i,j,k) ...         

      u_fs_lsqm(iufs) = q(iufs+1, i, j, k) 

      do tsum = 1,3 ! summation index from Taylor expansion: u(xs) = u(xs+α) - ∂u/∂x*α - ∂u/∂y*β - ∂u/∂z*γ
         
         u_fs_lsqm(iufs) = u_fs_lsqm(iufs) - du_dx_fs_lsqm(iufs,tsum)*alpha_vec(tsum)

      end do
   end do

end subroutine velocity_extrapolation_free_surface_lsqm