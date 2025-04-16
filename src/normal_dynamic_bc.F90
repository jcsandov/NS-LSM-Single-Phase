subroutine normal_dynamic_bc( du_dx_fs_lsqm , xnut_fs , nvec, pfs )

   !--------------------------------------------------------------------------
   ! Description: 
   ! ------------
   ! Locally determines the pressure at the free surface using the dynamic 
   ! boundary condition:
   !
   !
   !                /  1           \   /  ∂ui     ∂uj  \   
   ! -pfs ni ni  + ( ----- +  νsgs  ) (  ----- + -----  ) ni nj = 0     
   !                \  Re          /   \  ∂xj     ∂xi  /   
   !
   !--------------------------------------------------------------------------

   implicit none

   ! Input
   real(kind = rdf), dimension(1:3,1:3), intent(in) :: du_dx_fs_lsqm ! vel grad at the fs calculated
                                                                     ! with the TDBC
   real(kind = rdf), dimension(1:3)    , intent(in) :: nvec ! normal vector at the free surface
   real(kind = rdf)                    , intent(in) :: xnut_fs ! xnut interpolated at the free-surface


   ! Output
   real(kind = rdf), intent(out) :: pfs ! pressure at the free surface obtained from the NDBC


   ! Local variables
   real(kind = rdf) :: aux , v1

   integer :: ii,jj

   aux = zero


   if ( zero_pressure_fs ) then
   
      pfs = zero

   else

      do jj = 1,3
      do ii = 1,3
      
         aux = aux + ( du_dx_fs_lsqm(ii,jj) + du_dx_fs_lsqm(jj,ii) ) * nvec(ii) * nvec(jj) 
      
      end do
      end do
      
      aux = aux * (one / ren + xnut_fs )

      pfs = aux / ( nvec(1)*nvec(1) + nvec(2)*nvec(2) + nvec(3)*nvec(3) )

   end if

end subroutine normal_dynamic_bc

