subroutine p_no_slip_wall_bcond(b)

   implicit none

   integer, intent(in) :: b

   ! local variables
   real (kind = rdf) :: c1 , c2 , c3      ! coeff of each component of the extrapolation
   real (kind = rdf) :: g0 , g1           ! g22 at the wall and 1 node away from it
   real (kind = rdf) :: ret0 , ret1       ! (1/Re + xnut) at the wall and 1 node away from it
   real (kind = rdf) :: f1aux , f2aux     ! internal derivative of the viscous flux
   real (kind = rdf) :: fv1star , fv2star ! external derivative of the viscous flux

   integer :: i,j,k

   eta_dir: select case(b)
      
      case(3)

         j  = j_mysta
         c1 = four / three
         c2 = -one / three
         c3 = -two / three

         do k = k_mysta , k_myend
         do i = i_mysta , i_myend

            ! g0 = g22 at the wall (0 nodes away from the wall)
            g0 = eta(1,i,j,k)**two + eta(2,i,j,k)**two 

            ! g1 = g22 1 node away from the wall
            g1 = eta(1,i,j+1,k)**two + eta(2,i,j+1,k)**two 

            ! ret0 = ( 1/Re + xnut ) at the wall
            ret0 = one / ren + xnut(i,j,k) 

            ! ret1 = ( 1/Re + xnut ) 1 node away from the wall
            ret1 = one / ren + xnut(i,j+1,k) 

            ! u-derivatives
            f1aux = de * g0 / ( two * aj(i,j,k)   ) * ret0 * ( four * q(2,i,j+1,k) - q(2,i,j+2,k) ) 
            f2aux = de * g1 / ( two * aj(i,j+1,k) ) * ret1 * q(2,i,j+2,k) 

            fv1star = eta(1,i,j,k) / sqrt(g0) * aj(i,j,k) * ( f2aux - f1aux )

            ! v-derivatives
            f1aux = de * g0 / ( two * aj(i,j,k)   ) * ret0 * ( four * q(3,i,j+1,k) - q(3,i,j+2,k) ) 
            f2aux = de * g1 / ( two * aj(i,j+1,k) ) * ret1 * q(3,i,j+2,k) 

            fv2star = eta(2,i,j,k) / sqrt(g0) * aj(i,j,k) * ( f2aux - f1aux )

            ! pressure extrapolation
            q(1,i,j,k) = c1 * q(1,i,j+1,k) + c2 * q(1,i,j+2,k) + c3 * ( fv1star + fv2star )

         end do
         end do

      case(4)

         j  = j_myend
         c1 = four / three
         c2 = -one / three
         c3 =  two / three

         do k = k_mysta , k_myend
         do i = i_mysta , i_myend

            ! g0 = g22 at the wall (0 nodes away from the wall)
            g0 = eta(1,i,j,k)**two + eta(2,i,j,k)**two 

            ! g1 = g22 1 node away from the wall
            g1 = eta(1,i,j-1,k)**two + eta(2,i,j-1,k)**two 

            ! ret0 = ( 1/Re + xnut ) at the wall
            ret0 = one / ren + xnut(i,j,k) 

            ! ret1 = ( 1/Re + xnut ) 1 node away from the wall
            ret1 = one / ren + xnut(i,j-1,k) 

            ! u-derivatives
            f1aux = de * g1 / ( two * aj(i,j-1,k) ) * ret1 * (-q(2,i,j-2,k) ) 
            f2aux = de * g0 / ( two * aj(i,j,k)   ) * ret0 * (-four * q(2,i,j-1,k) + q(2,i,j-2,k) ) 

            fv1star = eta(1,i,j,k) / sqrt(g0) * aj(i,j,k) * ( f2aux - f1aux )

            ! v-derivatives
            f1aux = de * g1 / ( two * aj(i,j-1,k) ) * ret1 * (-q(3,i,j-2,k) ) 
            f2aux = de * g0 / ( two * aj(i,j,k)   ) * ret0 * (-four * q(3,i,j+1,k) + q(3,i,j-2,k) ) 

            fv2star = eta(2,i,j,k) / sqrt(g0) * aj(i,j,k) * ( f2aux - f1aux )

            ! pressure extrapolation
            q(1,i,j,k) = c1 * q(1,i,j-1,k) + c2 * q(1,i,j-2,k) + c3 * ( fv1star + fv2star )

         end do
         end do

   end select eta_dir

end subroutine p_no_slip_wall_bcond