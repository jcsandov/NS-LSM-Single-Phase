subroutine ghost_nodes_extrapolation_old(i,j,k, xs, ys, zs, u_fs_lsqm, du_dx_fs_lsqm, dp_dx_fs_lsqm)

   implicit none

   integer, intent(in) :: i,j,k
   real (kind = rdf), dimension(1:3), intent(in) :: u_fs_lsqm
   real (kind = rdf), dimension(1:3,1:3), intent(in) :: du_dx_fs_lsqm
   real (kind = rdf), dimension(1:3), intent(in) :: dp_dx_fs_lsqm
   real (kind = rdf), intent(in) :: xs, ys, zs ! free-surface position closest to (i,j,k) node

   ! local variables
   real (kind = rdf), dimension(1:3) :: alpha_local
   real (kind = rdf) :: uextp, vextp, wextp ! extrapolated velocity using Taylor expansion
   real (kind = rdf) :: pextp
   real (kind = rdf) :: rdiff_norm, exsign  

   integer :: ii, jj, kk, extp_swept

   integer :: i_mysta, &
              j_mysta, &
              k_mysta, &
              i_myend, &
              j_myend, &
              k_myend

  !Interior nodes including domain boundaries
   i_mysta = il + igp
   j_mysta = jl + jgp
   k_mysta = kl + kgp
   
   i_myend = iu - igp
   j_myend = ju - jgp
   k_myend = ku - kgp

   ! extp_swept = # of nodes around i,j,k to be candidates for extrapolation
   extp_swept = 3 ! TO DO: it could be set from control.dat

   ! max(k_mysta-1,k-extp_swept),min(k_myend+1,k+extp_swept) means that it even extrapolates
   ! the velocity to the air nodes at the boundaries of the domain. This is gonna be useful to
   ! compute velocity gradients in the air phase for the normal boundary condition.

   !do kk = max( k_mysta , k-extp_swept ) , min( k_myend , k+extp_swept ) 
   !do jj = max( j_mysta , j-extp_swept ) , min( j_myend , j+extp_swept )
   !do ii = max( i_mysta , i-extp_swept ) , min( i_myend , i+extp_swept )

   !do kk = max( kl , k-extp_swept ) , min( ku , k+extp_swept ) 
   !do jj = max( jl , j-extp_swept ) , min( ju , j+extp_swept )
   !do ii = max( il , i-extp_swept ) , min( iu , i+extp_swept )

   do kk = max( k_mysta , k-extp_swept ) , min( k_myend , k+extp_swept ) 
   do jj = max( j_mysta , j-extp_swept ) , min( j_myend , j+extp_swept )
   do ii = max( i_mysta , i-extp_swept ) , min( i_myend , i+extp_swept )

      ! if the (ii,jj,kk) is in the air-phase, it's extrapolated using
      ! a Taylor expansion
      
      if ( rsign(ii,jj,kk) < one_half ) then ! (air-phase)

         ! position vector from the free-surface to the extrapolated node 

         alpha_local = zero

         alpha_local(1) = x(ii,jj,kk) - xs
         alpha_local(2) = y(ii,jj,kk) - ys
         alpha_local(3) = z(ii,jj,kk) - zs

         ! distance between the extrapolation candidate (ii,jj,kk) node and the
         ! free - surface
         rdiff_norm = norm2( alpha_local )

         ! exsign: flag variable
         !
         ! exsign = 1 ,if current distance between extrapolation candidate and free 
         ! surface is less than a previous stored one (previous calls of the 
         ! ghost_nodes_velocity_extrapolation subroutine).
         !
         ! exsign = 0 otherwise

         exsign =  (sign(one , least_dis_extp(ii,jj,kk) - rdiff_norm) + one)/two
         
         ! if exsign = 1, the new least distance is stored in least_dis_extp(ii,jj,kk)
         least_dis_extp(ii,jj,kk) =            least_dis_extp(ii,jj,kk) + & 
                                    exsign * (-least_dis_extp(ii,jj,kk) + rdiff_norm)

         ! pressure extrapolation assuming pfs = 0
         pextp = zero

         pextp = dp_dx_fs_lsqm(1) * alpha_local(1) + &
                 dp_dx_fs_lsqm(2) * alpha_local(2) + &
                 dp_dx_fs_lsqm(3) * alpha_local(3)

         ! extrapolated velocities using Taylor expansion
         uextp = zero; vextp = zero; wextp = zero;

         uextp = u_fs_lsqm(1) +  du_dx_fs_lsqm(1,1) * alpha_local(1) + & 
                                 du_dx_fs_lsqm(1,2) * alpha_local(2) + & 
                                 du_dx_fs_lsqm(1,3) * alpha_local(3)

         vextp = u_fs_lsqm(2) +  du_dx_fs_lsqm(2,1) * alpha_local(1) + & 
                                 du_dx_fs_lsqm(2,2) * alpha_local(2) + & 
                                 du_dx_fs_lsqm(2,3) * alpha_local(3)

         wextp = u_fs_lsqm(3) +  du_dx_fs_lsqm(3,1) * alpha_local(1) + & 
                                 du_dx_fs_lsqm(3,2) * alpha_local(2) + & 
                                 du_dx_fs_lsqm(3,3) * alpha_local(3)

         ! if exsign = 1 (least distanced node), the velocity at ii,jj,kk
         ! is replaced by the extrapolated one

         if ( exsign > one_half ) then
      
            ! if rsign = -1, then I keep the geometrically extrapolated value
            ! if rsign =  0, then I use the LSQM extrapolation
            q(1,ii,jj,kk) =   q(1,ii,jj,kk)  &
                            + ( one + rsign(ii,jj,kk) ) * exsign * ( - q(1,ii,jj,kk) + pextp )

            ! Different from the one above, I will keep the LSQM extrapolation even for the 
            ! nondes next to the free surface
            !q(1,ii,jj,kk) = q(1,ii,jj,kk) + exsign * ( - q(1,ii,jj,kk) + pextp )
            q(2,ii,jj,kk) = q(2,ii,jj,kk) + exsign * ( - q(2,ii,jj,kk) + uextp )
            q(3,ii,jj,kk) = q(3,ii,jj,kk) + exsign * ( - q(3,ii,jj,kk) + vextp )
            q(4,ii,jj,kk) = q(4,ii,jj,kk) + exsign * ( - q(4,ii,jj,kk) + wextp )

         end if
      
      end if

   end do
   end do
   end do

end subroutine ghost_nodes_extrapolation_old