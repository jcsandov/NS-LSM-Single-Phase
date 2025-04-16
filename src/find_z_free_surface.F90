subroutine find_z_free_surface ( z_array                 , &  
                                 phi_array               , &  
                                 FreeSurfaceWithinProc   , & 
                                 z_free_surface            &
                                )

   real (kind = rdf), dimension(ksta:kend), intent(in) :: z_array , phi_array
   logical                                             :: FreeSurfaceWithinProc
   real (kind = rdf), intent(out)                      :: z_free_surface

   integer :: k

   FreeSurfaceWithinProc = .false.
   z_free_surface        = zero

   ! Linear search for the interval where y changes sign
   do k = ksta, kend-1
   
      ! The free surface is between k and k+1  
      if ( phi_array(k) * phi_array(k+1) <= zero ) then
   
         ! CASE 1: Free surface at node k
         if ( abs( phi_array(k) )   < eps_sims .and. &
              abs( phi_array(k+1) ) > eps_sims               ) then
   
            z_free_surface = z_array(k)
            FreeSurfaceWithinProc = .true.
            return
   
         ! CASE 2: Free surface at node k+1
         else if ( abs( phi_array(k+1) ) < eps_sims .and. &
                   abs( phi_array(k) )   > eps_sims          ) then   
   
            z_free_surface = z_array(k+1)
            FreeSurfaceWithinProc = .true.
            return
   
         ! CASE 3: Free surface almost vertical, I set it in the middle 
         else if ( abs( phi_array(k+1) ) < eps_sims .and. &
                   abs( phi_array(k) )   < eps_sims          ) then
   
            z_free_surface = pt5 * ( z_array(k+1) + z_array(k) )
            FreeSurfaceWithinProc = .true.
            return
   
         ! CASE 4 (default): A linearly interpolated value between z(k) and z(k+1)              
         else 
               
            z_free_surface =   z_array(k) +                                       &
                             ( z_array(k+1) - z_array(k) ) * ( -phi_array(k) ) /  &
                             ( phi_array(k+1) - phi_array(k) )
               
            FreeSurfaceWithinProc = .true.
            return
   
         end if ! cases for the position of the free-surface
   
      end if ! The free surface is between k and k+1

   end do ! k = ksta, kend-1

end subroutine find_z_free_surface