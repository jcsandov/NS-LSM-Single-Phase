subroutine free_surface_pressure_gradient( i , j , k , xs, ys, zs, pfs , nvec , pgrad_fs )

   !--------------------------------------------------------------------------
   ! Description: 
   ! ------------
   ! 
   ! Calculates the pressure gradient at the free surface (xs,ys,zs) using a 
   ! Weighted Least Square Method (WLSM). 
   ! 
   ! the weights are a combination of distance and alignment with the normal 
   ! direction, such that nodes closer and/or more normal to the free surface
   ! will contribute more to the computation of the gradient
   !
   !--------------------------------------------------------------------------

   use InterpolationMethods 

   implicit none

   ! Input
   integer , intent(in) :: i,j,k ! Node in the air phase to have the pressure extrapolated
   real(kind = rdf) , intent(in) :: xs, ys , zs ! Free surface location
   real(kind = rdf) , intent(in) :: pfs ! Pressure at the free surface obtained from the NDBC
   real(kind = rdf), dimension(1:3), intent(in) :: nvec ! normal vector at the free surface

   ! Output
   real(kind = rdf), dimension(1:3), intent(out) :: pgrad_fs ! pressure gradient at the free surface

   ! Local variables

   integer :: ii,jj,kk
   integer :: inst, jnst, knst ! i,j,k nearest
   real(kind = rdf), dimension(1:3,1:3) :: Gmatrix , Ginv  ! M-matrix WLSM
   real(kind = rdf), dimension(1:3)     :: bwvec ! rhs vector
   real(kind = rdf) :: wlocal ! local weight
   real(kind = rdf) :: dx , dy , dz , dp 
   real(kind = rdf), dimension(1:3) :: dr ! Δr = rii - rfs
   real(kind = rdf) :: num , denom 
   real(kind = rdf) :: aux 
   real(kind = rdf) :: exsign
   real(kind = rdf) :: dsmin 
   real(kind=rdf), parameter :: eps_denom = 1.0E-8
   real(kind = rdf) :: dp_dx , dp_dy , dp_dz 
   
   real (kind = rdf) :: max_pgrad_norm_neighbour   ! 
   real(kind = rdf), dimension(1:3) :: pgrad_nearest_point
   integer :: counter

   logical :: BlankingFlag
   logical :: OK_FLAG

   Gmatrix  = zero
   Ginv     = zero
   bwvec    = zero
   
   pgrad_fs = zero

   ! Initialisation
   dsmin = ten * ten

   ! Loop around i,j,k to determine the smallest grid distance dsmin with
   ! its neighbours
   !do kk = max( ksta , k-1 ), min( kend , k+1 )
   !do jj = max( jsta , j-1 ), min( jend , j+1 )
   !do ii = max( ista , i-1 ), min( iend , i+1 )

   !   if ( rsign( phi(ii,jj,kk) ) < one_half ) cycle
   !
   !   BlankingFlag = .false.
   !
   !   if ( nblke /= 0 ) then
   !
   !      do nb = 1,nblke
   !         if ( ii > le_blk_ia(1,nb) .and. ii < le_blk_ib(1,nb) .and. &  
   !              jj > le_blk_ja(1,nb) .and. jj < le_blk_jb(1,nb) .and. &
   !              kk > le_blk_ka(1,nb) .and. kk < le_blk_kb(1,nb) ) then
   !            
   !            BlankingFlag = .true.
   !
   !         end if
   !      end do
   !   end if
   !
   !   if ( BlankingFlag ) cycle
   !
   !   dx = x(i,j,k) - x(ii,jj,kk)
   !   dy = y(i,j,k) - y(ii,jj,kk)
   !   dz = z(i,j,k) - z(ii,jj,kk)
   !
   !  dsmin = minval( (/ dsmin , abs(dx) , abs(dy) , abs(dz) /) )
   !
   !end do
   !end do
   !end do

   ! Gradient limiting parameters
   max_pgrad_norm_neighbour = zero
   !pgrad_fs_avg             = zero
   counter                  = 0 ! counter of valid water nodes to calculate the average gradient

   do kk = max( ksta , k-sweep_lsqm ), min( kend , k+sweep_lsqm )
   do jj = max( jsta , j-sweep_lsqm ), min( jend , j+sweep_lsqm )
   do ii = max( ista , i-sweep_lsqm ), min( iend , i+sweep_lsqm )

      if ( rsign( phi(ii,jj,kk) ) < one_half ) cycle

      BlankingFlag = .false.

      if ( nblke /= 0 ) then
         do nb = 1,nblke
            if ( ii >= le_blk_ia(1,nb) .and. ii <= le_blk_ib(1,nb) .and. &  
                 jj >= le_blk_ja(1,nb) .and. jj <= le_blk_jb(1,nb) .and. &
                 kk >= le_blk_ka(1,nb) .and. kk <= le_blk_kb(1,nb) ) then
               
               BlankingFlag = .true.

            end if
         end do
      end if

      if ( BlankingFlag ) cycle

      dx = xs - x(ii,jj,kk)
      dy = ys - y(ii,jj,kk)
      dz = zs - z(ii,jj,kk)

      dr = (/dx,dy,dz/)

      ! Set the nearest node
      if ( norm2(dr) < dsmin ) then

         dsmin = norm2(dr)
         inst = ii
         jnst = jj
         knst = kk

      end if

      ! If the current water node is too close to the free surface
      ! then I don't use it
      !if ( norm2(dr) < 0.1 * dsmin ) cycle

      dp = pfs - q(1,ii,jj,kk) 

      !if ( abs(dx) > eps_sims .and. &
      !     abs(dy) > eps_sims .and. &
      !     abs(dz) > eps_sims ) then

      !   max_pgrad_norm_neighbour = max ( max_pgrad_norm_neighbour , &
      !                                    norm2( (/ dp/dx , dp/dy , dp/dz /) ) )

      !   pgrad_fs_avg = pgrad_fs_avg +  (/ dp/dx , dp/dy , dp/dz /)
      !   counter      = counter      + 1
      !
      !end if
      
      num   = abs( dot_product( dr/norm2(dr) , nvec ) )**(three/two) 
      denom = one !norm2(dr)**(three/two) + eps_denom

      wlocal = num / denom

      Gmatrix(1,1) = Gmatrix(1,1) + ( dx * dx ) * ( wlocal**two ) 
      Gmatrix(1,2) = Gmatrix(1,2) + ( dx * dy ) * ( wlocal**two ) 
      Gmatrix(1,3) = Gmatrix(1,3) + ( dx * dz ) * ( wlocal**two ) 

      Gmatrix(2,1) = Gmatrix(2,1) + ( dy * dx ) * ( wlocal**two ) 
      Gmatrix(2,2) = Gmatrix(2,2) + ( dy * dy ) * ( wlocal**two ) 
      Gmatrix(2,3) = Gmatrix(2,3) + ( dy * dz ) * ( wlocal**two ) 

      Gmatrix(3,1) = Gmatrix(3,1) + ( dz * dx ) * ( wlocal**two ) 
      Gmatrix(3,2) = Gmatrix(3,2) + ( dz * dy ) * ( wlocal**two ) 
      Gmatrix(3,3) = Gmatrix(3,3) + ( dz * dz ) * ( wlocal**two ) 

      bwvec(1)     = bwvec(1)     + ( dx * dp ) * ( wlocal**two )
      bwvec(2)     = bwvec(2)     + ( dy * dp ) * ( wlocal**two )
      bwvec(3)     = bwvec(3)     + ( dz * dp ) * ( wlocal**two )

      counter = counter + 1

   end do
   end do
   end do

   ! Gmatrix inversion
   call M33INV( Gmatrix , Ginv , OK_FLAG )

   if ( OK_FLAG ) then 

      pgrad_fs = matmul( Ginv , bwvec )   
   
   else

      print *, 'problem inverting the G matrix in WLSM for pressure at node: ',myid,i,j,k
      print *, 'counter = : ', counter
      print *, ' '

      call pressure_curv_gradient_vector(inst, jnst, knst, pgrad_nearest_point, exsign)

      ! ∂p      ∂p   ∂ξ     ∂p   ∂η     ∂p   ∂ζ 
      !---- =  ---- ---- + ---- ---- + ---- ----   
      ! ∂x      ∂ξ   ∂x     ∂η   ∂x     ∂ζ   ∂x 
      pgrad_fs(1) = pgrad_nearest_point(1) * csi(1,inst,jnst,knst) + &
                    pgrad_nearest_point(2) * eta(1,inst,jnst,knst) + &
                    pgrad_nearest_point(3) * zet(1,inst,jnst,knst)
    
      ! ∂p      ∂p   ∂ξ     ∂p   ∂η     ∂p   ∂ζ 
      !---- =  ---- ---- + ---- ---- + ---- ----   
      ! ∂y      ∂ξ   ∂y     ∂η   ∂y     ∂ζ   ∂y 
      pgrad_fs(2) = pgrad_nearest_point(1) * csi(2,inst,jnst,knst) + &
                    pgrad_nearest_point(2) * eta(2,inst,jnst,knst) + &
                    pgrad_nearest_point(3) * zet(2,inst,jnst,knst)
    
      ! ∂p      ∂p   ∂ξ     ∂p   ∂η     ∂p   ∂ζ 
      !---- =  ---- ---- + ---- ---- + ---- ----   
      ! ∂z      ∂ξ   ∂z     ∂η   ∂z     ∂ζ   ∂z 
      pgrad_fs(3) = pgrad_nearest_point(1) * csi(3,inst,jnst,knst) + &
                    pgrad_nearest_point(2) * eta(3,inst,jnst,knst) + &
                    pgrad_nearest_point(3) * zet(3,inst,jnst,knst)
    
   end if

   !if ( norm2( pgrad_fs ) > max_pgrad_norm_neighbour .or. .not. OK_FLAG ) then
   !
   !   print *, 'pgrad_fs = '                   , pgrad_fs       , &
   !            'applied pavg at myid,i,j,k = ' , myid,i,j,k     , &
   !            'pgrad_avg = '                  , pgrad_fs_avg / real((counter-1),kind=rdf)
   !
   !   if ( counter > 1 ) pgrad_fs = pgrad_fs_avg / real( (counter-1) , kind=rdf )
   !
   !end if

end subroutine free_surface_pressure_gradient







