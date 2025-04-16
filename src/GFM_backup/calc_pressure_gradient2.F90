


subroutine calc_pressure_gradient2(i, j, k, pressure_gradient, exsign)
   
   use AdvectionMethods

   ! compute derivatives of cartesian velocities in curvilinear directions
   ! It uses adaptative stencils depending on whether the node is within the 
   ! water phase or not
   ! velocity_curv_gradient(i,l) = ∂u_i/∂ξ^l

   implicit none

   integer, intent(in) :: i,j,k
   real (kind = rdf), dimension(1:3), intent(out) :: pressure_gradient 
   real (kind = rdf), intent(inout) :: exsign ! flag variable. exsign = 0, the gradient is set to zero
                               !                exsign = 1, the computed gradient is returned

   real (kind = rdf), dimension(1:3) :: pressure_curv_gradient    
   real ( kind = rdf ), dimension(3,3) :: MetricsTensor
   integer :: m,p 

   ! local variables

   real ( kind = rdf ) :: rLL , rL , rC , rR , rRR
   real ( kind = rdf ) :: pLL , pL , pC , pR , pRR

   real ( kind = rdf ) :: ucon_1 , ucon_2 , ucon_3
   real ( kind = rdf ) :: up1 , um1 , up2 , um2 , up3 , um3

   integer :: bias_csi , bias_eta , bias_zet


   pressure_gradient = zero
   
   ! Set the bias of the derivative if I'm at one of the boundaries
   bias_csi = 0
   bias_eta = 0
   bias_zet = 0


   !  Contravariant velocity
   ucon_1 = csi(1,i,j,k) * q(2,i,j,k) + &
            csi(2,i,j,k) * q(3,i,j,k) + &
            csi(3,i,j,k) * q(4,i,j,k)    

   ucon_2 = eta(1,i,j,k) * q(2,i,j,k) + &
            eta(2,i,j,k) * q(3,i,j,k) + &
            eta(3,i,j,k) * q(4,i,j,k)    

   ucon_3 = zet(1,i,j,k) * q(2,i,j,k) + &
            zet(2,i,j,k) * q(3,i,j,k) + &
            zet(3,i,j,k) * q(4,i,j,k)    


   up1 = zero
   um1 = zero

   up2 = zero
   um2 = zero

   up3 = zero
   um3 = zero

   if ( ucon_1 > 0 ) then
      up1 = one
   else
      um1 = one
   end if

   if ( ucon_2 > 0 ) then
      up2 = one
   else
      um2 = one
   end if

   if ( ucon_3 > 0 ) then
      up3 = one
   else
      um3 = one
   end if

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !                                       ξ - direction
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if ( i == ista ) then

      bias_csi = 1

      rLL = zero                    
      rL  = zero                   
      rC  = rsign ( i   , j , k )  
      rR  = rsign ( i+1 , j , k )  
      rRR = rsign ( i+2 , j , k )  

      pLL = zero                
      pL  = zero                
      pC  = q ( 1, i   , j , k )
      pR  = q ( 1, i+1 , j , k )
      pRR = q ( 1, i+2 , j , k )
 
      ! ∂p/∂ξ
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        pLL, pL , pC, pR, pRR        , &
                                        dc                           , &
                                        bias_csi                     , &
                                        exsign                       , &
                                        pressure_curv_gradient(1)      &
                                       )
   
      if ( exsign < one_half ) return                                    
                                                                      

   else if (i == iend) then

      bias_csi = -1

      rLL      = rsign ( i-2 , j , k )
      rL       = rsign ( i-1 , j , k )
      rC       = rsign ( i   , j , k )
      rR       = zero ! dummy
      rRR      = zero ! dummy

      pLL = q ( 1, i-2 , j , k ) 
      pL  = q ( 1, i-1 , j , k ) 
      pC  = q ( 1, i   , j , k ) 
      pR  = zero                 
      pRR = zero                 

      ! ∂u/∂ξ
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        pLL, pL , pC, pR, pRR        , &
                                        dc                           , &
                                        bias_csi                     , &
                                        exsign                       , &
                                        pressure_curv_gradient(1)      &
                                       )
   
      if ( exsign < one_half ) return                                    
                                       
   else

      rL = ( rsign ( i-1 , j , k ) + abs(rsign ( i-1 , j , k )) ) / two
      rC = ( rsign ( i   , j , k ) + abs(rsign ( i   , j , k )) ) / two
      rR = ( rsign ( i+1 , j , k ) + abs(rsign ( i+1 , j , k )) ) / two

      exsign = up1 * ( rC * rL ) + um1 * ( rC * rR )
      
      if ( exsign < one_half ) return

      pressure_curv_gradient(1) = up1 * dc * ( q(1,i,j,k)   - q(1,i-1,j,k)  ) + &
                                  um1 * dc * ( q(1,i+1,j,k) - q(1,i,j,k)    )

   end if


   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !                                       η - direction
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if ( j == jsta ) then

      bias_eta = 1

      rLL = zero                    
      rL  = zero                   
      rC  = rsign ( i, j    , k )  
      rR  = rsign ( i, j+1  , k )  
      rRR = rsign ( i, j+2  , k )  

      pLL = zero                 
      pL  = zero                 
      pC  = q ( 1, i , j   , k ) 
      pR  = q ( 1, i , j+1 , k ) 
      pRR = q ( 1, i , j+2 , k ) 

      ! ∂p/∂η
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        pLL, pL , pC, pR, pRR        , &
                                        de                           , &
                                        bias_eta                     , &
                                        exsign                       , &
                                        pressure_curv_gradient(2)      &
                                       )
   
      if ( exsign < one_half ) return                                    
                                       
   
   else if ( j == jend ) then

      bias_eta = -1

      rLL      = rsign ( i , j-2 , k )
      rL       = rsign ( i , j-1 , k )
      rC       = rsign ( i , j   , k )
      rR       = zero ! dummy
      rRR      = zero ! dummy

      pLL = q ( 1, i , j-2 , k ) 
      pL  = q ( 1, i , j-1 , k ) 
      pC  = q ( 1, i , j   , k ) 
      pR  = zero                 
      pRR = zero                 

      ! ∂p/∂η
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        pLL, pL , pC, pR, pRR        , &
                                        de                           , &
                                        bias_eta                     , &
                                        exsign                       , &
                                        pressure_curv_gradient(2)      &
                                       )
   
      if ( exsign < one_half ) return                                    


   else

      rL = ( rsign ( i , j-1 , k ) + abs(rsign ( i , j-1 , k )) ) / two
      rC = ( rsign ( i , j   , k ) + abs(rsign ( i , j   , k )) ) / two
      rR = ( rsign ( i , j+1 , k ) + abs(rsign ( i , j+1 , k )) ) / two

      exsign = up2 * ( rC * rL ) + um2 * ( rC * rR )
      
      if ( exsign < one_half ) return

      pressure_curv_gradient(2) = up2 * de * ( q(1,i,j,k)   - q(1,i,j-1,k)  ) + &
                                  um2 * de * ( q(1,i,j+1,k) - q(1,i,j,k)    )

   end if


   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !                                       ζ - direction
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if ( k == ksta ) then

      bias_zet = 1

      rLL = zero                    
      rL  = zero                   
      rC  = rsign ( i , j , k   )  
      rR  = rsign ( i , j , k+1 )  
      rRR = rsign ( i , j , k+2 )  

      pLL = zero                 
      pL  = zero                 
      pC  = q ( 1, i , j , k   ) 
      pR  = q ( 1, i , j , k+1 ) 
      pRR = q ( 1, i , j , k+2 ) 
 
      ! ∂p/∂ζ
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        pLL, pL , pC, pR, pRR        , &
                                        dz                           , &
                                        bias_zet                     , &
                                        exsign                       , &
                                        pressure_curv_gradient(3)      &
                                       )
   
      if ( exsign < one_half ) return                                    
                                       
   else if ( k == kend ) then

      bias_zet = -1

      rLL      = rsign ( i , j , k-2 )
      rL       = rsign ( i , j , k-1 )
      rC       = rsign ( i , j , k   )
      rR       = zero ! dummy
      rRR      = zero ! dummy

      pLL = q ( 1, i , j , k-2 ) 
      pL  = q ( 1, i , j , k-1 ) 
      pC  = q ( 1, i , j , k   ) 
      pR  = zero                 
      pRR = zero                 

      ! ∂p/∂ζ
      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                        pLL, pL , pC, pR, pRR        , &
                                        dz                           , &
                                        bias_zet                     , &
                                        exsign                       , &
                                        pressure_curv_gradient(3)      &
                                       )
   
      if ( exsign < one_half ) return                                    

   else

      rL = ( rsign ( i , j , k-1 ) + abs(rsign ( i , j , k-1 )) ) / two
      rC = ( rsign ( i , j , k   ) + abs(rsign ( i , j , k   )) ) / two
      rR = ( rsign ( i , j , k+1 ) + abs(rsign ( i , j , k+1 )) ) / two

      exsign = up3 * ( rC * rL ) + um3 * ( rC * rR )
      
      if ( exsign < one_half ) return

      pressure_curv_gradient(3) = up3 * dz * ( q(1,i,j,k)   - q(1,i,j,k-1)  ) + &
                                  um3 * dz * ( q(1,i,j,k+1) - q(1,i,j,k)  )

   end if

   MetricsTensor(1,1) = csi(1,i,j,k)
   MetricsTensor(1,2) = csi(2,i,j,k)
   MetricsTensor(1,3) = csi(3,i,j,k)

   MetricsTensor(2,1) = eta(1,i,j,k)
   MetricsTensor(2,2) = eta(2,i,j,k)
   MetricsTensor(2,3) = eta(3,i,j,k)

   MetricsTensor(3,1) = zet(1,i,j,k)
   MetricsTensor(3,2) = zet(2,i,j,k)
   MetricsTensor(3,3) = zet(3,i,j,k)

   ! PressureGradient(m) = ∂p/∂xm = ( ∂p/∂ξp ) * ( ∂ξp/∂xm )

   do m = 1,3
   do p = 1,3
       
      pressure_gradient(m) =   pressure_gradient(m) &
                             + pressure_curv_gradient(p) * MetricsTensor(p,m)
   end do
   end do



end subroutine calc_pressure_gradient2