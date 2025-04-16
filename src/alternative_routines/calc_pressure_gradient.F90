subroutine calc_pressure_gradient( i, j, k, pressure_gradient , exsign )
   
   use AdvectionMethods

   ! compute derivatives of cartesian velocities in curvilinear directions
   ! It uses adaptative stencils depending on whether the node is within the 
   ! water phase or not

   implicit none

   integer, intent(in) :: i,j,k
   real (kind = rdf) , dimension(1:3) , intent(out) :: pressure_gradient
   real (kind = rdf), intent(inout) :: exsign ! flag variable. exsign = 0, the gradient is set to zero
                               !                exsign = 1, the computed gradient is returned

   ! local variables

   real ( kind = rdf ) :: rLL , rL , rC , rR , rRR
   real ( kind = rdf ) :: pLL , pL , pC , pR , pRR
   real ( kind = rdf ) , dimension(3) :: pressure_curv_gradient 

   integer :: bias_csi , bias_eta , bias_zet

   ! Set the bias of the derivative if I'm at one of the boundaries
   bias_csi = 0
   bias_eta = 0
   bias_zet = 0

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

      pLL = zero                   !;  vLL = zero                  ;  wLL = zero                 ;
      pL  = zero                   !;  vL  = zero                  ;  wL  = zero                 ;
      pC  = q ( 1 , i   , j , k )  !;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
      pR  = q ( 1 , i+1 , j , k )  !;  vR  = q ( 3, i+1 , j , k )  ;  wR  = q ( 4, i+1 , j , k ) ;
      pRR = q ( 1 , i+2 , j , k )  !;  vRR = q ( 3, i+2 , j , k )  ;  wRR = q ( 4, i+2 , j , k ) ;
      
   else if ( i == ista + 1 ) then
         
      rLL      = zero ! dummy
      rL       = rsign ( i-1 , j , k )
      rC       = rsign ( i   , j , k )
      rR       = rsign ( i+1 , j , k )
      rRR      = rsign ( i+2 , j , k )

      pLL = zero                   !;  vLL = zero                  ;  wLL = zero                 ;
      pL  = q ( 1 , i-1 , j , k )  !;  vL  = q ( 3, i-1 , j , k )  ;  wL  = q ( 4, i-1 , j , k ) ;
      pC  = q ( 1 , i   , j , k )  !;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
      pR  = q ( 1 , i+1 , j , k )  !;  vR  = q ( 3, i+1 , j , k )  ;  wR  = q ( 4, i+1 , j , k ) ;
      pRR = q ( 1 , i+2 , j , k )  !;  vRR = q ( 3, i+2 , j , k )  ;  wRR = q ( 4, i+2 , j , k ) ;

   else if ( i == iend - 1 ) then
         
      rLL      = rsign ( i-2 , j , k )
      rL       = rsign ( i-1 , j , k )
      rC       = rsign ( i   , j , k )
      rR       = rsign ( i+1 , j , k )
      rRR      = zero ! dummy

      pLL = q ( 1 , i-2 , j , k )  !;  vLL = q ( 3, i-2 , j , k )  ;  wLL = q ( 4, i-2 , j , k ) ;
      pL  = q ( 1 , i-1 , j , k )  !;  vL  = q ( 3, i-1 , j , k )  ;  wL  = q ( 4, i-1 , j , k ) ;
      pC  = q ( 1 , i   , j , k )  !;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
      pR  = q ( 1 , i+1 , j , k )  !;  vR  = q ( 3, i+1 , j , k )  ;  wR  = q ( 4, i+1 , j , k ) ;
      pRR = zero                   !;  vRR = zero                  ;  wRR = zero                 ;


   else if ( i == iend ) then
         
      bias_csi = -1

      rLL      = rsign ( i-2 , j , k )
      rL       = rsign ( i-1 , j , k )
      rC       = rsign ( i   , j , k )
      rR       = zero ! dummy
      rRR      = zero ! dummy

      pLL = q ( 1 , i-2 , j , k )  !;  vLL = q ( 3, i-2 , j , k )  ;  wLL = q ( 4, i-2 , j , k ) ;
      pL  = q ( 1 , i-1 , j , k )  !;  vL  = q ( 3, i-1 , j , k )  ;  wL  = q ( 4, i-1 , j , k ) ;
      pC  = q ( 1 , i   , j , k )  !;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
      pR  = zero                   !;  vR  = zero                  ;  wR  = zero                 ;
      pRR = zero                   !;  vRR = zero                  ;  wRR = zero                 ;

   else

      rLL      = rsign ( i-2 , j , k )
      rL       = rsign ( i-1 , j , k )
      rC       = rsign ( i   , j , k )
      rR       = rsign ( i+1 , j , k )
      rRR      = rsign ( i+2 , j , k )

      pLL = q ( 1 , i-2 , j , k )  !;  vLL = q ( 3, i-2 , j , k )  ;  wLL = q ( 4, i-2 , j , k ) ;
      pL  = q ( 1 , i-1 , j , k )  !;  vL  = q ( 3, i-1 , j , k )  ;  wL  = q ( 4, i-1 , j , k ) ;
      pC  = q ( 1 , i   , j , k )  !;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
      pR  = q ( 1 , i+1 , j , k )  !;  vR  = q ( 3, i+1 , j , k )  ;  wR  = q ( 4, i+1 , j , k ) ;
      pRR = q ( 1 , i+2 , j , k )  !;  vRR = q ( 3, i+2 , j , k )  ;  wRR = q ( 4, i+2 , j , k ) ;

   end if

   exsign = zero

   ! ∂p/∂ξ
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     pLL, pL , pC, pR, pRR        , &
                                     dc                           , &
                                     bias_csi                     , &
                                     exsign                       , &
                                     pressure_curv_gradient(1)      &
                                    )

   if ( exsign < one_half ) return                                    

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !                                       η - direction
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if ( j == jsta ) then

      bias_eta = 1

      rLL = zero                    
      rL  = zero                   
      rC  = rsign ( i , j   , k )  
      rR  = rsign ( i , j+1 , k )  
      rRR = rsign ( i , j+2 , k )  

      pLL = zero                   !;  vLL = zero                  ;  wLL = zero                 ;
      pL  = zero                   !;  vL  = zero                  ;  wL  = zero                 ;
      pC  = q ( 1 , i , j   , k )  !;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
      pR  = q ( 1 , i , j+1 , k )  !;  vR  = q ( 3, i+1 , j , k )  ;  wR  = q ( 4, i+1 , j , k ) ;
      pRR = q ( 1 , i , j+2 , k )  !;  vRR = q ( 3, i+2 , j , k )  ;  wRR = q ( 4, i+2 , j , k ) ;
      
   else if ( j == jsta + 1 ) then
         
      rLL      = zero ! dummy
      rL       = rsign ( i , j-1 , k )
      rC       = rsign ( i , j   , k )
      rR       = rsign ( i , j+1 , k )
      rRR      = rsign ( i , j+2 , k )

      pLL = zero                   !;  vLL = zero                  ;  wLL = zero                 ;
      pL  = q ( 1 , i , j-1 , k )  !;  vL  = q ( 3, i-1 , j , k )  ;  wL  = q ( 4, i-1 , j , k ) ;
      pC  = q ( 1 , i , j   , k )  !;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
      pR  = q ( 1 , i , j+1 , k )  !;  vR  = q ( 3, i+1 , j , k )  ;  wR  = q ( 4, i+1 , j , k ) ;
      pRR = q ( 1 , i , j+2 , k )  !;  vRR = q ( 3, i+2 , j , k )  ;  wRR = q ( 4, i+2 , j , k ) ;

   else if ( j == jend - 1 ) then
         
      rLL      = rsign ( i , j-2 , k )
      rL       = rsign ( i , j-1 , k )
      rC       = rsign ( i , j   , k )
      rR       = rsign ( i , j+1 , k )
      rRR      = zero ! dummy

      pLL = q ( 1 , i , j-2 , k )  !;  vLL = q ( 3, i-2 , j , k )  ;  wLL = q ( 4, i-2 , j , k ) ;
      pL  = q ( 1 , i , j-1 , k )  !;  vL  = q ( 3, i-1 , j , k )  ;  wL  = q ( 4, i-1 , j , k ) ;
      pC  = q ( 1 , i , j   , k )  !;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
      pR  = q ( 1 , i , j+1 , k )  !;  vR  = q ( 3, i+1 , j , k )  ;  wR  = q ( 4, i+1 , j , k ) ;
      pRR = zero                   !;  vRR = zero                  ;  wRR = zero                 ;


   else if ( j == jend ) then
         
      bias_eta = -1

      rLL      = rsign ( i , j-2 , k )
      rL       = rsign ( i , j-1 , k )
      rC       = rsign ( i , j   , k )
      rR       = zero ! dummy
      rRR      = zero ! dummy

      pLL = q ( 1 , i , j-2 , k )  !;  vLL = q ( 3, i-2 , j , k )  ;  wLL = q ( 4, i-2 , j , k ) ;
      pL  = q ( 1 , i , j-1 , k )  !;  vL  = q ( 3, i-1 , j , k )  ;  wL  = q ( 4, i-1 , j , k ) ;
      pC  = q ( 1 , i , j   , k )  !;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
      pR  = zero                   !;  vR  = zero                  ;  wR  = zero                 ;
      pRR = zero                   !;  vRR = zero                  ;  wRR = zero                 ;

   else

      rLL      = rsign ( i , j-2 , k )
      rL       = rsign ( i , j-1 , k )
      rC       = rsign ( i , j   , k )
      rR       = rsign ( i , j+1 , k )
      rRR      = rsign ( i , j+2 , k )

      pLL = q ( 1 , i , j-2 , k )  !;  vLL = q ( 3, i-2 , j , k )  ;  wLL = q ( 4, i-2 , j , k ) ;
      pL  = q ( 1 , i , j-1 , k )  !;  vL  = q ( 3, i-1 , j , k )  ;  wL  = q ( 4, i-1 , j , k ) ;
      pC  = q ( 1 , i , j   , k )  !;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
      pR  = q ( 1 , i , j+1 , k )  !;  vR  = q ( 3, i+1 , j , k )  ;  wR  = q ( 4, i+1 , j , k ) ;
      pRR = q ( 1 , i , j+2 , k )  !;  vRR = q ( 3, i+2 , j , k )  ;  wRR = q ( 4, i+2 , j , k ) ;

   end if

   exsign = zero

   ! ∂p/∂η
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     pLL, pL , pC, pR, pRR        , &
                                     de                           , &
                                     bias_eta                     , &
                                     exsign                       , &
                                     pressure_curv_gradient(2)      &
                                    )

   if ( exsign < one_half ) return                                    

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

      pLL = zero                   !;  vLL = zero                  ;  wLL = zero                 ;
      pL  = zero                   !;  vL  = zero                  ;  wL  = zero                 ;
      pC  = q ( 1 , i , j , k   )  !;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
      pR  = q ( 1 , i , j , k+1 )  !;  vR  = q ( 3, i+1 , j , k )  ;  wR  = q ( 4, i+1 , j , k ) ;
      pRR = q ( 1 , i , j , k+2 )  !;  vRR = q ( 3, i+2 , j , k )  ;  wRR = q ( 4, i+2 , j , k ) ;
      
   else if ( k == ksta + 1 ) then
         
      rLL      = zero ! dummy
      rL       = rsign ( i , j , k-1 )
      rC       = rsign ( i , j , k   )
      rR       = rsign ( i , j , k+1 )
      rRR      = rsign ( i , j , k+2 )

      pLL = zero                   !;  vLL = zero                  ;  wLL = zero                 ;
      pL  = q ( 1 , i , j , k-1 )  !;  vL  = q ( 3, i-1 , j , k )  ;  wL  = q ( 4, i-1 , j , k ) ;
      pC  = q ( 1 , i , j , k   )  !;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
      pR  = q ( 1 , i , j , k+1 )  !;  vR  = q ( 3, i+1 , j , k )  ;  wR  = q ( 4, i+1 , j , k ) ;
      pRR = q ( 1 , i , j , k+2 )  !;  vRR = q ( 3, i+2 , j , k )  ;  wRR = q ( 4, i+2 , j , k ) ;

   else if ( k == kend - 1 ) then
         
      rLL      = rsign ( i , j , k-2 )
      rL       = rsign ( i , j , k-1 )
      rC       = rsign ( i , j , k   )
      rR       = rsign ( i , j , k+1 )
      rRR      = zero ! dummy

      pLL = q ( 1 , i , j , k-2 )  !;  vLL = q ( 3, i-2 , j , k )  ;  wLL = q ( 4, i-2 , j , k ) ;
      pL  = q ( 1 , i , j , k-1 )  !;  vL  = q ( 3, i-1 , j , k )  ;  wL  = q ( 4, i-1 , j , k ) ;
      pC  = q ( 1 , i , j , k   )  !;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
      pR  = q ( 1 , i , j , k+1 )  !;  vR  = q ( 3, i+1 , j , k )  ;  wR  = q ( 4, i+1 , j , k ) ;
      pRR = zero                   !;  vRR = zero                  ;  wRR = zero                 ;


   else if ( k == kend ) then
         
      bias_zet = -1

      rLL      = rsign ( i , j , k-2 )
      rL       = rsign ( i , j , k-1 )
      rC       = rsign ( i , j , k   )
      rR       = zero ! dummy
      rRR      = zero ! dummy

      pLL = q ( 1 , i , j , k-2 )  !;  vLL = q ( 3, i-2 , j , k )  ;  wLL = q ( 4, i-2 , j , k ) ;
      pL  = q ( 1 , i , j , k-1 )  !;  vL  = q ( 3, i-1 , j , k )  ;  wL  = q ( 4, i-1 , j , k ) ;
      pC  = q ( 1 , i , j , k   )  !;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
      pR  = zero                   !;  vR  = zero                  ;  wR  = zero                 ;
      pRR = zero                   !;  vRR = zero                  ;  wRR = zero                 ;

   else

      rLL      = rsign ( i , j , k-2 )
      rL       = rsign ( i , j , k-1 )
      rC       = rsign ( i , j , k   )
      rR       = rsign ( i , j , k+1 )
      rRR      = rsign ( i , j , k+2 )

      pLL = q ( 1 , i , j , k-2 )  !;  vLL = q ( 3, i-2 , j , k )  ;  wLL = q ( 4, i-2 , j , k ) ;
      pL  = q ( 1 , i , j , k-1 )  !;  vL  = q ( 3, i-1 , j , k )  ;  wL  = q ( 4, i-1 , j , k ) ;
      pC  = q ( 1 , i , j , k   )  !;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
      pR  = q ( 1 , i , j , k+1 )  !;  vR  = q ( 3, i+1 , j , k )  ;  wR  = q ( 4, i+1 , j , k ) ;
      pRR = q ( 1 , i , j , k+2 )  !;  vRR = q ( 3, i+2 , j , k )  ;  wRR = q ( 4, i+2 , j , k ) ;

   end if

   exsign = zero

   ! ∂p/∂ζ
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     pLL, pL , pC, pR, pRR        , &
                                     dz                           , &
                                     bias_zet                     , &
                                     exsign                       , &
                                     pressure_curv_gradient(3)      &
                                    )

   if ( exsign < one_half ) return                                    


   pressure_gradient(1) = pressure_curv_gradient(1) * csi(1,i,j,k) + &
                          pressure_curv_gradient(2) * eta(1,i,j,k) + & 
                          pressure_curv_gradient(3) * zet(1,i,j,k)

   pressure_gradient(2) = pressure_curv_gradient(1) * csi(2,i,j,k) + &
                          pressure_curv_gradient(2) * eta(2,i,j,k) + & 
                          pressure_curv_gradient(3) * zet(2,i,j,k)

   pressure_gradient(3) = pressure_curv_gradient(1) * csi(3,i,j,k) + &
                          pressure_curv_gradient(2) * eta(3,i,j,k) + & 
                          pressure_curv_gradient(3) * zet(3,i,j,k)



end subroutine calc_pressure_gradient