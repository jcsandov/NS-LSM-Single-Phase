subroutine velocity_curv_gradient_tensor1(i, j, k, velocity_curv_gradient, exsign)
   
   use AdvectionMethods

   ! compute derivatives of cartesian velocities in curvilinear directions
   ! It uses adaptative stencils depending on whether the node is within the 
   ! water phase or not
   ! velocity_curv_gradient(i,l) = ∂u_i/∂ξ^l

   implicit none

   integer, intent(in) :: i,j,k
   real (kind = rdf), dimension(1:3,1:3), intent(inout) :: velocity_curv_gradient 
   
   integer :: ll , mm 
   real (kind = rdf), intent(inout) :: exsign ! flag variable. exsign = 0, the gradient is set to zero
                               !                exsign = 1, the computed gradient is returned

   ! local variables

   real ( kind = rdf ) :: rLL , rL , rC , rR , rRR
   real ( kind = rdf ) :: uLL , uL , uC , uR , uRR
   real ( kind = rdf ) :: vLL , vL , vC , vR , vRR
   real ( kind = rdf ) :: wLL , wL , wC , wR , wRR

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

      uLL = zero                  ;  vLL = zero                  ;  wLL = zero                 ;
      uL  = zero                  ;  vL  = zero                  ;  wL  = zero                 ;
      uC  = q ( 2, i   , j , k )  ;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
      uR  = q ( 2, i+1 , j , k )  ;  vR  = q ( 3, i+1 , j , k )  ;  wR  = q ( 4, i+1 , j , k ) ;
      uRR = q ( 2, i+2 , j , k )  ;  vRR = q ( 3, i+2 , j , k )  ;  wRR = q ( 4, i+2 , j , k ) ;
      
   else if ( i == ista + 1 ) then
         
      rLL      = zero ! dummy
      rL       = rsign ( i-1 , j , k )
      rC       = rsign ( i   , j , k )
      rR       = rsign ( i+1 , j , k )
      rRR      = rsign ( i+2 , j , k )

      uLL = zero                  ;  vLL = zero                  ;  wLL = zero                 ;
      uL  = q ( 2, i-1 , j , k )  ;  vL  = q ( 3, i-1 , j , k )  ;  wL  = q ( 4, i-1 , j , k ) ;
      uC  = q ( 2, i   , j , k )  ;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
      uR  = q ( 2, i+1 , j , k )  ;  vR  = q ( 3, i+1 , j , k )  ;  wR  = q ( 4, i+1 , j , k ) ;
      uRR = q ( 2, i+2 , j , k )  ;  vRR = q ( 3, i+2 , j , k )  ;  wRR = q ( 4, i+2 , j , k ) ;

   else if ( i == iend - 1 ) then
         
      rLL      = rsign ( i-2 , j , k )
      rL       = rsign ( i-1 , j , k )
      rC       = rsign ( i   , j , k )
      rR       = rsign ( i+1 , j , k )
      rRR      = zero ! dummy

      uLL = q ( 2, i-2 , j , k )  ;  vLL = q ( 3, i-2 , j , k )  ;  wLL = q ( 4, i-2 , j , k ) ;
      uL  = q ( 2, i-1 , j , k )  ;  vL  = q ( 3, i-1 , j , k )  ;  wL  = q ( 4, i-1 , j , k ) ;
      uC  = q ( 2, i   , j , k )  ;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
      uR  = q ( 2, i+1 , j , k )  ;  vR  = q ( 3, i+1 , j , k )  ;  wR  = q ( 4, i+1 , j , k ) ;
      uRR = zero                  ;  vRR = zero                  ;  wRR = zero                 ;


   else if ( i == iend ) then
         
      bias_csi = -1

      rLL      = rsign ( i-2 , j , k )
      rL       = rsign ( i-1 , j , k )
      rC       = rsign ( i   , j , k )
      rR       = zero ! dummy
      rRR      = zero ! dummy

      uLL = q ( 2, i-2 , j , k )  ;  vLL = q ( 3, i-2 , j , k )  ;  wLL = q ( 4, i-2 , j , k ) ;
      uL  = q ( 2, i-1 , j , k )  ;  vL  = q ( 3, i-1 , j , k )  ;  wL  = q ( 4, i-1 , j , k ) ;
      uC  = q ( 2, i   , j , k )  ;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
      uR  = zero                  ;  vR  = zero                  ;  wR  = zero                 ;
      uRR = zero                  ;  vRR = zero                  ;  wRR = zero                 ;

   else

      rLL      = rsign ( i-2 , j , k )
      rL       = rsign ( i-1 , j , k )
      rC       = rsign ( i   , j , k )
      rR       = rsign ( i+1 , j , k )
      rRR      = rsign ( i+2 , j , k )

      uLL = q ( 2, i-2 , j , k )  ;  vLL = q ( 3, i-2 , j , k )  ;  wLL = q ( 4, i-2 , j , k ) ;
      uL  = q ( 2, i-1 , j , k )  ;  vL  = q ( 3, i-1 , j , k )  ;  wL  = q ( 4, i-1 , j , k ) ;
      uC  = q ( 2, i   , j , k )  ;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
      uR  = q ( 2, i+1 , j , k )  ;  vR  = q ( 3, i+1 , j , k )  ;  wR  = q ( 4, i+1 , j , k ) ;
      uRR = q ( 2, i+2 , j , k )  ;  vRR = q ( 3, i+2 , j , k )  ;  wRR = q ( 4, i+2 , j , k ) ;

   end if

   exsign = zero

   ! ∂u/∂ξ
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     uLL, uL , uC, uR, uRR        , &
                                     dc                           , &
                                     bias_csi                     , &
                                     exsign                       , &
                                     velocity_curv_gradient(1,1)    &
                                    )

   if ( exsign < one_half ) return                                    
                                    

   ! ∂v/∂ξ
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     vLL, vL , vC, vR, vRR        , &
                                     dc                           , &
                                     bias_csi                     , &
                                     exsign                       , &
                                     velocity_curv_gradient(2,1)    &
                                    )

   if ( exsign < one_half ) return                                    


   ! ∂w/∂ξ
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     wLL, wL , wC, wR, wRR        , &
                                     dc                           , &
                                     bias_csi                     , &
                                     exsign                       , &
                                     velocity_curv_gradient(3,1)    &
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

      uLL = zero                  ;  vLL = zero                  ;  wLL = zero                 ;
      uL  = zero                  ;  vL  = zero                  ;  wL  = zero                 ;
      uC  = q ( 2, i , j   , k )  ;  vC  = q ( 3, i , j   , k )  ;  wC  = q ( 4, i , j   , k ) ;
      uR  = q ( 2, i , j+1 , k )  ;  vR  = q ( 3, i , j+1 , k )  ;  wR  = q ( 4, i , j+1 , k ) ;
      uRR = q ( 2, i , j+2 , k )  ;  vRR = q ( 3, i , j+2 , k )  ;  wRR = q ( 4, i , j+2 , k ) ;
   
   else if ( j == jsta + 1 ) then
         
      rLL      = zero ! dummy
      rL       = rsign ( i , j-1 , k )
      rC       = rsign ( i , j   , k )
      rR       = rsign ( i , j+1 , k )
      rRR      = rsign ( i , j+2 , k )

      uLL = zero                  ;  vLL = zero                  ;  wLL = zero                 ;
      uL  = q ( 2, i , j-1 , k )  ;  vL  = q ( 3, i , j-1 , k )  ;  wL  = q ( 4, i , j-1 , k ) ;
      uC  = q ( 2, i , j   , k )  ;  vC  = q ( 3, i , j   , k )  ;  wC  = q ( 4, i , j   , k ) ;
      uR  = q ( 2, i , j+1 , k )  ;  vR  = q ( 3, i , j+1 , k )  ;  wR  = q ( 4, i , j+1 , k ) ;
      uRR = q ( 2, i , j+2 , k )  ;  vRR = q ( 3, i , j+2 , k )  ;  wRR = q ( 4, i , j+2 , k ) ;

   else if ( j == jend - 1 ) then
         
      rLL      = rsign ( i , j-2 , k )
      rL       = rsign ( i , j-1 , k )
      rC       = rsign ( i , j   , k )
      rR       = rsign ( i , j+1 , k )
      rRR      = zero ! dummy

      uLL = q ( 2, i , j-2 , k )  ;  vLL = q ( 3, i , j-2 , k )  ;  wLL = q ( 4, i , j-2 , k ) ;
      uL  = q ( 2, i , j-1 , k )  ;  vL  = q ( 3, i , j-1 , k )  ;  wL  = q ( 4, i , j-1 , k ) ;
      uC  = q ( 2, i , j   , k )  ;  vC  = q ( 3, i , j   , k )  ;  wC  = q ( 4, i , j   , k ) ;
      uR  = q ( 2, i , j+1 , k )  ;  vR  = q ( 3, i , j+1 , k )  ;  wR  = q ( 4, i , j+1 , k ) ;
      uRR = zero                  ;  vRR = zero                  ;  wRR = zero                 ;


   else if ( j == jend ) then
         
      bias_eta = -1

      rLL      = rsign ( i , j-2 , k )
      rL       = rsign ( i , j-1 , k )
      rC       = rsign ( i , j   , k )
      rR       = zero ! dummy
      rRR      = zero ! dummy

      uLL = q ( 2, i , j-2 , k )  ;  vLL = q ( 3, i , j-2 , k )  ;  wLL = q ( 4, i , j-2 , k ) ;
      uL  = q ( 2, i , j-1 , k )  ;  vL  = q ( 3, i , j-1 , k )  ;  wL  = q ( 4, i , j-1 , k ) ;
      uC  = q ( 2, i , j   , k )  ;  vC  = q ( 3, i , j   , k )  ;  wC  = q ( 4, i , j   , k ) ;
      uR  = zero                  ;  vR  = zero                  ;  wR  = zero                 ;
      uRR = zero                  ;  vRR = zero                  ;  wRR = zero                 ;

   else

      rLL      = rsign ( i , j-2 , k )
      rL       = rsign ( i , j-1 , k )
      rC       = rsign ( i , j   , k )
      rR       = rsign ( i , j+1 , k )
      rRR      = rsign ( i , j+2 , k )

      uLL = q ( 2, i , j-2 , k )  ;  vLL = q ( 3, i , j-2 , k )  ;  wLL = q ( 4, i , j-2 , k ) ;
      uL  = q ( 2, i , j-1 , k )  ;  vL  = q ( 3, i , j-1 , k )  ;  wL  = q ( 4, i , j-1 , k ) ;
      uC  = q ( 2, i , j   , k )  ;  vC  = q ( 3, i , j   , k )  ;  wC  = q ( 4, i , j   , k ) ;
      uR  = q ( 2, i , j+1 , k )  ;  vR  = q ( 3, i , j+1 , k )  ;  wR  = q ( 4, i , j+1 , k ) ;
      uRR = q ( 2, i , j+2 , k )  ;  vRR = q ( 3, i , j+2 , k )  ;  wRR = q ( 4, i , j+2 , k ) ;

   end if

   exsign = zero

   ! 

   ! ∂u/∂η
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     uLL, uL , uC, uR, uRR        , &
                                     de                           , &
                                     bias_eta                     , &
                                     exsign                       , &
                                     velocity_curv_gradient(1,2)    &
                                    )

   if ( exsign < one_half ) return                                    
                                    

   ! ∂v/∂η
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     vLL, vL , vC, vR, vRR        , &
                                     de                           , &
                                     bias_eta                     , &
                                     exsign                       , &
                                     velocity_curv_gradient(2,2)    &
                                    )

   if ( exsign < one_half ) return                                    


   ! ∂w/∂η
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     wLL, wL , wC, wR, wRR        , &
                                     de                           , &
                                     bias_eta                     , &
                                     exsign                       , &
                                     velocity_curv_gradient(3,2)    &
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

      uLL = zero                  ;  vLL = zero                  ;  wLL = zero                 ;
      uL  = zero                  ;  vL  = zero                  ;  wL  = zero                 ;
      uC  = q ( 2, i , j , k   )  ;  vC  = q ( 3, i , j , k   )  ;  wC  = q ( 4, i , j , k   ) ;
      uR  = q ( 2, i , j , k+1 )  ;  vR  = q ( 3, i , j , k+1 )  ;  wR  = q ( 4, i , j , k+1 ) ;
      uRR = q ( 2, i , j , k+2 )  ;  vRR = q ( 3, i , j , k+2 )  ;  wRR = q ( 4, i , j , k+2 ) ;
      
   else if ( k == ksta + 1 ) then
         
      rLL      = zero ! dummy
      rL       = rsign ( i , j , k-1 )
      rC       = rsign ( i , j , k   )
      rR       = rsign ( i , j , k+1 )
      rRR      = rsign ( i , j , k+2 )

      uLL = zero                  ;  vLL = zero                  ;  wLL = zero                 ;
      uL  = q ( 2, i , j , k-1 )  ;  vL  = q ( 3, i , j , k-1 )  ;  wL  = q ( 4, i , j , k-1 ) ;
      uC  = q ( 2, i , j , k   )  ;  vC  = q ( 3, i , j , k   )  ;  wC  = q ( 4, i , j , k   ) ;
      uR  = q ( 2, i , j , k+1 )  ;  vR  = q ( 3, i , j , k+1 )  ;  wR  = q ( 4, i , j , k+1 ) ;
      uRR = q ( 2, i , j , k+2 )  ;  vRR = q ( 3, i , j , k+2 )  ;  wRR = q ( 4, i , j , k+2 ) ;

   else if ( k == kend - 1 ) then
         
      rLL      = rsign ( i , j , k-2 )
      rL       = rsign ( i , j , k-1 )
      rC       = rsign ( i , j , k   )
      rR       = rsign ( i , j , k+1 )
      rRR      = zero ! dummy

      uLL = q ( 2, i , j , k-2 )  ;  vLL = q ( 3, i , j , k-2 )  ;  wLL = q ( 4, i , j , k-2 ) ;
      uL  = q ( 2, i , j , k-1 )  ;  vL  = q ( 3, i , j , k-1 )  ;  wL  = q ( 4, i , j , k-1 ) ;
      uC  = q ( 2, i , j , k   )  ;  vC  = q ( 3, i , j , k   )  ;  wC  = q ( 4, i , j , k   ) ;
      uR  = q ( 2, i , j , k+1 )  ;  vR  = q ( 3, i , j , k+1 )  ;  wR  = q ( 4, i , j , k+1 ) ;
      uRR = zero                  ;  vRR = zero                  ;  wRR = zero                 ;


   else if ( k == kend ) then
         
      bias_zet = -1

      rLL      = rsign ( i , j , k-2 )
      rL       = rsign ( i , j , k-1 )
      rC       = rsign ( i , j , k   )
      rR       = zero ! dummy
      rRR      = zero ! dummy

      uLL = q ( 2, i , j , k-2 )  ;  vLL = q ( 3, i , j , k-2 )  ;  wLL = q ( 4, i , j , k-2 ) ;
      uL  = q ( 2, i , j , k-1 )  ;  vL  = q ( 3, i , j , k-1 )  ;  wL  = q ( 4, i , j , k-1 ) ;
      uC  = q ( 2, i , j , k   )  ;  vC  = q ( 3, i , j , k   )  ;  wC  = q ( 4, i , j , k   ) ;
      uR  = zero                  ;  vR  = zero                  ;  wR  = zero                 ;
      uRR = zero                  ;  vRR = zero                  ;  wRR = zero                 ;

   else

      rLL      = rsign ( i , j , k-2 )
      rL       = rsign ( i , j , k-1 )
      rC       = rsign ( i , j , k   )
      rR       = rsign ( i , j , k+1 )
      rRR      = rsign ( i , j , k+2 )

      uLL = q ( 2, i , j , k-2 )  ;  vLL = q ( 3, i , j , k-2 )  ;  wLL = q ( 4, i , j , k-2 ) ;
      uL  = q ( 2, i , j , k-1 )  ;  vL  = q ( 3, i , j , k-1 )  ;  wL  = q ( 4, i , j , k-1 ) ;
      uC  = q ( 2, i , j , k   )  ;  vC  = q ( 3, i , j , k   )  ;  wC  = q ( 4, i , j , k   ) ;
      uR  = q ( 2, i , j , k+1 )  ;  vR  = q ( 3, i , j , k+1 )  ;  wR  = q ( 4, i , j , k+1 ) ;
      uRR = q ( 2, i , j , k+2 )  ;  vRR = q ( 3, i , j , k+2 )  ;  wRR = q ( 4, i , j , k+2 ) ;

   end if

   exsign = zero

   ! 

   ! ∂u/∂ζ
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     uLL, uL , uC, uR, uRR        , &
                                     dz                           , &
                                     bias_zet                     , &
                                     exsign                       , &
                                     velocity_curv_gradient(1,3)    &
                                    )

   if ( exsign < one_half ) return                                    
                                    

   ! ∂v/∂ζ
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     vLL, vL , vC, vR, vRR        , &
                                     dz                           , &
                                     bias_zet                     , &
                                     exsign                       , &
                                     velocity_curv_gradient(2,3)    &
                                    )

   if ( exsign < one_half ) return                                    


   ! ∂w/∂ζ
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     wLL, wL , wC, wR, wRR        , &
                                     dz                           , &
                                     bias_zet                     , &
                                     exsign                       , &
                                     velocity_curv_gradient(3,3)    &
                                    )

   if ( exsign < one_half ) return                                    


end subroutine velocity_curv_gradient_tensor1


