
subroutine calc_pressure_gradient_across_interface( i, j, k, pressure_gradient , exsign )
   
   use AdvectionMethods

   ! compute derivatives of cartesian velocities in curvilinear directions
   ! It uses adaptative stencils depending on whether the node is within the 
   ! water phase or not

   implicit none

   integer, intent(in) :: i,j,k
   real (kind = rdf), dimension(1:3) , intent(out)   :: pressure_gradient 
   real (kind = rdf), intent(inout) :: exsign ! flag variable. exsign = 0, the gradient is set to zero
                               !                exsign = 1, the computed gradient is returned

   ! local variables

   real ( kind = rdf ) :: phiLL , phiL , phiC , phiR , phiRR
   real ( kind = rdf ) :: pLL   , pL   , pC   , pR   , pRR
   real ( kind = rdf ) , dimension(3) :: pressure_curv_gradient 

   integer :: bias_csi , bias_eta , bias_zet

   ! Set the bias of the derivative if I'm at one of the boundaries
   bias_csi = 0
   bias_eta = 0
   bias_zet = 0

   pressure_gradient = zero

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !                                       ξ - direction
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if ( i == ista ) then

      bias_csi = 1

      phiLL = zero                 ;  pLL   = zero ! dummy                      
      phiL  = zero                 ;  pL    = zero ! dummy                    
      phiC  = phi ( i   , j , k )  ;  pC    = q ( 1, i   , j , k ) 
      phiR  = phi ( i+1 , j , k )  ;  pR    = q ( 1, i+1 , j , k ) 
      phiRR = phi ( i+2 , j , k )  ;  pRR   = q ( 1, i+2 , j , k ) 
      
   else if ( i == ista + 1 ) then
   
      phiLL = zero                 ;  pLL   = zero ! dummy                       
      phiL  = phi ( i-1 , j , k )  ;  pL    = q ( 1, i-1 , j , k )
      phiC  = phi ( i   , j , k )  ;  pC    = q ( 1, i   , j , k )
      phiR  = phi ( i+1 , j , k )  ;  pR    = q ( 1, i+1 , j , k )
      phiRR = phi ( i+2 , j , k )  ;  pRR   = q ( 1, i+2 , j , k )

   else if ( i == iend - 1 ) then

      phiLL = phi ( i-2 , j , k )  ;  pLL   = q ( 1, i-2 , j , k )
      phiL  = phi ( i-1 , j , k )  ;  pL    = q ( 1, i-1 , j , k )
      phiC  = phi ( i   , j , k )  ;  pC    = q ( 1, i   , j , k )
      phiR  = phi ( i+1 , j , k )  ;  pR    = q ( 1, i+1 , j , k )
      phiRR = zero                 ;  pRR   = zero                             

   else if ( i == iend ) then
         
      bias_csi = -1

      phiLL  = phi ( i-2 , j , k ) ;  pLL    = q ( 1, i-2 , j , k )
      phiL   = phi ( i-1 , j , k ) ;  pL     = q ( 1, i-1 , j , k )
      phiC   = phi ( i   , j , k ) ;  pC     = q ( 1, i   , j , k )
      phiR   = zero                ;  pR     = zero                 
      phiRR  = zero                ;  pRR    = zero                 

   else

      phiLL = phi ( i-2 , j , k )  ;  pLL   = q ( 1 , i-2 , j , k ) 
      phiL  = phi ( i-1 , j , k )  ;  pL    = q ( 1 , i-1 , j , k ) 
      phiC  = phi ( i   , j , k )  ;  pC    = q ( 1 , i   , j , k ) 
      phiR  = phi ( i+1 , j , k )  ;  pR    = q ( 1 , i+1 , j , k ) 
      phiRR = phi ( i+2 , j , k )  ;  pRR   = q ( 1 , i+2 , j , k ) 

   end if

   exsign = zero

   ! ∂p/∂ξ
   call NearSurface_P_FirstDerivative ( phiLL, phiL , phiC , phiR , phiRR , &
                                        pLL  , pL   , pC   , pR   , pRR   , &
                                        dc                                , &
                                        bias_csi                          , &
                                        exsign                            , &
                                        pressure_curv_gradient(1)           &
                                      )

   if ( exsign < one_half ) return                                    


   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !                                       η - direction
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


   if ( j == jsta ) then

      bias_eta = 1

      phiLL = zero                 ;  pLL   = zero ! dummy                      
      phiL  = zero                 ;  pL    = zero ! dummy                    
      phiC  = phi ( i , j   , k )  ;  pC    = q ( 1, i , j   , k ) 
      phiR  = phi ( i , j+1 , k )  ;  pR    = q ( 1, i , j+1 , k ) 
      phiRR = phi ( i , j+2 , k )  ;  pRR   = q ( 1, i , j+2 , k ) 
      
   else if ( j == jsta + 1 ) then
         
      phiLL = zero                 ;  pLL   = zero ! dummy                       
      phiL  = phi ( i , j-1 , k )  ;  pL    = q ( 1, i , j-1 , k )
      phiC  = phi ( i , j   , k )  ;  pC    = q ( 1, i , j   , k )
      phiR  = phi ( i , j+1 , k )  ;  pR    = q ( 1, i , j+1 , k )
      phiRR = phi ( i , j+2 , k )  ;  pRR   = q ( 1, i , j+2 , k )

   else if ( j == jend - 1 ) then

      phiLL = phi ( i , j-2 , k )  ;  pLL   = q ( 1, i , j-2 , k )
      phiL  = phi ( i , j-1 , k )  ;  pL    = q ( 1, i , j-1 , k )
      phiC  = phi ( i , j   , k )  ;  pC    = q ( 1, i , j   , k )
      phiR  = phi ( i , j+1 , k )  ;  pR    = q ( 1, i , j+1 , k )
      phiRR = zero                 ;  pRR   = zero                             

   else if ( j == jend ) then

      bias_eta = -1

      phiLL  = phi ( i , j-2 , k ) ;  pLL    = q ( 1, i , j-2 , k )
      phiL   = phi ( i , j-1 , k ) ;  pL     = q ( 1, i , j-1 , k )
      phiC   = phi ( i , j   , k ) ;  pC     = q ( 1, i , j   , k )
      phiR   = zero                ;  pR     = zero                 
      phiRR  = zero                ;  pRR    = zero                 

   else

      phiLL = phi ( i , j-2 , k )  ;  pLL   = q ( 1 , i , j-2 , k ) 
      phiL  = phi ( i , j-1 , k )  ;  pL    = q ( 1 , i , j-1 , k ) 
      phiC  = phi ( i , j   , k )  ;  pC    = q ( 1 , i , j   , k ) 
      phiR  = phi ( i , j+1 , k )  ;  pR    = q ( 1 , i , j+1 , k ) 
      phiRR = phi ( i , j+2 , k )  ;  pRR   = q ( 1 , i , j+2 , k ) 

   end if

   exsign = zero

   ! ∂p/∂η
   call NearSurface_P_FirstDerivative ( phiLL, phiL , phiC , phiR , phiRR , &
                                        pLL  , pL   , pC   , pR   , pRR   , &
                                        de                                , &
                                        bias_eta                          , &
                                        exsign                            , &
                                        pressure_curv_gradient(2)           &
                                      )

   if ( exsign < one_half ) return                                    
                              
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !                                       ζ - direction
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if ( k == ksta ) then

      bias_zet = 1

      phiLL = zero                 ;  pLL   = zero ! dummy                      
      phiL  = zero                 ;  pL    = zero ! dummy                    
      phiC  = phi ( i , j , k   )  ;  pC    = q ( 1, i , j , k   ) 
      phiR  = phi ( i , j , k+1 )  ;  pR    = q ( 1, i , j , k+1 ) 
      phiRR = phi ( i , j , k+2 )  ;  pRR   = q ( 1, i , j , k+2 ) 
      
   else if ( k == ksta + 1 ) then
      
      phiLL = zero                 ;  pLL   = zero ! dummy                       
      phiL  = phi ( i , j , k-1 )  ;  pL    = q ( 1, i , j , k-1 )
      phiC  = phi ( i , j , k   )  ;  pC    = q ( 1, i , j , k   )
      phiR  = phi ( i , j , k+1 )  ;  pR    = q ( 1, i , j , k+1 )
      phiRR = phi ( i , j , k+2 )  ;  pRR   = q ( 1, i , j , k+2 )

   else if ( k == kend - 1 ) then

      phiLL = phi ( i , j , k-2 )  ;  pLL   = q ( 1, i , j , k-2 )
      phiL  = phi ( i , j , k-1 )  ;  pL    = q ( 1, i , j , k-1 )
      phiC  = phi ( i , j , k   )  ;  pC    = q ( 1, i , j , k   )
      phiR  = phi ( i , j , k+1 )  ;  pR    = q ( 1, i , j , k+1 )
      phiRR = zero                 ;  pRR   = zero                             

   else if ( k == kend ) then
         
      bias_zet = -1

      phiLL  = phi ( i , j , k-2 ) ;  pLL    = q ( 1, i , j , k-2 )
      phiL   = phi ( i , j , k-1 ) ;  pL     = q ( 1, i , j , k-1 )
      phiC   = phi ( i , j , k   ) ;  pC     = q ( 1, i , j , k   )
      phiR   = zero                ;  pR     = zero                 
      phiRR  = zero                ;  pRR    = zero                 

   else

      phiLL = phi ( i , j , k-2 )  ;  pLL   = q ( 1 , i , j , k-2 ) 
      phiL  = phi ( i , j , k-1 )  ;  pL    = q ( 1 , i , j , k-1 ) 
      phiC  = phi ( i , j , k   )  ;  pC    = q ( 1 , i , j , k   ) 
      phiR  = phi ( i , j , k+1 )  ;  pR    = q ( 1 , i , j , k+1 ) 
      phiRR = phi ( i , j , k+2 )  ;  pRR   = q ( 1 , i , j , k+2 ) 

   end if

   exsign = zero

   ! ∂p/∂ζ
   call NearSurface_P_FirstDerivative ( phiLL, phiL , phiC , phiR , phiRR , &
                                        pLL  , pL   , pC   , pR   , pRR   , &
                                        dz                                , &
                                        bias_zet                          , &
                                        exsign                            , &
                                        pressure_curv_gradient(3)           &
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



end subroutine calc_pressure_gradient_across_interface