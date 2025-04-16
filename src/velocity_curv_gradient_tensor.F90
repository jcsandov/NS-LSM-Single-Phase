subroutine velocity_curv_gradient_tensor(i, j, k, velocity_curv_gradient, exsign)


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


   logical :: BlankingFlagFront , BlankingFlagBack , BlankingFlagLeft , BlankingFlagRight


   ! Set the bias of the derivative if I'm at one of the boundaries
   bias_csi = 0
   bias_eta = 0
   bias_zet = 0

   BlankingFlagFront  = .false.
   BlankingFlagBack   = .false.

   BlankingFlagLeft   = .false.
   BlankingFlagRight  = .false.

   if ( nblk /= 0 ) then
      
      do nb = 1,nblk 

         if ( i == li_blk_ia(1,nb) .and. le_blk_ia(1,nb) > ista .and. &
              j >= li_blk_ja(1,nb) .and. j <= li_blk_jb(1,nb) ) then
            
            BlankingFlagFront = .true.

         end if   

         if ( i == li_blk_ib(1,nb) .and. le_blk_ib(1,nb) < iend .and. &
              j >= li_blk_ja(1,nb) .and. j <= li_blk_jb(1,nb) ) then
            
            BlankingFlagBack = .true.

         end if   

         if ( j == li_blk_ja(1,nb) .and. le_blk_ja(1,nb) > jsta .and. &
              i >= li_blk_ia(1,nb) .and. i <= li_blk_ib(1,nb) ) then
            
            BlankingFlagLeft = .true.

         end if   

         if ( j == li_blk_jb(1,nb) .and. le_blk_jb(1,nb) < jend .and. &
              i >= li_blk_ia(1,nb) .and. i <= li_blk_ib(1,nb) ) then
            
            BlankingFlagRight = .true.

         end if   

      end do

   end if

   if ( nblke /= 0 ) then
      
      do nb = 1,nblke 

         if ( i == le_blk_ia(1,nb) .and. le_blk_ia(1,nb) > ista .and. &
              j >= le_blk_ja(1,nb) .and. j <= le_blk_jb(1,nb) ) then
            
            BlankingFlagFront = .true.

         end if   

         if ( i == le_blk_ib(1,nb) .and. le_blk_ib(1,nb) < iend .and. &
              j >= le_blk_ja(1,nb) .and. j <= le_blk_jb(1,nb) ) then
            
            BlankingFlagBack = .true.

         end if   

         if ( j == le_blk_ja(1,nb) .and. le_blk_ja(1,nb) > jsta .and. &
              i >= le_blk_ia(1,nb) .and. i <= le_blk_ib(1,nb) ) then
            
            BlankingFlagLeft = .true.

         end if   

         if ( j == le_blk_jb(1,nb) .and. le_blk_jb(1,nb) < jend .and. &
              i >= le_blk_ia(1,nb) .and. i <= le_blk_ib(1,nb) ) then
            
            BlankingFlagRight = .true.

         end if   

      end do

   end if

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !                                       ξ - direction
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if ( i == ista .or. BlankingFlagBack ) then

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


   else if ( i == iend .or. BlankingFlagFront ) then
         
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


   if ( j == jsta .or. BlankingFlagRight ) then

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


   else if ( j == jend .or. BlankingFlagLeft ) then
         
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


! The version below was an upwind version of the derivatives. I 
! implemented it to check if that was the problem when I was imposing
! p = 0 at the free-surface. It was too restrictive and was generating
! problems near the blanking region
   
!   use AdvectionMethods
!
!   ! compute derivatives of cartesian velocities in curvilinear directions
!   ! It uses adaptative stencils depending on whether the node is within the 
!   ! water phase or not
!   ! velocity_curv_gradient(i,l) = ∂u_i/∂ξ^l
!   ! 
!   ! If possible, it uses and upwind derivative for the velocity
!   
!   implicit none
!
!   integer, intent(in) :: i,j,k
!   real (kind = rdf), dimension(1:3,1:3), intent(inout) :: velocity_curv_gradient 
!   
!   integer :: ll , mm 
!   real (kind = rdf), intent(inout) :: exsign ! flag variable. exsign = 0, the gradient is set to zero
!                               !                exsign = 1, the computed gradient is returned
!
!   ! local variables
!
!   real ( kind = rdf ) :: rLL , rL , rC , rR , rRR
!   real ( kind = rdf ) :: uLL , uL , uC , uR , uRR
!   real ( kind = rdf ) :: vLL , vL , vC , vR , vRR
!   real ( kind = rdf ) :: wLL , wL , wC , wR , wRR
!
!   real ( kind = rdf ) :: ucon_1 , ucon_2 , ucon_3
!   real ( kind = rdf ) :: up1 , um1 , up2 , um2 , up3 , um3
!
!   logical :: BlankingFlagFront , BlankingFlagBack , BlankingFlagLeft , BlankingFlagRight
!
!
!   integer :: bias_csi , bias_eta , bias_zet
!
!   ! Set the bias of the derivative if I'm at one of the boundaries
!   bias_csi = 0
!   bias_eta = 0
!   bias_zet = 0
!
!   !  Contravariant velocity
!   ucon_1 = csi(1,i,j,k) * q(2,i,j,k) + &
!            csi(2,i,j,k) * q(3,i,j,k) + &
!            csi(3,i,j,k) * q(4,i,j,k)    
!
!   ucon_2 = eta(1,i,j,k) * q(2,i,j,k) + &
!            eta(2,i,j,k) * q(3,i,j,k) + &
!            eta(3,i,j,k) * q(4,i,j,k)    
!
!   ucon_3 = zet(1,i,j,k) * q(2,i,j,k) + &
!            zet(2,i,j,k) * q(3,i,j,k) + &
!            zet(3,i,j,k) * q(4,i,j,k)    
!
!
!   up1 = zero
!   um1 = zero
!
!   up2 = zero
!   um2 = zero
!
!   up3 = zero
!   um3 = zero
!
!   if ( ucon_1 > 0 ) then
!      up1 = one
!   else
!      um1 = one
!   end if
!
!   if ( ucon_2 > 0 ) then
!      up2 = one
!   else
!      um2 = one
!   end if
!
!   if ( ucon_3 > 0 ) then
!      up3 = one
!   else
!      um3 = one
!   end if
!
!   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!   !                                       ξ - direction
!   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!    
!    BlankingFlagFront  = .false.
!    BlankingFlagBack   = .false.
!
!    BlankingFlagLeft   = .false.
!    BlankingFlagRight  = .false.
!
!   if ( nblk /= 0 ) then
!      
!      do nb = 1,nblk 
!
!         if ( i == li_blk_ia(1,nb) .and. &
!              j >= li_blk_ja(1,nb) .and. j <= li_blk_jb(1,nb) ) then
!            
!            BlankingFlagFront = .true.
!
!         end if   
!
!         if ( i == li_blk_ib(1,nb) .and. &
!              j >= li_blk_ja(1,nb) .and. j <= li_blk_jb(1,nb) ) then
!            
!            BlankingFlagBack = .true.
!
!         end if   
!
!         if ( j == li_blk_ja(1,nb) .and. &
!              i >= li_blk_ia(1,nb) .and. i <= li_blk_ib(1,nb) ) then
!            
!            BlankingFlagLeft = .true.
!
!         end if   
!
!         if ( j == li_blk_jb(1,nb) .and. &
!              i >= li_blk_ia(1,nb) .and. i <= li_blk_ib(1,nb) ) then
!            
!            BlankingFlagRight = .true.
!
!         end if   
!
!      end do
!
!   end if
!      
!   if ( i == ista .or. BlankingFlagBack ) then
!
!      bias_csi = 1
!
!      rLL = zero                    
!      rL  = zero                   
!      rC  = rsign ( i   , j , k )  
!      rR  = rsign ( i+1 , j , k )  
!      rRR = rsign ( i+2 , j , k )  
!
!      uLL = zero                  ;  vLL = zero                  ;  wLL = zero                 ;
!      uL  = zero                  ;  vL  = zero                  ;  wL  = zero                 ;
!      uC  = q ( 2, i   , j , k )  ;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
!      uR  = q ( 2, i+1 , j , k )  ;  vR  = q ( 3, i+1 , j , k )  ;  wR  = q ( 4, i+1 , j , k ) ;
!      uRR = q ( 2, i+2 , j , k )  ;  vRR = q ( 3, i+2 , j , k )  ;  wRR = q ( 4, i+2 , j , k ) ;
! 
!      ! ∂u/∂ξ
!      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
!                                        uLL, uL , uC, uR, uRR        , &
!                                        dc                           , &
!                                        bias_csi                     , &
!                                        exsign                       , &
!                                        velocity_curv_gradient(1,1)    &
!                                       )
!   
!      if ( exsign < one_half .and. i==15 .and. j==23 .and. k==26 ) print *, 'l153'                                    
!      if ( exsign < one_half ) return                                    
!                                       
!   
!      ! ∂v/∂ξ
!      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
!                                        vLL, vL , vC, vR, vRR        , &
!                                        dc                           , &
!                                        bias_csi                     , &
!                                        exsign                       , &
!                                        velocity_curv_gradient(2,1)    &
!                                       )
!   
!      if ( exsign < one_half ) print *, 'l166'                                    
!      if ( exsign < one_half ) return                                    
!   
!   
!      ! ∂w/∂ξ
!      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
!                                        wLL, wL , wC, wR, wRR        , &
!                                        dc                           , &
!                                        bias_csi                     , &
!                                        exsign                       , &
!                                        velocity_curv_gradient(3,1)    &
!                                       )
!   
!      if ( exsign < one_half .and. i==15 .and. j==23 .and. k==26 ) print *, 'l179'                                    
!      if ( exsign < one_half ) return                                    
!
!   else if ( i == iend .or. BlankingFlagFront ) then
!
!      bias_csi = -1
!
!      rLL      = rsign ( i-2 , j , k )
!      rL       = rsign ( i-1 , j , k )
!      rC       = rsign ( i   , j , k )
!      rR       = zero ! dummy
!      rRR      = zero ! dummy
!
!      uLL = q ( 2, i-2 , j , k )  ;  vLL = q ( 3, i-2 , j , k )  ;  wLL = q ( 4, i-2 , j , k ) ;
!      uL  = q ( 2, i-1 , j , k )  ;  vL  = q ( 3, i-1 , j , k )  ;  wL  = q ( 4, i-1 , j , k ) ;
!      uC  = q ( 2, i   , j , k )  ;  vC  = q ( 3, i   , j , k )  ;  wC  = q ( 4, i   , j , k ) ;
!      uR  = zero                  ;  vR  = zero                  ;  wR  = zero                 ;
!      uRR = zero                  ;  vRR = zero                  ;  wRR = zero                 ;
!
!      ! ∂u/∂ξ
!      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
!                                        uLL, uL , uC, uR, uRR        , &
!                                        dc                           , &
!                                        bias_csi                     , &
!                                        exsign                       , &
!                                        velocity_curv_gradient(1,1)    &
!                                       )
!   
!      if ( exsign < one_half .and. i==15 .and. j==23 .and. k==26 ) print *, 'l207'                                    
!      if ( exsign < one_half ) return                                    
!                                       
!   
!      ! ∂v/∂ξ
!      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
!                                        vLL, vL , vC, vR, vRR        , &
!                                        dc                           , &
!                                        bias_csi                     , &
!                                        exsign                       , &
!                                        velocity_curv_gradient(2,1)    &
!                                       )
!   
!      if ( exsign < one_half .and. i==15 .and. j==23 .and. k==26 ) print *, 'l220'                                    
!      if ( exsign < one_half ) return                                    
!   
!   
!      ! ∂w/∂ξ
!      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
!                                        wLL, wL , wC, wR, wRR        , &
!                                        dc                           , &
!                                        bias_csi                     , &
!                                        exsign                       , &
!                                        velocity_curv_gradient(3,1)    &
!                                       )
!   
!      if ( exsign < one_half .and. i==15 .and. j==23 .and. k==26 ) print *, 'l233'                                    
!      if ( exsign < one_half ) return     
!
!   else
!
!      rL = ( rsign ( i-1 , j , k ) + abs(rsign ( i-1 , j , k )) ) / two
!      rC = ( rsign ( i   , j , k ) + abs(rsign ( i   , j , k )) ) / two
!      rR = ( rsign ( i+1 , j , k ) + abs(rsign ( i+1 , j , k )) ) / two
!
!      exsign = up1 * ( rC * rL ) + &
!               um1 * ( rR * rC )
!      
!      if ( exsign < one_half .and. i==15 .and. j==23 .and. k==26 ) print *, 'l245'                                    
!      if ( exsign < one_half ) return
!
!      velocity_curv_gradient(1,1) = up1 * dc * ( q(2,i,j,k)   - q(2,i-1,j,k)  ) + &
!                                    um1 * dc * ( q(2,i+1,j,k) - q(2,i,j,k)  )
!
!      velocity_curv_gradient(2,1) = up1 * dc * ( q(3,i,j,k)   - q(3,i-1,j,k)  ) + &
!                                    um1 * dc * ( q(3,i+1,j,k) - q(3,i,j,k)  )
!
!      velocity_curv_gradient(3,1) = up1 * dc * ( q(4,i,j,k)   - q(4,i-1,j,k)  ) + &
!                                    um1 * dc * ( q(4,i+1,j,k) - q(4,i,j,k)  )
!   end if
!
!
!   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!   !                                       η - direction
!   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!   if ( j == jsta .or. BlankingFlagRight ) then
!
!      bias_eta = 1
!
!      rLL = zero                    
!      rL  = zero                   
!      rC  = rsign ( i, j    , k )  
!      rR  = rsign ( i, j+1  , k )  
!      rRR = rsign ( i, j+2  , k )  
!
!      uLL = zero                  ;  vLL = zero                  ;  wLL = zero                 ;
!      uL  = zero                  ;  vL  = zero                  ;  wL  = zero                 ;
!      uC  = q ( 2, i , j   , k )  ;  vC  = q ( 3, i , j   , k )  ;  wC  = q ( 4, i , j   , k ) ;
!      uR  = q ( 2, i , j+1 , k )  ;  vR  = q ( 3, i , j+1 , k )  ;  wR  = q ( 4, i , j+1 , k ) ;
!      uRR = q ( 2, i , j+2 , k )  ;  vRR = q ( 3, i , j+2 , k )  ;  wRR = q ( 4, i , j+2 , k ) ;
!
!      ! ∂u/∂η
!      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
!                                        uLL, uL , uC, uR, uRR        , &
!                                        de                           , &
!                                        bias_eta                     , &
!                                        exsign                       , &
!                                        velocity_curv_gradient(1,2)    &
!                                       )
!   
!      if ( exsign < one_half .and. i==15 .and. j==23 .and. k==26 ) print *, 'l288'                                    
!      if ( exsign < one_half ) return                                    
!                                       
!   
!      ! ∂v/∂η
!      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
!                                        vLL, vL , vC, vR, vRR        , &
!                                        de                           , &
!                                        bias_eta                     , &
!                                        exsign                       , &
!                                        velocity_curv_gradient(2,2)    &
!                                       )
!   
!      if ( exsign < one_half .and. i==15 .and. j==23 .and. k==26 ) print *, 'l301'                                    
!      if ( exsign < one_half ) return                                    
!   
!   
!      ! ∂w/∂η
!      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
!                                        wLL, wL , wC, wR, wRR        , &
!                                        dc                           , &
!                                        bias_eta                     , &
!                                        exsign                       , &
!                                        velocity_curv_gradient(3,2)    &
!                                       )
!   
!      if ( exsign < one_half .and. i==15 .and. j==23 .and. k==26 ) print *, 'l314'                                    
!      if ( exsign < one_half ) return      
!   
!   else if ( j == jend .or. BlankingFlagLeft ) then
!
!      bias_eta = -1
!
!      rLL      = rsign ( i , j-2 , k )
!      rL       = rsign ( i , j-1 , k )
!      rC       = rsign ( i , j   , k )
!      rR       = zero ! dummy
!      rRR      = zero ! dummy
!
!      uLL = q ( 2, i , j-2 , k )  ;  vLL = q ( 3, i , j-2 , k )  ;  wLL = q ( 4, i , j-2 , k ) ;
!      uL  = q ( 2, i , j-1 , k )  ;  vL  = q ( 3, i , j-1 , k )  ;  wL  = q ( 4, i , j-1 , k ) ;
!      uC  = q ( 2, i , j   , k )  ;  vC  = q ( 3, i , j   , k )  ;  wC  = q ( 4, i , j   , k ) ;
!      uR  = zero                  ;  vR  = zero                  ;  wR  = zero                 ;
!      uRR = zero                  ;  vRR = zero                  ;  wRR = zero                 ;
!
!      ! ∂u/∂η
!      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
!                                        uLL, uL , uC, uR, uRR        , &
!                                        de                           , &
!                                        bias_eta                     , &
!                                        exsign                       , &
!                                        velocity_curv_gradient(1,2)    &
!                                       )
!   
!      if ( exsign < one_half .and. i==15 .and. j==23 .and. k==26 ) print *, 'l342'                                    
!      if ( exsign < one_half ) return                                    
!                                       
!   
!      ! ∂v/∂η
!      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
!                                        vLL, vL , vC, vR, vRR        , &
!                                        de                           , &
!                                        bias_eta                     , &
!                                        exsign                       , &
!                                        velocity_curv_gradient(2,2)    &
!                                       )
!   
!      if ( exsign < one_half .and. i==15 .and. j==23 .and. k==26 ) print *, 'l355'                                    
!      if ( exsign < one_half ) return                                    
!   
!   
!      ! ∂w/∂η
!      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
!                                        wLL, wL , wC, wR, wRR        , &
!                                        dc                           , &
!                                        bias_eta                     , &
!                                        exsign                       , &
!                                        velocity_curv_gradient(3,2)    &
!                                       )
!   
!      if ( exsign < one_half .and. i==15 .and. j==23 .and. k==26 ) print *, 'l368'                                    
!      if ( exsign < one_half ) return   
!
!   else
!
!      rL = ( rsign ( i , j-1 , k ) + abs(rsign ( i , j-1 , k )) ) / two
!      rC = ( rsign ( i , j   , k ) + abs(rsign ( i , j   , k )) ) / two
!      rR = ( rsign ( i , j+1 , k ) + abs(rsign ( i , j+1 , k )) ) / two
!
!      exsign = up2 * ( rC * rL ) + um2 * ( rC * rR )
!      
!      if ( exsign < one_half .and. i==15 .and. j==23 .and. k==26 ) print *, 'l379'                                    
!      if ( exsign < one_half ) return
!
!      velocity_curv_gradient(1,2) = up2 * de * ( q(2,i,j,k)   - q(2,i,j-1,k)  ) + &
!                                    um2 * de * ( q(2,i,j+1,k) - q(2,i,j,k)  )
!
!      velocity_curv_gradient(2,2) = up2 * de * ( q(3,i,j,k)   - q(3,i,j-1,k)  ) + &
!                                    um2 * de * ( q(3,i,j+1,k) - q(3,i,j,k)  )
!
!      velocity_curv_gradient(3,2) = up2 * de * ( q(4,i,j,k)   - q(4,i,j-1,k)  ) + &
!                                    um2 * de * ( q(4,i,j+1,k) - q(4,i,j,k)  )
!   end if
!
!
!   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!   !                                       ζ - direction
!   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!   
!   if ( k == ksta ) then
!
!      bias_zet = 1
!
!      rLL = zero                    
!      rL  = zero                   
!      rC  = rsign ( i , j , k   )  
!      rR  = rsign ( i , j , k+1 )  
!      rRR = rsign ( i , j , k+2 )  
!
!      uLL = zero                  ;  vLL = zero                  ;  wLL = zero                 ;
!      uL  = zero                  ;  vL  = zero                  ;  wL  = zero                 ;
!      uC  = q ( 2, i , j , k   )  ;  vC  = q ( 3, i , j , k   )  ;  wC  = q ( 4, i , j , k   ) ;
!      uR  = q ( 2, i , j , k+1 )  ;  vR  = q ( 3, i , j , k+1 )  ;  wR  = q ( 4, i , j , k+1 ) ;
!      uRR = q ( 2, i , j , k+2 )  ;  vRR = q ( 3, i , j , k+2 )  ;  wRR = q ( 4, i , j , k+2 ) ;
! 
!      ! ∂u/∂ζ
!      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
!                                        uLL, uL , uC, uR, uRR        , &
!                                        dz                           , &
!                                        bias_zet                     , &
!                                        exsign                       , &
!                                        velocity_curv_gradient(1,3)    &
!                                       )
!   
!      if ( exsign < one_half .and. i==15 .and. j==23 .and. k==26 ) print *, 'l422'                                    
!      if ( exsign < one_half ) return                                    
!                                       
!   
!      ! ∂v/∂ζ
!      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
!                                        vLL, vL , vC, vR, vRR        , &
!                                        dz                           , &
!                                        bias_zet                     , &
!                                        exsign                       , &
!                                        velocity_curv_gradient(2,3)    &
!                                       )
!   
!      if ( exsign < one_half .and. i==15 .and. j==23 .and. k==26 ) print *, 'l435'                                    
!      if ( exsign < one_half ) return                                    
!   
!   
!      ! ∂w/∂ζ
!      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
!                                        wLL, wL , wC, wR, wRR        , &
!                                        dz                           , &
!                                        bias_zet                     , &
!                                        exsign                       , &
!                                        velocity_curv_gradient(3,3)    &
!                                       )
!   
!      if ( exsign < one_half .and. i==15 .and. j==23 .and. k==26 ) print *, 'l448'                                    
!      if ( exsign < one_half ) return   
!
!   else if ( k == kend ) then
!
!      bias_zet = -1
!
!      rLL      = rsign ( i , j , k-2 )
!      rL       = rsign ( i , j , k-1 )
!      rC       = rsign ( i , j , k   )
!      rR       = zero ! dummy
!      rRR      = zero ! dummy
!
!      uLL = q ( 2, i , j , k-2 )  ;  vLL = q ( 3, i , j , k-2 )  ;  wLL = q ( 4, i , j , k-2 ) ;
!      uL  = q ( 2, i , j , k-1 )  ;  vL  = q ( 3, i , j , k-1 )  ;  wL  = q ( 4, i , j , k-1 ) ;
!      uC  = q ( 2, i , j , k   )  ;  vC  = q ( 3, i , j , k   )  ;  wC  = q ( 4, i , j , k   ) ;
!      uR  = zero                  ;  vR  = zero                  ;  wR  = zero                 ;
!      uRR = zero                  ;  vRR = zero                  ;  wRR = zero                 ;
!
!      ! ∂u/∂ζ
!      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
!                                        uLL, uL , uC, uR, uRR        , &
!                                        dz                           , &
!                                        bias_zet                     , &
!                                        exsign                       , &
!                                        velocity_curv_gradient(1,3)    &
!                                       )
!   
!      if ( exsign < one_half .and. i==15 .and. j==23 .and. k==26 ) print *, 'l476'                                    
!      if ( exsign < one_half ) return                                    
!                                       
!   
!      ! ∂v/∂ζ
!      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
!                                        vLL, vL , vC, vR, vRR        , &
!                                        dz                           , &
!                                        bias_zet                     , &
!                                        exsign                       , &
!                                        velocity_curv_gradient(2,3)    &
!                                       )
!   
!      if ( exsign < one_half .and. i==15 .and. j==23 .and. k==26 ) print *, 'l489'                                    
!      if ( exsign < one_half ) return                                    
!   
!   
!      ! ∂w/∂ζ
!      call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
!                                        wLL, wL , wC, wR, wRR        , &
!                                        dz                           , &
!                                        bias_zet                     , &
!                                        exsign                       , &
!                                        velocity_curv_gradient(3,3)    &
!                                       )
!   
!      if ( exsign < one_half .and. i==15 .and. j==23 .and. k==26 ) print *, 'l502'                                    
!      if ( exsign < one_half ) return  
!
!   else
!
!      rL = ( rsign ( i , j , k-1 ) + abs(rsign ( i , j , k-1 )) ) / two
!      rC = ( rsign ( i , j , k   ) + abs(rsign ( i , j , k   )) ) / two
!      rR = ( rsign ( i , j , k+1 ) + abs(rsign ( i , j , k+1 )) ) / two
!
!      exsign = up3 * ( rC * rL ) + um3 * ( rC * rR )
!      
!      if ( exsign < one_half .and. i==15 .and. j==23 .and. k==26 ) print *, 'l513'                                    
!      if ( exsign < one_half ) return
!
!      velocity_curv_gradient(1,3) = up3 * dz * ( q(2,i,j,k)   - q(2,i,j,k-1)  ) + &
!                                    um3 * dz * ( q(2,i,j,k+1) - q(2,i,j,k)  )
!
!      velocity_curv_gradient(2,3) = up3 * dz * ( q(3,i,j,k)   - q(3,i,j,k-1)  ) + &
!                                    um3 * dz * ( q(3,i,j,k+1) - q(3,i,j,k)  )
!
!      velocity_curv_gradient(3,3) = up3 * dz * ( q(4,i,j,k)   - q(4,i,j,k-1)  ) + &
!                                    um3 * dz * ( q(4,i,j,k+1) - q(4,i,j,k)  )
!   end if


end subroutine velocity_curv_gradient_tensor