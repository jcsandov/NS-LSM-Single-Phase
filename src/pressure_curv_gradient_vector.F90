subroutine pressure_curv_gradient_vector(i, j, k, pressure_curv_gradient, exsign)


   use AdvectionMethods

   ! compute derivatives of cartesian velocities in curvilinear directions
   ! It uses adaptative stencils depending on whether the node is within the 
   ! water phase or not
   ! velocity_curv_gradient(i,l) = ∂u_i/∂ξ^l

   implicit none

   integer, intent(in) :: i,j,k
   real (kind = rdf), dimension(1:3), intent(out) :: pressure_curv_gradient 
   
   integer :: ll , mm 
   real (kind = rdf), intent(inout) :: exsign ! flag variable. exsign = 0, the gradient is set to zero
                               !                exsign = 1, the computed gradient is returned

   ! local variables

   real ( kind = rdf ) :: rLL , rL , rC , rR , rRR
   real ( kind = rdf ) :: pLL , pL , pC , pR , pRR

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
      rC  = rsign ( phi( i   , j , k ) )  
      rR  = rsign ( phi( i+1 , j , k ) )  
      rRR = rsign ( phi( i+2 , j , k ) )  

      pLL = zero                  
      pL  = zero                  
      pC  = q ( 1, i   , j , k )  
      pR  = q ( 1, i+1 , j , k )  
      pRR = q ( 1, i+2 , j , k )  
      
   else if ( i == ista + 1 ) then
         
      rLL      = zero ! dummy
      rL       = rsign ( phi( i-1 , j , k ) )
      rC       = rsign ( phi( i   , j , k ) )
      rR       = rsign ( phi( i+1 , j , k ) )
      rRR      = rsign ( phi( i+2 , j , k ) )

      pLL = zero                 
      pL  = q ( 1, i-1 , j , k ) 
      pC  = q ( 1, i   , j , k ) 
      pR  = q ( 1, i+1 , j , k ) 
      pRR = q ( 1, i+2 , j , k ) 

   else if ( i == iend - 1 ) then
         
      rLL      = rsign ( phi( i-2 , j , k ) )
      rL       = rsign ( phi( i-1 , j , k ) )
      rC       = rsign ( phi( i   , j , k ) )
      rR       = rsign ( phi( i+1 , j , k ) )
      rRR      = zero ! dummy

      pLL = q ( 1, i-2 , j , k ) 
      pL  = q ( 1, i-1 , j , k ) 
      pC  = q ( 1, i   , j , k ) 
      pR  = q ( 1, i+1 , j , k ) 
      pRR = zero                 


   else if ( i == iend .or. BlankingFlagFront ) then
         
      bias_csi = -1

      rLL      = rsign ( phi( i-2 , j , k ) )
      rL       = rsign ( phi( i-1 , j , k ) )
      rC       = rsign ( phi( i   , j , k ) )
      rR       = zero ! dummy
      rRR      = zero ! dummy

      pLL = q ( 1, i-2 , j , k ) 
      pL  = q ( 1, i-1 , j , k ) 
      pC  = q ( 1, i   , j , k ) 
      pR  = zero                 
      pRR = zero                 

   else

      rLL      = rsign ( phi( i-2 , j , k ) )
      rL       = rsign ( phi( i-1 , j , k ) )
      rC       = rsign ( phi( i   , j , k ) )
      rR       = rsign ( phi( i+1 , j , k ) )
      rRR      = rsign ( phi( i+2 , j , k ) )

      pLL = q ( 1, i-2 , j , k ) 
      pL  = q ( 1, i-1 , j , k ) 
      pC  = q ( 1, i   , j , k ) 
      pR  = q ( 1, i+1 , j , k ) 
      pRR = q ( 1, i+2 , j , k ) 

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


   if ( j == jsta .or. BlankingFlagRight ) then

      bias_eta = 1

      rLL = zero                    
      rL  = zero                   
      rC  = rsign ( phi (i , j   , k ) )  
      rR  = rsign ( phi (i , j+1 , k ) )  
      rRR = rsign ( phi (i , j+2 , k ) )  

      pLL = zero                 
      pL  = zero                 
      pC  = q ( 1, i , j   , k ) 
      pR  = q ( 1, i , j+1 , k ) 
      pRR = q ( 1, i , j+2 , k ) 
   
   else if ( j == jsta + 1 ) then
         
      rLL      = zero ! dummy
      rL       = rsign ( phi( i , j-1 , k ) )
      rC       = rsign ( phi( i , j   , k ) )
      rR       = rsign ( phi( i , j+1 , k ) )
      rRR      = rsign ( phi( i , j+2 , k ) )

      pLL = zero                 
      pL  = q ( 1, i , j-1 , k ) 
      pC  = q ( 1, i , j   , k ) 
      pR  = q ( 1, i , j+1 , k ) 
      pRR = q ( 1, i , j+2 , k ) 

   else if ( j == jend - 1 ) then
         
      rLL      = rsign ( phi( i , j-2 , k ) )
      rL       = rsign ( phi( i , j-1 , k ) )
      rC       = rsign ( phi( i , j   , k ) )
      rR       = rsign ( phi( i , j+1 , k ) )
      rRR      = zero ! dummy

      pLL = q ( 1, i , j-2 , k ) 
      pL  = q ( 1, i , j-1 , k ) 
      pC  = q ( 1, i , j   , k ) 
      pR  = q ( 1, i , j+1 , k ) 
      pRR = zero                 


   else if ( j == jend .or. BlankingFlagLeft ) then
         
      bias_eta = -1

      rLL      = rsign ( phi( i , j-2 , k ) )
      rL       = rsign ( phi( i , j-1 , k ) )
      rC       = rsign ( phi( i , j   , k ) )
      rR       = zero ! dummy
      rRR      = zero ! dummy

      pLL = q ( 1, i , j-2 , k ) 
      pL  = q ( 1, i , j-1 , k ) 
      pC  = q ( 1, i , j   , k ) 
      pR  = zero                 
      pRR = zero                 

   else

      rLL      = rsign ( phi( i , j-2 , k ) )
      rL       = rsign ( phi( i , j-1 , k ) )
      rC       = rsign ( phi( i , j   , k ) )
      rR       = rsign ( phi( i , j+1 , k ) )
      rRR      = rsign ( phi( i , j+2 , k ) )

      pLL = q ( 1, i , j-2 , k ) 
      pL  = q ( 1, i , j-1 , k ) 
      pC  = q ( 1, i , j   , k ) 
      pR  = q ( 1, i , j+1 , k ) 
      pRR = q ( 1, i , j+2 , k ) 

   end if

   exsign = zero

   ! 

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
      rC  = rsign ( phi( i , j , k   ) )  
      rR  = rsign ( phi( i , j , k+1 ) )  
      rRR = rsign ( phi( i , j , k+2 ) )  

      pLL = zero                 
      pL  = zero                 
      pC  = q ( 1, i , j , k   ) 
      pR  = q ( 1, i , j , k+1 ) 
      pRR = q ( 1, i , j , k+2 ) 
      
   else if ( k == ksta + 1 ) then
         
      rLL      = zero ! dummy
      rL       = rsign ( phi( i , j , k-1 ) )
      rC       = rsign ( phi( i , j , k   ) )
      rR       = rsign ( phi( i , j , k+1 ) )
      rRR      = rsign ( phi( i , j , k+2 ) )

      pLL = zero                  
      pL  = q ( 1, i , j , k-1 )  
      pC  = q ( 1, i , j , k   )  
      pR  = q ( 1, i , j , k+1 )  
      pRR = q ( 1, i , j , k+2 )  

   else if ( k == kend - 1 ) then
         
      rLL      = rsign ( phi( i , j , k-2 ) )
      rL       = rsign ( phi( i , j , k-1 ) )
      rC       = rsign ( phi( i , j , k   ) )
      rR       = rsign ( phi( i , j , k+1 ) )
      rRR      = zero ! dummy

      pLL = q ( 1, i , j , k-2 ) 
      pL  = q ( 1, i , j , k-1 ) 
      pC  = q ( 1, i , j , k   ) 
      pR  = q ( 1, i , j , k+1 ) 
      pRR = zero                 


   else if ( k == kend ) then
         
      bias_zet = -1

      rLL      = rsign ( phi( i , j , k-2 ) )
      rL       = rsign ( phi( i , j , k-1 ) )
      rC       = rsign ( phi( i , j , k   ) )
      rR       = zero ! dummy
      rRR      = zero ! dummy

      pLL = q ( 1, i , j , k-2 ) 
      pL  = q ( 1, i , j , k-1 ) 
      pC  = q ( 1, i , j , k   ) 
      pR  = zero                 
      pRR = zero                 

   else

      rLL      = rsign ( phi( i , j , k-2 ) )
      rL       = rsign ( phi( i , j , k-1 ) )
      rC       = rsign ( phi( i , j , k   ) )
      rR       = rsign ( phi( i , j , k+1 ) )
      rRR      = rsign ( phi( i , j , k+2 ) )

      pLL = q ( 1, i , j , k-2 ) 
      pL  = q ( 1, i , j , k-1 ) 
      pC  = q ( 1, i , j , k   ) 
      pR  = q ( 1, i , j , k+1 ) 
      pRR = q ( 1, i , j , k+2 ) 

   end if

   exsign = zero

   ! 

   ! ∂p/∂ζ
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     pLL, pL , pC, pR, pRR        , &
                                     dz                           , &
                                     bias_zet                     , &
                                     exsign                       , &
                                     pressure_curv_gradient(3)      &
                                    )

   if ( exsign < one_half ) return                                    
                                    

end subroutine pressure_curv_gradient_vector