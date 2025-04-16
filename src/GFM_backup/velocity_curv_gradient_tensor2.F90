subroutine velocity_curv_gradient_tensor2(i, j, k, velocity_gradient, exsign)
   
   use AdvectionMethods

   ! compute derivatives of cartesian velocities in curvilinear directions
   ! It uses adaptative stencils depending on whether the node is within the 
   ! water phase or not
   ! velocity_curv_gradient(i,l) = ∂u_i/∂ξ^l

   implicit none

   integer, intent(in) :: i,j,k
   real (kind = rdf), dimension(1:3,1:3), intent(inout) :: velocity_gradient 
   
   integer :: ll , mm 
   real (kind = rdf), intent(inout) :: exsign ! flag variable. exsign = 0, the gradient is set to zero
                               !                exsign = 1, the computed gradient is returned

   ! local variables

   real ( kind = rdf ) :: rLL , rL , rC , rR , rRR
   real ( kind = rdf ) :: u1LL , u1L , u1C , u1R , u1RR
   real ( kind = rdf ) :: u2LL , u2L , u2C , u2R , u2RR
   real ( kind = rdf ) :: u3LL , u3L , u3C , u3R , u3RR
   real ( kind = rdf ) :: u4LL , u4L , u4C , u4R , u4RR
   real ( kind = rdf ) :: u5LL , u5L , u5C , u5R , u5RR
   real ( kind = rdf ) :: u6LL , u6L , u6C , u6R , u6RR
   real ( kind = rdf ) :: u7LL , u7L , u7C , u7R , u7RR
   real ( kind = rdf ) :: u8LL , u8L , u8C , u8R , u8RR
   real ( kind = rdf ) :: u9LL , u9L , u9C , u9R , u9RR

   real ( kind = rdf ) :: v1LL , v1L , v1C , v1R , v1RR
   real ( kind = rdf ) :: v2LL , v2L , v2C , v2R , v2RR
   real ( kind = rdf ) :: v3LL , v3L , v3C , v3R , v3RR
   real ( kind = rdf ) :: v4LL , v4L , v4C , v4R , v4RR
   real ( kind = rdf ) :: v5LL , v5L , v5C , v5R , v5RR
   real ( kind = rdf ) :: v6LL , v6L , v6C , v6R , v6RR
   real ( kind = rdf ) :: v7LL , v7L , v7C , v7R , v7RR
   real ( kind = rdf ) :: v8LL , v8L , v8C , v8R , v8RR
   real ( kind = rdf ) :: v9LL , v9L , v9C , v9R , v9RR

   real ( kind = rdf ) :: w1LL , w1L , w1C , w1R , w1RR
   real ( kind = rdf ) :: w2LL , w2L , w2C , w2R , w2RR
   real ( kind = rdf ) :: w3LL , w3L , w3C , w3R , w3RR
   real ( kind = rdf ) :: w4LL , w4L , w4C , w4R , w4RR
   real ( kind = rdf ) :: w5LL , w5L , w5C , w5R , w5RR
   real ( kind = rdf ) :: w6LL , w6L , w6C , w6R , w6RR
   real ( kind = rdf ) :: w7LL , w7L , w7C , w7R , w7RR
   real ( kind = rdf ) :: w8LL , w8L , w8C , w8R , w8RR
   real ( kind = rdf ) :: w9LL , w9L , w9C , w9R , w9RR

   real ( kind = rdf ) :: du_dx , du1_dcsi , du4_deta , du7_dzet  
   real ( kind = rdf ) :: du_dy , du2_dcsi , du5_deta , du8_dzet  
   real ( kind = rdf ) :: du_dz , du3_dcsi , du6_deta , du9_dzet  

   real ( kind = rdf ) :: dv_dx , dv1_dcsi , dv4_deta , dv7_dzet  
   real ( kind = rdf ) :: dv_dy , dv2_dcsi , dv5_deta , dv8_dzet  
   real ( kind = rdf ) :: dv_dz , dv3_dcsi , dv6_deta , dv9_dzet  

   real ( kind = rdf ) :: dw_dx , dw1_dcsi , dw4_deta , dw7_dzet  
   real ( kind = rdf ) :: dw_dy , dw2_dcsi , dw5_deta , dw8_dzet  
   real ( kind = rdf ) :: dw_dz , dw3_dcsi , dw6_deta , dw9_dzet  


   integer :: bias_csi , bias_eta , bias_zet

   ! Set the bias of the derivative if I'm at one of the boundaries
   bias_csi = 0
   bias_eta = 0
   bias_zet = 0

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   !                                       ξ - direction
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   velocity_gradient = zero
   
   if ( i == ista ) then

      bias_csi = 1

      rLL = zero                    
      rL  = zero                   
      rC  = rsign ( i   , j , k )  
      rR  = rsign ( i+1 , j , k )  
      rRR = rsign ( i+2 , j , k )  

      u1LL = zero                                    
      u1L  = zero                                    
      u1C  = q(2,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      u1R  = q(2,i+1,j,k) * csi(1,i+1,j,k) / aj(i+1,j,k) 
      u1RR = q(2,i+2,j,k) * csi(1,i+2,j,k) / aj(i+2,j,k) 

      u2LL = zero                                    
      u2L  = zero                                    
      u2C  = q(2,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      u2R  = q(2,i+1,j,k) * csi(2,i+1,j,k) / aj(i+1,j,k) 
      u2RR = q(2,i+2,j,k) * csi(2,i+2,j,k) / aj(i+2,j,k) 

      u3LL = zero                                    
      u3L  = zero                                    
      u3C  = q(2,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      u3R  = q(2,i+1,j,k) * csi(3,i+1,j,k) / aj(i+1,j,k) 
      u3RR = q(2,i+2,j,k) * csi(3,i+2,j,k) / aj(i+2,j,k) 

      !-------------------------------------------------------

      v1LL = zero                                    
      v1L  = zero                                    
      v1C  = q(3,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      v1R  = q(3,i+1,j,k) * csi(1,i+1,j,k) / aj(i+1,j,k) 
      v1RR = q(3,i+2,j,k) * csi(1,i+2,j,k) / aj(i+2,j,k) 

      v2LL = zero                                    
      v2L  = zero                                    
      v2C  = q(3,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      v2R  = q(3,i+1,j,k) * csi(2,i+1,j,k) / aj(i+1,j,k) 
      v2RR = q(3,i+2,j,k) * csi(2,i+2,j,k) / aj(i+2,j,k) 

      v3LL = zero                                    
      v3L  = zero                                    
      v3C  = q(3,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      v3R  = q(3,i+1,j,k) * csi(3,i+1,j,k) / aj(i+1,j,k) 
      v3RR = q(3,i+2,j,k) * csi(3,i+2,j,k) / aj(i+2,j,k) 

      !-------------------------------------------------------
   
      w1LL = zero                                    
      w1L  = zero                                    
      w1C  = q(4,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      w1R  = q(4,i+1,j,k) * csi(1,i+1,j,k) / aj(i+1,j,k) 
      w1RR = q(4,i+2,j,k) * csi(1,i+2,j,k) / aj(i+2,j,k) 

      w2LL = zero                                    
      w2L  = zero                                    
      w2C  = q(4,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      w2R  = q(4,i+1,j,k) * csi(2,i+1,j,k) / aj(i+1,j,k) 
      w2RR = q(4,i+2,j,k) * csi(2,i+2,j,k) / aj(i+2,j,k) 

      w3LL = zero                                    
      w3L  = zero                                    
      w3C  = q(4,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      w3R  = q(4,i+1,j,k) * csi(3,i+1,j,k) / aj(i+1,j,k) 
      w3RR = q(4,i+2,j,k) * csi(3,i+2,j,k) / aj(i+2,j,k) 


   else if ( i == ista + 1 ) then
         
      rLL      = zero ! dummy
      rL       = rsign ( i-1 , j , k )
      rC       = rsign ( i   , j , k )
      rR       = rsign ( i+1 , j , k )
      rRR      = rsign ( i+2 , j , k )

      u1LL = zero                                    
      u1L  = q(2,i-1,j,k) * csi(1,i-1,j,k) / aj(i-1,j,k)                                    
      u1C  = q(2,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      u1R  = q(2,i+1,j,k) * csi(1,i+1,j,k) / aj(i+1,j,k) 
      u1RR = q(2,i+2,j,k) * csi(1,i+2,j,k) / aj(i+2,j,k) 

      u2LL = zero                                    
      u2L  = q(2,i-1,j,k) * csi(2,i-1,j,k) / aj(i-1,j,k) 
      u2C  = q(2,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      u2R  = q(2,i+1,j,k) * csi(2,i+1,j,k) / aj(i+1,j,k) 
      u2RR = q(2,i+2,j,k) * csi(2,i+2,j,k) / aj(i+2,j,k) 

      u3LL = zero                                    
      u3L  = q(2,i-1,j,k) * csi(3,i-1,j,k) / aj(i-1,j,k) 
      u3C  = q(2,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      u3R  = q(2,i+1,j,k) * csi(3,i+1,j,k) / aj(i+1,j,k) 
      u3RR = q(2,i+2,j,k) * csi(3,i+2,j,k) / aj(i+2,j,k) 

      !-------------------------------------------------------

      v1LL = zero                                    
      v1L  = q(3,i-1,j,k) * csi(1,i-1,j,k) / aj(i-1,j,k) 
      v1C  = q(3,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      v1R  = q(3,i+1,j,k) * csi(1,i+1,j,k) / aj(i+1,j,k) 
      v1RR = q(3,i+2,j,k) * csi(1,i+2,j,k) / aj(i+2,j,k) 

      v2LL = zero                                    
      v2L  = q(3,i-1,j,k) * csi(2,i-1,j,k) / aj(i-1,j,k) 
      v2C  = q(3,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      v2R  = q(3,i+1,j,k) * csi(2,i+1,j,k) / aj(i+1,j,k) 
      v2RR = q(3,i+2,j,k) * csi(2,i+2,j,k) / aj(i+2,j,k) 

      v3LL = zero                                    
      v3L  = q(3,i-1,j,k) * csi(3,i-1,j,k) / aj(i-1,j,k) 
      v3C  = q(3,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      v3R  = q(3,i+1,j,k) * csi(3,i+1,j,k) / aj(i+1,j,k) 
      v3RR = q(3,i+2,j,k) * csi(3,i+2,j,k) / aj(i+2,j,k) 

      !-------------------------------------------------------
   
      w1LL = zero                                    
      w1L  = q(4,i-1,j,k) * csi(1,i-1,j,k) / aj(i-1,j,k) 
      w1C  = q(4,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      w1R  = q(4,i+1,j,k) * csi(1,i+1,j,k) / aj(i+1,j,k) 
      w1RR = q(4,i+2,j,k) * csi(1,i+2,j,k) / aj(i+2,j,k) 

      w2LL = zero                                    
      w2L  = q(4,i-1,j,k) * csi(2,i-1,j,k) / aj(i-1,j,k) 
      w2C  = q(4,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      w2R  = q(4,i+1,j,k) * csi(2,i+1,j,k) / aj(i+1,j,k) 
      w2RR = q(4,i+2,j,k) * csi(2,i+2,j,k) / aj(i+2,j,k) 

      w3LL = zero                                    
      w3L  = q(4,i-1,j,k) * csi(3,i-1,j,k) / aj(i-1,j,k) 
      w3C  = q(4,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      w3R  = q(4,i+1,j,k) * csi(3,i+1,j,k) / aj(i+1,j,k) 
      w3RR = q(4,i+2,j,k) * csi(3,i+2,j,k) / aj(i+2,j,k) 
   
   else if ( i == iend - 1 ) then
         
      rLL      = rsign ( i-2 , j , k )
      rL       = rsign ( i-1 , j , k )
      rC       = rsign ( i   , j , k )
      rR       = rsign ( i+1 , j , k )
      rRR      = zero ! dummy

      u1LL = q(2,i-2,j,k) * csi(1,i-2,j,k) / aj(i-2,j,k)                                    
      u1L  = q(2,i-1,j,k) * csi(1,i-1,j,k) / aj(i-1,j,k)                                    
      u1C  = q(2,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      u1R  = q(2,i+1,j,k) * csi(1,i+1,j,k) / aj(i+1,j,k) 
      u1RR = zero

      u2LL = q(2,i-2,j,k) * csi(2,i-2,j,k) / aj(i-2,j,k) 
      u2L  = q(2,i-1,j,k) * csi(2,i-1,j,k) / aj(i-1,j,k) 
      u2C  = q(2,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      u2R  = q(2,i+1,j,k) * csi(2,i+1,j,k) / aj(i+1,j,k) 
      u2RR = zero 

      u3LL = q(2,i-2,j,k) * csi(3,i-2,j,k) / aj(i-2,j,k) 
      u3L  = q(2,i-1,j,k) * csi(3,i-1,j,k) / aj(i-1,j,k) 
      u3C  = q(2,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      u3R  = q(2,i+1,j,k) * csi(3,i+1,j,k) / aj(i+1,j,k) 
      u3RR = zero 

      !-------------------------------------------------------

      v1LL = q(3,i-2,j,k) * csi(1,i-2,j,k) / aj(i-2,j,k) 
      v1L  = q(3,i-1,j,k) * csi(1,i-1,j,k) / aj(i-1,j,k) 
      v1C  = q(3,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      v1R  = q(3,i+1,j,k) * csi(1,i+1,j,k) / aj(i+1,j,k) 
      v1RR = zero 

      v2LL = q(3,i-2,j,k) * csi(2,i-2,j,k) / aj(i-2,j,k) 
      v2L  = q(3,i-1,j,k) * csi(2,i-1,j,k) / aj(i-1,j,k) 
      v2C  = q(3,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      v2R  = q(3,i+1,j,k) * csi(2,i+1,j,k) / aj(i+1,j,k) 
      v2RR = zero

      v3LL = q(3,i-2,j,k) * csi(3,i-2,j,k) / aj(i-2,j,k) 
      v3L  = q(3,i-1,j,k) * csi(3,i-1,j,k) / aj(i-1,j,k) 
      v3C  = q(3,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      v3R  = q(3,i+1,j,k) * csi(3,i+1,j,k) / aj(i+1,j,k) 
      v3RR = zero 

      !-------------------------------------------------------
   
      w1LL = q(4,i-2,j,k) * csi(1,i-2,j,k) / aj(i-2,j,k) 
      w1L  = q(4,i-1,j,k) * csi(1,i-1,j,k) / aj(i-1,j,k) 
      w1C  = q(4,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      w1R  = q(4,i+1,j,k) * csi(1,i+1,j,k) / aj(i+1,j,k) 
      w1RR = zero 

      w2LL = q(4,i-2,j,k) * csi(2,i-2,j,k) / aj(i-2,j,k) 
      w2L  = q(4,i-1,j,k) * csi(2,i-1,j,k) / aj(i-1,j,k) 
      w2C  = q(4,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      w2R  = q(4,i+1,j,k) * csi(2,i+1,j,k) / aj(i+1,j,k) 
      w2RR = zero 

      w3LL = q(4,i-2,j,k) * csi(3,i-2,j,k) / aj(i-2,j,k) 
      w3L  = q(4,i-1,j,k) * csi(3,i-1,j,k) / aj(i-1,j,k) 
      w3C  = q(4,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      w3R  = q(4,i+1,j,k) * csi(3,i+1,j,k) / aj(i+1,j,k) 
      w3RR = zero 


   else if ( i == iend ) then
         
      bias_csi = -1

      rLL      = rsign ( i-2 , j , k )
      rL       = rsign ( i-1 , j , k )
      rC       = rsign ( i   , j , k )
      rR       = zero ! dummy
      rRR      = zero ! dummy

      u1LL = q(2,i-2,j,k) * csi(1,i-2,j,k) / aj(i-2,j,k)                                    
      u1L  = q(2,i-1,j,k) * csi(1,i-1,j,k) / aj(i-1,j,k)                                    
      u1C  = q(2,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      u1R  = zero 
      u1RR = zero 

      u2LL = q(2,i-2,j,k) * csi(2,i-2,j,k) / aj(i-2,j,k) 
      u2L  = q(2,i-1,j,k) * csi(2,i-1,j,k) / aj(i-1,j,k) 
      u2C  = q(2,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      u2R  = zero 
      u2RR = zero 

      u3LL = q(2,i-2,j,k) * csi(3,i-2,j,k) / aj(i-2,j,k) 
      u3L  = q(2,i-1,j,k) * csi(3,i-1,j,k) / aj(i-1,j,k) 
      u3C  = q(2,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      u3R  = zero 
      u3RR = zero 

      !-------------------------------------------------------

      v1LL = q(3,i-2,j,k) * csi(1,i-2,j,k) / aj(i-2,j,k) 
      v1L  = q(3,i-1,j,k) * csi(1,i-1,j,k) / aj(i-1,j,k) 
      v1C  = q(3,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      v1R  = zero 
      v1RR = zero 

      v2LL = q(3,i-2,j,k) * csi(2,i-2,j,k) / aj(i-2,j,k) 
      v2L  = q(3,i-1,j,k) * csi(2,i-1,j,k) / aj(i-1,j,k) 
      v2C  = q(3,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      v2R  = zero 
      v2RR = zero 

      v3LL = q(3,i-2,j,k) * csi(3,i-2,j,k) / aj(i-2,j,k) 
      v3L  = q(3,i-1,j,k) * csi(3,i-1,j,k) / aj(i-1,j,k) 
      v3C  = q(3,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      v3R  = zero 
      v3RR = zero 

      !-------------------------------------------------------
   
      w1LL = q(4,i-2,j,k) * csi(1,i-2,j,k) / aj(i-2,j,k) 
      w1L  = q(4,i-1,j,k) * csi(1,i-1,j,k) / aj(i-1,j,k) 
      w1C  = q(4,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      w1R  = zero 
      w1RR = zero 

      w2LL = q(4,i-2,j,k) * csi(2,i-2,j,k) / aj(i-2,j,k) 
      w2L  = q(4,i-1,j,k) * csi(2,i-1,j,k) / aj(i-1,j,k) 
      w2C  = q(4,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      w2R  = zero 
      w2RR = zero 

      w3LL = q(4,i-2,j,k) * csi(3,i-2,j,k) / aj(i-2,j,k) 
      w3L  = q(4,i-1,j,k) * csi(3,i-1,j,k) / aj(i-1,j,k) 
      w3C  = q(4,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      w3R  = zero 
      w3RR = zero 

   else

      rLL      = rsign ( i-2 , j , k )
      rL       = rsign ( i-1 , j , k )
      rC       = rsign ( i   , j , k )
      rR       = rsign ( i+1 , j , k )
      rRR      = rsign ( i+2 , j , k )

      u1LL = q(2,i-2,j,k) * csi(1,i-2,j,k) / aj(i-2,j,k)                                    
      u1L  = q(2,i-1,j,k) * csi(1,i-1,j,k) / aj(i-1,j,k)                                    
      u1C  = q(2,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      u1R  = q(2,i+1,j,k) * csi(1,i+1,j,k) / aj(i+1,j,k) 
      u1RR = q(2,i+2,j,k) * csi(1,i+2,j,k) / aj(i+2,j,k) 

      u2LL = q(2,i-2,j,k) * csi(2,i-2,j,k) / aj(i-2,j,k) 
      u2L  = q(2,i-1,j,k) * csi(2,i-1,j,k) / aj(i-1,j,k) 
      u2C  = q(2,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      u2R  = q(2,i+1,j,k) * csi(2,i+1,j,k) / aj(i+1,j,k) 
      u2RR = q(2,i+2,j,k) * csi(2,i+2,j,k) / aj(i+2,j,k) 

      u3LL = q(2,i-2,j,k) * csi(3,i-2,j,k) / aj(i-2,j,k) 
      u3L  = q(2,i-1,j,k) * csi(3,i-1,j,k) / aj(i-1,j,k) 
      u3C  = q(2,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      u3R  = q(2,i+1,j,k) * csi(3,i+1,j,k) / aj(i+1,j,k) 
      u3RR = q(2,i+2,j,k) * csi(3,i+2,j,k) / aj(i+2,j,k) 

      !-------------------------------------------------------

      v1LL = q(3,i-2,j,k) * csi(1,i-2,j,k) / aj(i-2,j,k) 
      v1L  = q(3,i-1,j,k) * csi(1,i-1,j,k) / aj(i-1,j,k) 
      v1C  = q(3,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      v1R  = q(3,i+1,j,k) * csi(1,i+1,j,k) / aj(i+1,j,k) 
      v1RR = q(3,i+2,j,k) * csi(1,i+2,j,k) / aj(i+2,j,k) 

      v2LL = q(3,i-2,j,k) * csi(2,i-2,j,k) / aj(i-2,j,k) 
      v2L  = q(3,i-1,j,k) * csi(2,i-1,j,k) / aj(i-1,j,k) 
      v2C  = q(3,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      v2R  = q(3,i+1,j,k) * csi(2,i+1,j,k) / aj(i+1,j,k) 
      v2RR = q(3,i+2,j,k) * csi(2,i+2,j,k) / aj(i+2,j,k) 

      v3LL = q(3,i-2,j,k) * csi(3,i-2,j,k) / aj(i-2,j,k)
      v3L  = q(3,i-1,j,k) * csi(3,i-1,j,k) / aj(i-1,j,k) 
      v3C  = q(3,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      v3R  = q(3,i+1,j,k) * csi(3,i+1,j,k) / aj(i+1,j,k) 
      v3RR = q(3,i+2,j,k) * csi(3,i+2,j,k) / aj(i+2,j,k) 

      !-------------------------------------------------------
   
      w1LL = q(4,i-2,j,k) * csi(1,i-2,j,k) / aj(i-2,j,k) 
      w1L  = q(4,i-1,j,k) * csi(1,i-1,j,k) / aj(i-1,j,k) 
      w1C  = q(4,i  ,j,k) * csi(1,i  ,j,k) / aj(i  ,j,k) 
      w1R  = q(4,i+1,j,k) * csi(1,i+1,j,k) / aj(i+1,j,k) 
      w1RR = q(4,i+2,j,k) * csi(1,i+2,j,k) / aj(i+2,j,k) 

      w2LL = q(4,i-2,j,k) * csi(2,i-2,j,k) / aj(i-2,j,k) 
      w2L  = q(4,i-1,j,k) * csi(2,i-1,j,k) / aj(i-1,j,k) 
      w2C  = q(4,i  ,j,k) * csi(2,i  ,j,k) / aj(i  ,j,k) 
      w2R  = q(4,i+1,j,k) * csi(2,i+1,j,k) / aj(i+1,j,k) 
      w2RR = q(4,i+2,j,k) * csi(2,i+2,j,k) / aj(i+2,j,k) 

      w3LL = q(4,i-2,j,k) * csi(3,i-2,j,k) / aj(i-2,j,k) 
      w3L  = q(4,i-1,j,k) * csi(3,i-1,j,k) / aj(i-1,j,k) 
      w3C  = q(4,i  ,j,k) * csi(3,i  ,j,k) / aj(i  ,j,k) 
      w3R  = q(4,i+1,j,k) * csi(3,i+1,j,k) / aj(i+1,j,k) 
      w3RR = q(4,i+2,j,k) * csi(3,i+2,j,k) / aj(i+2,j,k) 


   end if

   exsign = zero

   ! ∂u/∂ξ
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     u1LL, u1L , u1C, u1R, u1RR   , &
                                     dc                           , &
                                     bias_csi                     , &
                                     exsign                       , &
                                     du1_dcsi                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     u2LL, u2L , u2C, u2R, u2RR   , &
                                     dc                           , &
                                     bias_csi                     , &
                                     exsign                       , &
                                     du2_dcsi                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     u3LL, u3L , u3C, u3R, u3RR   , &
                                     dc                           , &
                                     bias_csi                     , &
                                     exsign                       , &
                                     du3_dcsi                       &
                                    )

   if ( exsign < one_half ) return                                    
                                    


   ! ∂v/∂ξ
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     v1LL, v1L , v1C, v1R, v1RR   , &
                                     dc                           , &
                                     bias_csi                     , &
                                     exsign                       , &
                                     dv1_dcsi                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     v2LL, v2L , v2C, v2R, v2RR   , &
                                     dc                           , &
                                     bias_csi                     , &
                                     exsign                       , &
                                     dv2_dcsi                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     v3LL, v3L , v3C, v3R, v3RR   , &
                                     dc                           , &
                                     bias_csi                     , &
                                     exsign                       , &
                                     dv3_dcsi                       &
                                    )

   if ( exsign < one_half ) return   

   ! ∂w/∂ξ
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     w1LL, w1L , w1C, w1R, w1RR   , &
                                     dc                           , &
                                     bias_csi                     , &
                                     exsign                       , &
                                     dw1_dcsi                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     w2LL, w2L , w2C, w2R, w2RR   , &
                                     dc                           , &
                                     bias_csi                     , &
                                     exsign                       , &
                                     dw2_dcsi                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     w3LL, w3L , w3C, w3R, w3RR   , &
                                     dc                           , &
                                     bias_csi                     , &
                                     exsign                       , &
                                     dw3_dcsi                       &
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

      u4LL = zero                                    
      u4L  = zero                                    
      u4C  = q(2,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      u4R  = q(2,i,j+1,k) * eta(1,i,j+1,k) / aj(i,j+1,k) 
      u4RR = q(2,i,j+2,k) * eta(1,i,j+2,k) / aj(i,j+2,k) 

      u5LL = zero 
      u5L  = zero 
      u5C  = q(2,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      u5R  = q(2,i,j+1,k) * eta(2,i,j+1,k) / aj(i,j+1,k) 
      u5RR = q(2,i,j+2,k) * eta(2,i,j+2,k) / aj(i,j+2,k) 

      u6LL = zero 
      u6L  = zero 
      u6C  = q(2,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      u6R  = q(2,i,j+1,k) * eta(3,i,j+1,k) / aj(i,j+1,k) 
      u6RR = q(2,i,j+2,k) * eta(3,i,j+2,k) / aj(i,j+2,k) 

      !-------------------------------------------------------

      v4LL = zero 
      v4L  = zero 
      v4C  = q(3,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      v4R  = q(3,i,j+1,k) * eta(1,i,j+1,k) / aj(i,j+1,k) 
      v4RR = q(3,i,j+2,k) * eta(1,i,j+2,k) / aj(i,j+2,k) 

      v5LL = zero 
      v5L  = zero 
      v5C  = q(3,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      v5R  = q(3,i,j+1,k) * eta(2,i,j+1,k) / aj(i,j+1,k) 
      v5RR = q(3,i,j+2,k) * eta(2,i,j+2,k) / aj(i,j+2,k) 

      v6LL = q(3,i,j-2,k) * eta(3,i,j-2,k) / aj(i,j-2,k) 
      v6L  = q(3,i,j-1,k) * eta(3,i,j-1,k) / aj(i,j-1,k) 
      v6C  = q(3,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      v6R  = q(3,i,j+1,k) * eta(3,i,j+1,k) / aj(i,j+1,k) 
      v6RR = q(3,i,j+2,k) * eta(3,i,j+2,k) / aj(i,j+2,k) 

      !-------------------------------------------------------
   
      w4LL = zero 
      w4L  = zero 
      w4C  = q(4,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      w4R  = q(4,i,j+1,k) * eta(1,i,j+1,k) / aj(i,j+1,k) 
      w4RR = q(4,i,j+2,k) * eta(1,i,j+2,k) / aj(i,j+2,k) 

      w5LL = zero 
      w5L  = zero 
      w5C  = q(4,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      w5R  = q(4,i,j+1,k) * eta(2,i,j+1,k) / aj(i,j+1,k) 
      w5RR = q(4,i,j+2,k) * eta(2,i,j+2,k) / aj(i,j+2,k) 

      w6LL = zero 
      w6L  = zero 
      w6C  = q(4,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      w6R  = q(4,i,j+1,k) * eta(3,i,j+1,k) / aj(i,j+1,k) 
      w6RR = q(4,i,j+2,k) * eta(3,i,j+2,k) / aj(i,j+2,k) 


   else if ( j == jsta + 1 ) then
         
      rLL      = zero ! dummy
      rL       = rsign ( i , j-1 , k )
      rC       = rsign ( i , j   , k )
      rR       = rsign ( i , j+1 , k )
      rRR      = rsign ( i , j+2 , k )

      u4LL = zero                                    
      u4L  = q(2,i,j-1,k) * eta(1,i,j-1,k) / aj(i,j-1,k)                                    
      u4C  = q(2,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      u4R  = q(2,i,j+1,k) * eta(1,i,j+1,k) / aj(i,j+1,k) 
      u4RR = q(2,i,j+2,k) * eta(1,i,j+2,k) / aj(i,j+2,k) 

      u5LL = zero 
      u5L  = q(2,i,j-1,k) * eta(2,i,j-1,k) / aj(i,j-1,k) 
      u5C  = q(2,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      u5R  = q(2,i,j+1,k) * eta(2,i,j+1,k) / aj(i,j+1,k) 
      u5RR = q(2,i,j+2,k) * eta(2,i,j+2,k) / aj(i,j+2,k) 

      u6LL = zero 
      u6L  = q(2,i,j-1,k) * eta(3,i,j-1,k) / aj(i,j-1,k) 
      u6C  = q(2,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      u6R  = q(2,i,j+1,k) * eta(3,i,j+1,k) / aj(i,j+1,k) 
      u6RR = q(2,i,j+2,k) * eta(3,i,j+2,k) / aj(i,j+2,k) 

      !-------------------------------------------------------

      v4LL = zero 
      v4L  = q(3,i,j-1,k) * eta(1,i,j-1,k) / aj(i,j-1,k) 
      v4C  = q(3,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      v4R  = q(3,i,j+1,k) * eta(1,i,j+1,k) / aj(i,j+1,k) 
      v4RR = q(3,i,j+2,k) * eta(1,i,j+2,k) / aj(i,j+2,k) 

      v5LL = zero 
      v5L  = q(3,i,j-1,k) * eta(2,i,j-1,k) / aj(i,j-1,k) 
      v5C  = q(3,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      v5R  = q(3,i,j+1,k) * eta(2,i,j+1,k) / aj(i,j+1,k) 
      v5RR = q(3,i,j+2,k) * eta(2,i,j+2,k) / aj(i,j+2,k) 

      v6LL = zero 
      v6L  = q(3,i,j-1,k) * eta(3,i,j-1,k) / aj(i,j-1,k) 
      v6C  = q(3,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      v6R  = q(3,i,j+1,k) * eta(3,i,j+1,k) / aj(i,j+1,k) 
      v6RR = q(3,i,j+2,k) * eta(3,i,j+2,k) / aj(i,j+2,k) 

      !-------------------------------------------------------
   
      w4LL = zero 
      w4L  = q(4,i,j-1,k) * eta(1,i,j-1,k) / aj(i,j-1,k) 
      w4C  = q(4,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      w4R  = q(4,i,j+1,k) * eta(1,i,j+1,k) / aj(i,j+1,k) 
      w4RR = q(4,i,j+2,k) * eta(1,i,j+2,k) / aj(i,j+2,k) 

      w5LL = zero 
      w5L  = q(4,i,j-1,k) * eta(2,i,j-1,k) / aj(i,j-1,k) 
      w5C  = q(4,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      w5R  = q(4,i,j+1,k) * eta(2,i,j+1,k) / aj(i,j+1,k) 
      w5RR = q(4,i,j+2,k) * eta(2,i,j+2,k) / aj(i,j+2,k) 

      w6LL = zero 
      w6L  = q(4,i,j-1,k) * eta(3,i,j-1,k) / aj(i,j-1,k) 
      w6C  = q(4,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      w6R  = q(4,i,j+1,k) * eta(3,i,j+1,k) / aj(i,j+1,k) 
      w6RR = q(4,i,j+2,k) * eta(3,i,j+2,k) / aj(i,j+2,k) 

   
   else if ( j == jend - 1 ) then
         
      rLL      = rsign ( i , j-2 , k )
      rL       = rsign ( i , j-1 , k )
      rC       = rsign ( i , j   , k )
      rR       = rsign ( i , j+1 , k )
      rRR      = zero ! dummy

      u4LL = q(2,i,j-2,k) * eta(1,i,j-2,k) / aj(i,j-2,k)                                    
      u4L  = q(2,i,j-1,k) * eta(1,i,j-1,k) / aj(i,j-1,k)                                    
      u4C  = q(2,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      u4R  = q(2,i,j+1,k) * eta(1,i,j+1,k) / aj(i,j+1,k) 
      u4RR = zero

      u5LL = q(2,i,j-2,k) * eta(2,i,j-2,k) / aj(i,j-2,k) 
      u5L  = q(2,i,j-1,k) * eta(2,i,j-1,k) / aj(i,j-1,k) 
      u5C  = q(2,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      u5R  = q(2,i,j+1,k) * eta(2,i,j+1,k) / aj(i,j+1,k) 
      u5RR = zero

      u6LL = q(2,i,j-2,k) * eta(3,i,j-2,k) / aj(i,j-2,k) 
      u6L  = q(2,i,j-1,k) * eta(3,i,j-1,k) / aj(i,j-1,k) 
      u6C  = q(2,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      u6R  = q(2,i,j+1,k) * eta(3,i,j+1,k) / aj(i,j+1,k) 
      u6RR = zero 

      !-------------------------------------------------------

      v4LL = q(3,i,j-2,k) * eta(1,i,j-2,k) / aj(i,j-2,k) 
      v4L  = q(3,i,j-1,k) * eta(1,i,j-1,k) / aj(i,j-1,k) 
      v4C  = q(3,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      v4R  = q(3,i,j+1,k) * eta(1,i,j+1,k) / aj(i,j+1,k) 
      v4RR = zero

      v5LL = q(3,i,j-2,k) * eta(2,i,j-2,k) / aj(i,j-2,k) 
      v5L  = q(3,i,j-1,k) * eta(2,i,j-1,k) / aj(i,j-1,k) 
      v5C  = q(3,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      v5R  = q(3,i,j+1,k) * eta(2,i,j+1,k) / aj(i,j+1,k) 
      v5RR = zero

      v6LL = q(3,i,j-2,k) * eta(3,i,j-2,k) / aj(i,j-2,k) 
      v6L  = q(3,i,j-1,k) * eta(3,i,j-1,k) / aj(i,j-1,k) 
      v6C  = q(3,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      v6R  = q(3,i,j+1,k) * eta(3,i,j+1,k) / aj(i,j+1,k) 
      v6RR = zero

      !-------------------------------------------------------
   
      w4LL = q(4,i,j-2,k) * eta(1,i,j-2,k) / aj(i,j-2,k) 
      w4L  = q(4,i,j-1,k) * eta(1,i,j-1,k) / aj(i,j-1,k) 
      w4C  = q(4,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      w4R  = q(4,i,j+1,k) * eta(1,i,j+1,k) / aj(i,j+1,k) 
      w4RR = zero

      w5LL = q(4,i,j-2,k) * eta(2,i,j-2,k) / aj(i,j-2,k) 
      w5L  = q(4,i,j-1,k) * eta(2,i,j-1,k) / aj(i,j-1,k) 
      w5C  = q(4,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      w5R  = q(4,i,j+1,k) * eta(2,i,j+1,k) / aj(i,j+1,k) 
      w5RR = zero

      w6LL = q(4,i,j-2,k) * eta(3,i,j-2,k) / aj(i,j-2,k) 
      w6L  = q(4,i,j-1,k) * eta(3,i,j-1,k) / aj(i,j-1,k) 
      w6C  = q(4,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      w6R  = q(4,i,j+1,k) * eta(3,i,j+1,k) / aj(i,j+1,k) 
      w6RR = zero

   else if ( j == jend ) then
         
      bias_eta = -1

      rLL      = rsign ( i , j-2 , k )
      rL       = rsign ( i , j-1 , k )
      rC       = rsign ( i , j   , k )
      rR       = zero ! dummy
      rRR      = zero ! dummy

      u4LL = q(2,i,j-2,k) * eta(1,i,j-2,k) / aj(i,j-2,k)                                    
      u4L  = q(2,i,j-1,k) * eta(1,i,j-1,k) / aj(i,j-1,k)                                    
      u4C  = q(2,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      u4R  = zero 
      u4RR = zero 

      u5LL = q(2,i,j-2,k) * eta(2,i,j-2,k) / aj(i,j-2,k) 
      u5L  = q(2,i,j-1,k) * eta(2,i,j-1,k) / aj(i,j-1,k) 
      u5C  = q(2,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      u5R  = zero 
      u5RR = zero 

      u6LL = q(2,i,j-2,k) * eta(3,i,j-2,k) / aj(i,j-2,k) 
      u6L  = q(2,i,j-1,k) * eta(3,i,j-1,k) / aj(i,j-1,k) 
      u6C  = q(2,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      u6R  = zero 
      u6RR = zero 

      !-------------------------------------------------------

      v4LL = q(3,i,j-2,k) * eta(1,i,j-2,k) / aj(i,j-2,k) 
      v4L  = q(3,i,j-1,k) * eta(1,i,j-1,k) / aj(i,j-1,k) 
      v4C  = q(3,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      v4R  = zero 
      v4RR = zero 

      v5LL = q(3,i,j-2,k) * eta(2,i,j-2,k) / aj(i,j-2,k) 
      v5L  = q(3,i,j-1,k) * eta(2,i,j-1,k) / aj(i,j-1,k) 
      v5C  = q(3,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      v5R  = zero 
      v5RR = zero 

      v6LL = q(3,i,j-2,k) * eta(3,i,j-2,k) / aj(i,j-2,k) 
      v6L  = q(3,i,j-1,k) * eta(3,i,j-1,k) / aj(i,j-1,k) 
      v6C  = q(3,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      v6R  = zero 
      v6RR = zero 

      !-------------------------------------------------------
   
      w4LL = q(4,i,j-2,k) * eta(1,i,j-2,k) / aj(i,j-2,k) 
      w4L  = q(4,i,j-1,k) * eta(1,i,j-1,k) / aj(i,j-1,k) 
      w4C  = q(4,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      w4R  = zero 
      w4RR = zero 

      w5LL = q(4,i,j-2,k) * eta(2,i,j-2,k) / aj(i,j-2,k) 
      w5L  = q(4,i,j-1,k) * eta(2,i,j-1,k) / aj(i,j-1,k) 
      w5C  = q(4,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      w5R  = zero 
      w5RR = zero 

      w6LL = q(4,i,j-2,k) * eta(3,i,j-2,k) / aj(i,j-2,k) 
      w6L  = q(4,i,j-1,k) * eta(3,i,j-1,k) / aj(i,j-1,k) 
      w6C  = q(4,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      w6R  = zero 
      w6RR = zero 

   else

      rLL      = rsign ( i , j-2 , k )
      rL       = rsign ( i , j-1 , k )
      rC       = rsign ( i , j   , k )
      rR       = rsign ( i , j+1 , k )
      rRR      = rsign ( i , j+2 , k )

      u4LL = q(2,i,j-2,k) * eta(1,i,j-2,k) / aj(i,j-2,k)                                    
      u4L  = q(2,i,j-1,k) * eta(1,i,j-1,k) / aj(i,j-1,k)                                    
      u4C  = q(2,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      u4R  = q(2,i,j+1,k) * eta(1,i,j+1,k) / aj(i,j+1,k) 
      u4RR = q(2,i,j+2,k) * eta(1,i,j+2,k) / aj(i,j+2,k) 

      u5LL = q(2,i,j-2,k) * eta(2,i,j-2,k) / aj(i,j-2,k) 
      u5L  = q(2,i,j-1,k) * eta(2,i,j-1,k) / aj(i,j-1,k) 
      u5C  = q(2,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      u5R  = q(2,i,j+1,k) * eta(2,i,j+1,k) / aj(i,j+1,k) 
      u5RR = q(2,i,j+2,k) * eta(2,i,j+2,k) / aj(i,j+2,k) 

      u6LL = q(2,i,j-2,k) * eta(3,i,j-2,k) / aj(i,j-2,k) 
      u6L  = q(2,i,j-1,k) * eta(3,i,j-1,k) / aj(i,j-1,k) 
      u6C  = q(2,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      u6R  = q(2,i,j+1,k) * eta(3,i,j+1,k) / aj(i,j+1,k) 
      u6RR = q(2,i,j+2,k) * eta(3,i,j+2,k) / aj(i,j+2,k) 

      !-------------------------------------------------------

      v4LL = q(3,i,j-2,k) * eta(1,i,j-2,k) / aj(i,j-2,k) 
      v4L  = q(3,i,j-1,k) * eta(1,i,j-1,k) / aj(i,j-1,k) 
      v4C  = q(3,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      v4R  = q(3,i,j+1,k) * eta(1,i,j+1,k) / aj(i,j+1,k) 
      v4RR = q(3,i,j+2,k) * eta(1,i,j+2,k) / aj(i,j+2,k) 

      v5LL = q(3,i,j-2,k) * eta(2,i,j-2,k) / aj(i,j-2,k) 
      v5L  = q(3,i,j-1,k) * eta(2,i,j-1,k) / aj(i,j-1,k) 
      v5C  = q(3,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      v5R  = q(3,i,j+1,k) * eta(2,i,j+1,k) / aj(i,j+1,k) 
      v5RR = q(3,i,j+2,k) * eta(2,i,j+2,k) / aj(i,j+2,k) 

      v6LL = q(3,i,j-2,k) * eta(3,i,j-2,k) / aj(i,j-2,k) 
      v6L  = q(3,i,j-1,k) * eta(3,i,j-1,k) / aj(i,j-1,k) 
      v6C  = q(3,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      v6R  = q(3,i,j+1,k) * eta(3,i,j+1,k) / aj(i,j+1,k) 
      v6RR = q(3,i,j+2,k) * eta(3,i,j+2,k) / aj(i,j+2,k) 

      !-------------------------------------------------------
   
      w4LL = q(4,i,j-2,k) * eta(1,i,j-2,k) / aj(i,j-2,k) 
      w4L  = q(4,i,j-1,k) * eta(1,i,j-1,k) / aj(i,j-1,k) 
      w4C  = q(4,i,j  ,k) * eta(1,i,j  ,k) / aj(i,j  ,k) 
      w4R  = q(4,i,j+1,k) * eta(1,i,j+1,k) / aj(i,j+1,k) 
      w4RR = q(4,i,j+2,k) * eta(1,i,j+2,k) / aj(i,j+2,k) 

      w5LL = q(4,i,j-2,k) * eta(2,i,j-2,k) / aj(i,j-2,k) 
      w5L  = q(4,i,j-1,k) * eta(2,i,j-1,k) / aj(i,j-1,k) 
      w5C  = q(4,i,j  ,k) * eta(2,i,j  ,k) / aj(i,j  ,k) 
      w5R  = q(4,i,j+1,k) * eta(2,i,j+1,k) / aj(i,j+1,k) 
      w5RR = q(4,i,j+2,k) * eta(2,i,j+2,k) / aj(i,j+2,k) 

      w6LL = q(4,i,j-2,k) * eta(3,i,j-2,k) / aj(i,j-2,k) 
      w6L  = q(4,i,j-1,k) * eta(3,i,j-1,k) / aj(i,j-1,k) 
      w6C  = q(4,i,j  ,k) * eta(3,i,j  ,k) / aj(i,j  ,k) 
      w6R  = q(4,i,j+1,k) * eta(3,i,j+1,k) / aj(i,j+1,k) 
      w6RR = q(4,i,j+2,k) * eta(3,i,j+2,k) / aj(i,j+2,k) 


   end if

   exsign = zero

   ! ∂u/∂η
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     u4LL, u4L , u4C, u4R, u4RR   , &
                                     de                           , &
                                     bias_eta                     , &
                                     exsign                       , &
                                     du4_deta                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     u5LL, u5L , u5C, u5R, u5RR   , &
                                     de                           , &
                                     bias_eta                     , &
                                     exsign                       , &
                                     du5_deta                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     u6LL, u6L , u6C, u6R, u6RR   , &
                                     de                           , &
                                     bias_eta                     , &
                                     exsign                       , &
                                     du6_deta                       &
                                    )

   if ( exsign < one_half ) return                                    
                                    


   ! ∂v/∂η
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     v4LL, v4L , v4C, v4R, v4RR   , &
                                     de                           , &
                                     bias_eta                     , &
                                     exsign                       , &
                                     dv4_deta                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     v5LL, v5L , v5C, v5R, v5RR   , &
                                     de                           , &
                                     bias_eta                     , &
                                     exsign                       , &
                                     dv5_deta                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     v6LL, v6L , v6C, v6R, v6RR   , &
                                     de                           , &
                                     bias_eta                     , &
                                     exsign                       , &
                                     dv6_deta                       &
                                    )

   if ( exsign < one_half ) return     

   ! ∂w/∂η
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     w4LL, w4L , w4C, w4R, w4RR   , &
                                     de                           , &
                                     bias_eta                     , &
                                     exsign                       , &
                                     dw4_deta                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     w5LL, w5L , w5C, w5R, w5RR   , &
                                     de                           , &
                                     bias_eta                     , &
                                     exsign                       , &
                                     dw5_deta                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     w6LL, w6L , w6C, w6R, w6RR   , &
                                     de                           , &
                                     bias_eta                     , &
                                     exsign                       , &
                                     dw6_deta                       &
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

      u7LL = zero                                    
      u7L  = zero                                    
      u7C  = q(2,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      u7R  = q(2,i,j,k+1) * zet(1,i,j,k+1) / aj(i,j,k+1) 
      u7RR = q(2,i,j,k+2) * zet(1,i,j,k+2) / aj(i,j,k+2) 

      u8LL = zero 
      u8L  = zero 
      u8C  = q(2,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      u8R  = q(2,i,j,k+1) * zet(2,i,j,k+1) / aj(i,j,k+1) 
      u8RR = q(2,i,j,k+2) * zet(2,i,j,k+2) / aj(i,j,k+2) 

      u9LL = zero 
      u9L  = zero 
      u9C  = q(2,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      u9R  = q(2,i,j,k+1) * zet(3,i,j,k+1) / aj(i,j,k+1) 
      u9RR = q(2,i,j,k+2) * zet(3,i,j,k+2) / aj(i,j,k+2) 

      !-------------------------------------------------------

      v7LL = zero 
      v7L  = zero 
      v7C  = q(3,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      v7R  = q(3,i,j,k+1) * zet(1,i,j,k+1) / aj(i,j,k+1) 
      v7RR = q(3,i,j,k+2) * zet(1,i,j,k+2) / aj(i,j,k+2) 

      v8LL = zero 
      v8L  = zero 
      v8C  = q(3,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      v8R  = q(3,i,j,k+1) * zet(2,i,j,k+1) / aj(i,j,k+1) 
      v8RR = q(3,i,j,k+2) * zet(2,i,j,k+2) / aj(i,j,k+2) 

      v9LL = zero 
      v9L  = zero 
      v9C  = q(3,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      v9R  = q(3,i,j,k+1) * zet(3,i,j,k+1) / aj(i,j,k+1) 
      v9RR = q(3,i,j,k+2) * zet(3,i,j,k+2) / aj(i,j,k+2) 

      !-------------------------------------------------------
   
      w7LL = zero 
      w7L  = zero 
      w7C  = q(4,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      w7R  = q(4,i,j,k+1) * zet(1,i,j,k+1) / aj(i,j,k+1) 
      w7RR = q(4,i,j,k+2) * zet(1,i,j,k+2) / aj(i,j,k+2) 

      w8LL = zero 
      w8L  = zero 
      w8C  = q(4,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      w8R  = q(4,i,j,k+1) * zet(2,i,j,k+1) / aj(i,j,k+1) 
      w8RR = q(4,i,j,k+2) * zet(2,i,j,k+2) / aj(i,j,k+2) 

      w9LL = zero 
      w9L  = zero 
      w9C  = q(4,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      w9R  = q(4,i,j,k+1) * zet(3,i,j,k+1) / aj(i,j,k+1) 
      w9RR = q(4,i,j,k+2) * zet(3,i,j,k+2) / aj(i,j,k+2) 



   else if ( k == ksta + 1 ) then
         
      rLL      = zero ! dummy
      rL       = rsign ( i , j , k-1 )
      rC       = rsign ( i , j , k   )
      rR       = rsign ( i , j , k+1 )
      rRR      = rsign ( i , j , k+2 )

      u7LL = zero                                    
      u7L  = q(2,i,j,k-1) * zet(1,i,j,k-1) / aj(i,j,k-1)                                    
      u7C  = q(2,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      u7R  = q(2,i,j,k+1) * zet(1,i,j,k+1) / aj(i,j,k+1) 
      u7RR = q(2,i,j,k+2) * zet(1,i,j,k+2) / aj(i,j,k+2) 

      u8LL = zero 
      u8L  = q(2,i,j,k-1) * zet(2,i,j,k-1) / aj(i,j,k-1) 
      u8C  = q(2,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      u8R  = q(2,i,j,k+1) * zet(2,i,j,k+1) / aj(i,j,k+1) 
      u8RR = q(2,i,j,k+2) * zet(2,i,j,k+2) / aj(i,j,k+2) 

      u9LL = zero 
      u9L  = q(2,i,j,k-1) * zet(3,i,j,k-1) / aj(i,j,k-1) 
      u9C  = q(2,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      u9R  = q(2,i,j,k+1) * zet(3,i,j,k+1) / aj(i,j,k+1) 
      u9RR = q(2,i,j,k+2) * zet(3,i,j,k+2) / aj(i,j,k+2) 

      !-------------------------------------------------------

      v7LL = zero 
      v7L  = q(3,i,j,k-1) * zet(1,i,j,k-1) / aj(i,j,k-1) 
      v7C  = q(3,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      v7R  = q(3,i,j,k+1) * zet(1,i,j,k+1) / aj(i,j,k+1) 
      v7RR = q(3,i,j,k+2) * zet(1,i,j,k+2) / aj(i,j,k+2) 

      v8LL = zero 
      v8L  = q(3,i,j,k-1) * zet(2,i,j,k-1) / aj(i,j,k-1) 
      v8C  = q(3,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      v8R  = q(3,i,j,k+1) * zet(2,i,j,k+1) / aj(i,j,k+1) 
      v8RR = q(3,i,j,k+2) * zet(2,i,j,k+2) / aj(i,j,k+2) 

      v9LL = zero 
      v9L  = q(3,i,j,k-1) * zet(3,i,j,k-1) / aj(i,j,k-1) 
      v9C  = q(3,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      v9R  = q(3,i,j,k+1) * zet(3,i,j,k+1) / aj(i,j,k+1) 
      v9RR = q(3,i,j,k+2) * zet(3,i,j,k+2) / aj(i,j,k+2) 

      !-------------------------------------------------------
   
      w7LL = zero 
      w7L  = q(4,i,j,k-1) * zet(1,i,j,k-1) / aj(i,j,k-1) 
      w7C  = q(4,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      w7R  = q(4,i,j,k+1) * zet(1,i,j,k+1) / aj(i,j,k+1) 
      w7RR = q(4,i,j,k+2) * zet(1,i,j,k+2) / aj(i,j,k+2) 

      w8LL = zero 
      w8L  = q(4,i,j,k-1) * zet(2,i,j,k-1) / aj(i,j,k-1) 
      w8C  = q(4,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      w8R  = q(4,i,j,k+1) * zet(2,i,j,k+1) / aj(i,j,k+1) 
      w8RR = q(4,i,j,k+2) * zet(2,i,j,k+2) / aj(i,j,k+2) 

      w9LL = zero 
      w9L  = q(4,i,j,k-1) * zet(3,i,j,k-1) / aj(i,j,k-1) 
      w9C  = q(4,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      w9R  = q(4,i,j,k+1) * zet(3,i,j,k+1) / aj(i,j,k+1) 
      w9RR = q(4,i,j,k+2) * zet(3,i,j,k+2) / aj(i,j,k+2) 

   
   else if ( k == kend - 1 ) then
         
      rLL      = rsign ( i , j , k-2 )
      rL       = rsign ( i , j , k-1 )
      rC       = rsign ( i , j , k   )
      rR       = rsign ( i , j , k+1 )
      rRR      = zero ! dummy

      u7LL = q(2,i,j,k-2) * zet(1,i,j,k-2) / aj(i,j,k-2)                                    
      u7L  = q(2,i,j,k-1) * zet(1,i,j,k-1) / aj(i,j,k-1)                                    
      u7C  = q(2,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      u7R  = q(2,i,j,k+1) * zet(1,i,j,k+1) / aj(i,j,k+1) 
      u7RR = zero 

      u8LL = q(2,i,j,k-2) * zet(2,i,j,k-2) / aj(i,j,k-2) 
      u8L  = q(2,i,j,k-1) * zet(2,i,j,k-1) / aj(i,j,k-1) 
      u8C  = q(2,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      u8R  = q(2,i,j,k+1) * zet(2,i,j,k+1) / aj(i,j,k+1) 
      u8RR = zero 

      u9LL = q(2,i,j,k-2) * zet(3,i,j,k-2) / aj(i,j,k-2) 
      u9L  = q(2,i,j,k-1) * zet(3,i,j,k-1) / aj(i,j,k-1) 
      u9C  = q(2,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      u9R  = q(2,i,j,k+1) * zet(3,i,j,k+1) / aj(i,j,k+1) 
      u9RR = zero 

      !-------------------------------------------------------

      v7LL = q(3,i,j,k-2) * zet(1,i,j,k-2) / aj(i,j,k-2) 
      v7L  = q(3,i,j,k-1) * zet(1,i,j,k-1) / aj(i,j,k-1) 
      v7C  = q(3,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      v7R  = q(3,i,j,k+1) * zet(1,i,j,k+1) / aj(i,j,k+1) 
      v7RR = zero 

      v8LL = q(3,i,j,k-2) * zet(2,i,j,k-2) / aj(i,j,k-2) 
      v8L  = q(3,i,j,k-1) * zet(2,i,j,k-1) / aj(i,j,k-1) 
      v8C  = q(3,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      v8R  = q(3,i,j,k+1) * zet(2,i,j,k+1) / aj(i,j,k+1) 
      v8RR = zero 

      v9LL = q(3,i,j,k-2) * zet(3,i,j,k-2) / aj(i,j,k-2) 
      v9L  = q(3,i,j,k-1) * zet(3,i,j,k-1) / aj(i,j,k-1) 
      v9C  = q(3,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      v9R  = q(3,i,j,k+1) * zet(3,i,j,k+1) / aj(i,j,k+1) 
      v9RR = zero

      !-------------------------------------------------------
   
      w7LL = q(4,i,j,k-2) * zet(1,i,j,k-2) / aj(i,j,k-2) 
      w7L  = q(4,i,j,k-1) * zet(1,i,j,k-1) / aj(i,j,k-1) 
      w7C  = q(4,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      w7R  = q(4,i,j,k+1) * zet(1,i,j,k+1) / aj(i,j,k+1) 
      w7RR = zero 

      w8LL = q(4,i,j,k-2) * zet(2,i,j,k-2) / aj(i,j,k-2) 
      w8L  = q(4,i,j,k-1) * zet(2,i,j,k-1) / aj(i,j,k-1) 
      w8C  = q(4,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      w8R  = q(4,i,j,k+1) * zet(2,i,j,k+1) / aj(i,j,k+1) 
      w8RR = zero 

      w9LL = q(4,i,j,k-2) * zet(3,i,j,k-2) / aj(i,j,k-2) 
      w9L  = q(4,i,j,k-1) * zet(3,i,j,k-1) / aj(i,j,k-1) 
      w9C  = q(4,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      w9R  = q(4,i,j,k+1) * zet(3,i,j,k+1) / aj(i,j,k+1) 
      w9RR = zero 

   else if ( k == kend ) then
         
      bias_zet = -1

      rLL      = rsign ( i , j , k-2 )
      rL       = rsign ( i , j , k-1 )
      rC       = rsign ( i , j , k   )
      rR       = zero ! dummy
      rRR      = zero ! dummy

      u7LL = q(2,i,j,k-2) * zet(1,i,j,k-2) / aj(i,j,k-2)                                    
      u7L  = q(2,i,j,k-1) * zet(1,i,j,k-1) / aj(i,j,k-1)                                    
      u7C  = q(2,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      u7R  = zero 
      u7RR = zero 

      u8LL = q(2,i,j,k-2) * zet(2,i,j,k-2) / aj(i,j,k-2) 
      u8L  = q(2,i,j,k-1) * zet(2,i,j,k-1) / aj(i,j,k-1) 
      u8C  = q(2,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      u8R  = zero 
      u8RR = zero 

      u9LL = q(2,i,j,k-2) * zet(3,i,j,k-2) / aj(i,j,k-2) 
      u9L  = q(2,i,j,k-1) * zet(3,i,j,k-1) / aj(i,j,k-1) 
      u9C  = q(2,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      u9R  = zero 
      u9RR = zero 

      !-------------------------------------------------------

      v7LL = q(3,i,j,k-2) * zet(1,i,j,k-2) / aj(i,j,k-2) 
      v7L  = q(3,i,j,k-1) * zet(1,i,j,k-1) / aj(i,j,k-1) 
      v7C  = q(3,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      v7R  = zero 
      v7RR = zero 

      v8LL = q(3,i,j,k-2) * zet(2,i,j,k-2) / aj(i,j,k-2) 
      v8L  = q(3,i,j,k-1) * zet(2,i,j,k-1) / aj(i,j,k-1) 
      v8C  = q(3,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      v8R  = zero 
      v8RR = zero 

      v9LL = q(3,i,j,k-2) * zet(3,i,j,k-2) / aj(i,j,k-2) 
      v9L  = q(3,i,j,k-1) * zet(3,i,j,k-1) / aj(i,j,k-1) 
      v9C  = q(3,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      v9R  = zero 
      v9RR = zero 

      !-------------------------------------------------------
   
      w7LL = q(4,i,j,k-2) * zet(1,i,j,k-2) / aj(i,j,k-2) 
      w7L  = q(4,i,j,k-1) * zet(1,i,j,k-1) / aj(i,j,k-1) 
      w7C  = q(4,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      w7R  = zero 
      w7RR = zero 

      w8LL = q(4,i,j,k-2) * zet(2,i,j,k-2) / aj(i,j,k-2) 
      w8L  = q(4,i,j,k-1) * zet(2,i,j,k-1) / aj(i,j,k-1) 
      w8C  = q(4,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      w8R  = zero 
      w8RR = zero 

      w9LL = q(4,i,j,k-2) * zet(3,i,j,k-2) / aj(i,j,k-2) 
      w9L  = q(4,i,j,k-1) * zet(3,i,j,k-1) / aj(i,j,k-1) 
      w9C  = q(4,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      w9R  = zero 
      w9RR = zero 

   else

      rLL      = rsign ( i , j , k-2 )
      rL       = rsign ( i , j , k-1 )
      rC       = rsign ( i , j , k   )
      rR       = rsign ( i , j , k+1 )
      rRR      = rsign ( i , j , k+2 )

      u7LL = q(2,i,j,k-2) * zet(1,i,j,k-2) / aj(i,j,k-2)                                    
      u7L  = q(2,i,j,k-1) * zet(1,i,j,k-1) / aj(i,j,k-1)                                    
      u7C  = q(2,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      u7R  = q(2,i,j,k+1) * zet(1,i,j,k+1) / aj(i,j,k+1) 
      u7RR = q(2,i,j,k+2) * zet(1,i,j,k+2) / aj(i,j,k+2) 

      u8LL = q(2,i,j,k-2) * zet(2,i,j,k-2) / aj(i,j,k-2) 
      u8L  = q(2,i,j,k-1) * zet(2,i,j,k-1) / aj(i,j,k-1) 
      u8C  = q(2,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      u8R  = q(2,i,j,k+1) * zet(2,i,j,k+1) / aj(i,j,k+1) 
      u8RR = q(2,i,j,k+2) * zet(2,i,j,k+2) / aj(i,j,k+2) 

      u9LL = q(2,i,j,k-2) * zet(3,i,j,k-2) / aj(i,j,k-2) 
      u9L  = q(2,i,j,k-1) * zet(3,i,j,k-1) / aj(i,j,k-1) 
      u9C  = q(2,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      u9R  = q(2,i,j,k+1) * zet(3,i,j,k+1) / aj(i,j,k+1) 
      u9RR = q(2,i,j,k+2) * zet(3,i,j,k+2) / aj(i,j,k+2) 

      !-------------------------------------------------------

      v7LL = q(3,i,j,k-2) * zet(1,i,j,k-2) / aj(i,j,k-2) 
      v7L  = q(3,i,j,k-1) * zet(1,i,j,k-1) / aj(i,j,k-1) 
      v7C  = q(3,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      v7R  = q(3,i,j,k+1) * zet(1,i,j,k+1) / aj(i,j,k+1) 
      v7RR = q(3,i,j,k+2) * zet(1,i,j,k+2) / aj(i,j,k+2) 

      v8LL = q(3,i,j,k-2) * zet(2,i,j,k-2) / aj(i,j,k-2) 
      v8L  = q(3,i,j,k-1) * zet(2,i,j,k-1) / aj(i,j,k-1) 
      v8C  = q(3,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      v8R  = q(3,i,j,k+1) * zet(2,i,j,k+1) / aj(i,j,k+1) 
      v8RR = q(3,i,j,k+2) * zet(2,i,j,k+2) / aj(i,j,k+2) 

      v9LL = q(3,i,j,k-2) * zet(3,i,j,k-2) / aj(i,j,k-2) 
      v9L  = q(3,i,j,k-1) * zet(3,i,j,k-1) / aj(i,j,k-1) 
      v9C  = q(3,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      v9R  = q(3,i,j,k+1) * zet(3,i,j,k+1) / aj(i,j,k+1) 
      v9RR = q(3,i,j,k+2) * zet(3,i,j,k+2) / aj(i,j,k+2) 

      !-------------------------------------------------------
   
      w7LL = q(4,i,j,k-2) * zet(1,i,j,k-2) / aj(i,j,k-2) 
      w7L  = q(4,i,j,k-1) * zet(1,i,j,k-1) / aj(i,j,k-1) 
      w7C  = q(4,i,j,k  ) * zet(1,i,j,k  ) / aj(i,j,k  ) 
      w7R  = q(4,i,j,k+1) * zet(1,i,j,k+1) / aj(i,j,k+1) 
      w7RR = q(4,i,j,k+2) * zet(1,i,j,k+2) / aj(i,j,k+2) 

      w8LL = q(4,i,j,k-2) * zet(2,i,j,k-2) / aj(i,j,k-2) 
      w8L  = q(4,i,j,k-1) * zet(2,i,j,k-1) / aj(i,j,k-1) 
      w8C  = q(4,i,j,k  ) * zet(2,i,j,k  ) / aj(i,j,k  ) 
      w8R  = q(4,i,j,k+1) * zet(2,i,j,k+1) / aj(i,j,k+1) 
      w8RR = q(4,i,j,k+2) * zet(2,i,j,k+2) / aj(i,j,k+2) 

      w9LL = q(4,i,j,k-2) * zet(3,i,j,k-2) / aj(i,j,k-2) 
      w9L  = q(4,i,j,k-1) * zet(3,i,j,k-1) / aj(i,j,k-1) 
      w9C  = q(4,i,j,k  ) * zet(3,i,j,k  ) / aj(i,j,k  ) 
      w9R  = q(4,i,j,k+1) * zet(3,i,j,k+1) / aj(i,j,k+1) 
      w9RR = q(4,i,j,k+2) * zet(3,i,j,k+2) / aj(i,j,k+2) 


   end if

   exsign = zero

   ! ∂u/∂ζ
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     u7LL, u7L , u7C, u7R, u7RR   , &
                                     dz                           , &
                                     bias_zet                     , &
                                     exsign                       , &
                                     du7_dzet                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     u8LL, u8L , u8C, u8R, u8RR   , &
                                     dz                           , &
                                     bias_zet                     , &
                                     exsign                       , &
                                     du8_dzet                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     u9LL, u9L , u9C, u9R, u9RR   , &
                                     dz                           , &
                                     bias_zet                     , &
                                     exsign                       , &
                                     du9_dzet                       &
                                    )

   if ( exsign < one_half ) return                                    
                                    


   ! ∂v/∂ζ
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     v7LL, v7L , v7C, v7R, v7RR   , &
                                     dz                           , &
                                     bias_zet                     , &
                                     exsign                       , &
                                     dv7_dzet                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     v8LL, v8L , v8C, v8R, v8RR   , &
                                     dz                           , &
                                     bias_zet                     , &
                                     exsign                       , &
                                     dv8_dzet                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     v9LL, v9L , v9C, v9R, v9RR   , &
                                     dz                           , &
                                     bias_zet                     , &
                                     exsign                       , &
                                     dv9_dzet                       &
                                    )

   if ( exsign < one_half ) return        

   ! ∂w/∂ζ
   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     w7LL, w7L , w7C, w7R, w7RR   , &
                                     dz                           , &
                                     bias_zet                     , &
                                     exsign                       , &
                                     dw7_dzet                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     w8LL, w8L , w8C, w8R, w8RR   , &
                                     dz                           , &
                                     bias_zet                     , &
                                     exsign                       , &
                                     dw8_dzet                       &
                                    )

   if ( exsign < one_half ) return                                    

   call WaterPhase_FirstDerivative ( rLL, rL , rC, rR, rRR        , &
                                     w9LL, w9L , w9C, w9R, w9RR   , &
                                     dz                           , &
                                     bias_zet                     , &
                                     exsign                       , &
                                     dw9_dzet                       &
                                    )

   if ( exsign < one_half ) return                                                                  



   du_dx = aj(i,j,k) * ( du1_dcsi + du4_deta + du7_dzet ) 
   du_dy = aj(i,j,k) * ( du2_dcsi + du5_deta + du8_dzet ) 
   du_dz = aj(i,j,k) * ( du3_dcsi + du6_deta + du9_dzet ) 

   dv_dx = aj(i,j,k) * ( dv1_dcsi + dv4_deta + dv7_dzet ) 
   dv_dy = aj(i,j,k) * ( dv2_dcsi + dv5_deta + dv8_dzet ) 
   dv_dz = aj(i,j,k) * ( dv3_dcsi + dv6_deta + dv9_dzet ) 

   dw_dx = aj(i,j,k) * ( dw1_dcsi + dw4_deta + dw7_dzet ) 
   dw_dy = aj(i,j,k) * ( dw2_dcsi + dw5_deta + dw8_dzet ) 
   dw_dz = aj(i,j,k) * ( dw3_dcsi + dw6_deta + dw9_dzet ) 

   velocity_gradient(1,1) = du_dx
   velocity_gradient(1,2) = du_dy
   velocity_gradient(1,3) = du_dz

   velocity_gradient(2,1) = dv_dx
   velocity_gradient(2,2) = dv_dy
   velocity_gradient(2,3) = dv_dz

   velocity_gradient(3,1) = dw_dx
   velocity_gradient(3,2) = dw_dy
   velocity_gradient(3,3) = dw_dz

end subroutine velocity_curv_gradient_tensor2