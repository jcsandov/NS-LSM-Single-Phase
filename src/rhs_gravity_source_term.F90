
subroutine rhs_gravity_source_term
  
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   ! It adds the gravity source term to the rhs
   !
   ! output
   !     rh(4,ijk)

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   ! local dummy reals
   !

   real (kind = rdf), dimension( il:iu , jl:ju , kl:ku ) :: GravitySource

   !DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
   real(kind = rdf) :: hgradx

   GravitySource = zero

   if ( hydraulic_mode ) then

      do k = k_mysta , k_myend
      do j = j_mysta , j_myend
      do i = i_mysta , i_myend
      
         rh(2,i,j,k) = rh(2,i,j,k) + h_gradient(1,i,j,k) / ( aj(i,j,k) * FrLSM**two ) 
         rh(3,i,j,k) = rh(3,i,j,k) + h_gradient(2,i,j,k) / ( aj(i,j,k) * FrLSM**two ) 

         !rh(2,i,j,k) = rh(2,i,j,k) + hgradx / ( aj(i,j,k) * FrLSM**two ) 
         !rh(3,i,j,k) = rh(3,i,j,k) + zero * h_gradient(2,i,j,k) / ( aj(i,j,k) * FrLSM**two ) 
   
      end do
      end do
      end do

   else

      do k = k_mysta , k_myend
      do j = j_mysta , j_myend
      do i = i_mysta , i_myend
   
         GravitySource(i,j,k) =  one / ( aj(i,j,k) * FrLSM**two )
   
         rh(4,i,j,k) = rh(4,i,j,k) + GravitySource(i,j,k)
   
      end do
      end do
      end do


   end if

   !write(debugname, fmt ='(a,i6.6)') 'GravitySource', iteraciontiempo
   !call outputD3_real( GravitySource, debugname ) 

end subroutine rhs_gravity_source_term









