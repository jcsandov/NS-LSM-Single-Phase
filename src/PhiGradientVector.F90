subroutine PhiGradientVector( i , j , k , PhiGradient )

   ! Input and output variables
   integer, intent(in) :: i,j,k
   real ( kind = rdf ), dimension(3), intent(out) :: PhiGradient
   real ( kind = rdf ) :: dc2, de2, dz2

   ! Local variables
   real ( kind = rdf ), dimension(3)   :: PhiCurvilinearGradient
   real ( kind = rdf ), dimension(3,3) :: MetricsTensor
   integer :: m,l,p

   dc2 = one_half * dc
   de2 = one_half * de
   dz2 = one_half * dz 

   PhiGradient            = zero
   PhiCurvilinearGradient = zero
   MetricsTensor          = zero

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                       
   !
   !  ∂ϕ     ∂ϕ     ∂ξ^k    
   ! ---- = ---- * ------ 
   ! ∂x_j   ∂ξ^k    ∂x_j  
   !
   !                           _     _                                          
   !                          | ∂ϕ/∂ξ |  
   ! PhiCurvilinearGradient = | ∂ϕ/∂η |                               
   !                          | ∂ϕ/∂ζ |  
   !                           -     - 
   !
   ! We calculate the derivatives using a second order centred
   ! difference scheme. If the stencil is bounded, then a second
   ! order biased derivative is applied.
   ! 
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                       
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                       
   !                          ξ - DIRECTION
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                       

   if ( myback == mpi_proc_null .and. i == il ) then

      ! ∂ϕ/∂ξ
      PhiCurvilinearGradient(1)   =  dc2 * (  - one   * phi( i+2,j,k ) &
                                              + four  * phi( i+1,j,k ) &
                                              - three * phi( i  ,j,k ) )

   else if ( myfront == mpi_proc_null .and. i == iu) then

      ! ∂ϕ/∂ξ
      PhiCurvilinearGradient(1)   =  dc2 * (    one   * phi( i-2,j,k ) &
                                              - four  * phi( i-1,j,k ) &
                                              + three * phi( i  ,j,k ) )
   else

      ! ∂ϕ/∂ξ
      PhiCurvilinearGradient(1)   = dc2 * ( phi(i+1,j,k) - phi(i-1,j,k) )

   end if
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                       
   !                          η - DIRECTION
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                       

   if ( myleft == mpi_proc_null .and. j == jl ) then

      ! ∂ϕ/∂η
      PhiCurvilinearGradient(2)   =  de2 * (  - one   * phi( i,j+2,k ) &
                                              + four  * phi( i,j+1,k ) &
                                              - three * phi( i,j  ,k ) )

   else if ( myright == mpi_proc_null .and. j == ju) then

      ! ∂ϕ/∂η
      PhiCurvilinearGradient(2)   =  de2 * (    one   * phi( i,j-2,k ) &
                                              - four  * phi( i,j-1,k ) &
                                              + three * phi( i,j  ,k ) )

   else

      ! ∂ϕ/∂η
      PhiCurvilinearGradient(2) = de2 * ( phi(i,j+1,k) - phi(i,j-1,k) )

   end if

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                       
   !                          ζ - DIRECTION
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                       

   if ( mydown == mpi_proc_null .and. k == kl ) then

      ! ∂ϕ/∂ζ
      PhiCurvilinearGradient(3)   =  dz2 * (  - one   * phi( i,j,k+2) &
                                              + four  * phi( i,j,k+1) &
                                              - three * phi( i,j,k  ) )

   else if ( myright == mpi_proc_null .and. j == ju) then

      ! ∂ϕ/∂ζ
      PhiCurvilinearGradient(3)   =  dz2 * (    one   * phi( i,j,k-2 ) &
                                              - four  * phi( i,j,k-1 ) &
                                              + three * phi( i,j,k   ) )
   else

      ! ∂ϕ/∂ζ
      PhiCurvilinearGradient(3) = dz2 * ( phi(i,j,k+1) - phi(i,j,k-1) )

   end if


   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                       
   !                             _                    _                                          
   !                            | ∂ξ/∂x  ∂ξ/∂y  ∂ξ/∂z  |  
   !            MetricsTensor = | ∂η/∂x  ∂η/∂y  ∂η/∂z  |                               
   !                            | ∂ζ/∂x  ∂ζ/∂y  ∂ζ/∂z  |  
   !                            •-                    -•  
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                       

   MetricsTensor(1,1) = csi(1,i,j,k)
   MetricsTensor(1,2) = csi(2,i,j,k)
   MetricsTensor(1,3) = csi(3,i,j,k)

   MetricsTensor(2,1) = eta(1,i,j,k)
   MetricsTensor(2,2) = eta(2,i,j,k)
   MetricsTensor(2,3) = eta(3,i,j,k)

   MetricsTensor(3,1) = zet(1,i,j,k)
   MetricsTensor(3,2) = zet(2,i,j,k)
   MetricsTensor(3,3) = zet(3,i,j,k)

   ! PhiGradient(m) = ∂ϕ/∂xm = ( ∂ϕ/∂ξp ) * ( ∂ξp/∂xm )

   do m = 1,3
      do p = 1,3
         
         PhiGradient(m) =    PhiGradient(m) &
                           + PhiCurvilinearGradient(p) * MetricsTensor(p,m)
      end do
   end do

end subroutine PhiGradientVector
