subroutine calc_pressure_gradient1(i, j, k, PressureGradient, exsign ) 
   
  ! Input and output variables
  integer, intent(in) :: i,j,k
  real ( kind = rdf ), intent(inout) :: exsign
  real ( kind = rdf ), dimension(3), intent(out) :: PressureGradient

   ! Local variables
   real ( kind = rdf ), dimension(3)   :: PressureCurvilinearGradient
   real ( kind = rdf ), dimension(3,3) :: MetricsTensor
   real ( kind = rdf ) :: dc2, de2, dz2
   integer :: m,l,p
   integer :: ista , jsta , ksta , iend , jend , kend


   ! Physical boundaries

   ista = il ; jsta = jl ; ksta = kl 
   iend = iu ; jend = ju ; kend = ku 

   if ( myback  == mpi_proc_null )  ista = il + igp 
   if ( myleft  == mpi_proc_null )  jsta = jl + jgp 
   if ( mydown  == mpi_proc_null )  ksta = kl + kgp 

   if ( myfront == mpi_proc_null )  iend = iu - igp
   if ( myright == mpi_proc_null )  jend = ju - jgp
   if ( myup    == mpi_proc_null )  kend = ku - kgp

   PressureGradient            = zero
   PressureCurvilinearGradient = zero
   MetricsTensor               = zero

   dc2 = one_half * dc
   de2 = one_half * de
   dz2 = one_half * dz 

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                         
   !                                      _         _                                          
   !                                     |   ∂p/∂ξ   |  
   !       PressureCurvilinearGradient = |   ∂p/∂η   |                               
   !                                     |   ∂p/∂ζ   |  
   !                                      •-        -• 
   !       We calculate the derivatives using a second order centred
   !       difference scheme. If the stencil is bounded, then a second
   !       order biased derivative is applied.
   ! 
   !       It assumes the routine that calls it has ista, jsta, ksta and
   !       iend, jend, kend, defined.
   !
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                         
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                         
   !                            ξ - DIRECTION
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                         

   if ( myback == mpi_proc_null .and. i == ista ) then

      ! ∂p/∂ξ
      PressureCurvilinearGradient(1)   =  dc2 * (  - one   * q(1,i+2,j,k) &
                                                   + four  * q(1,i+1,j,k) &
                                                   - three * q(1,i  ,j,k) )

      ! if rsign == 1 or -1, then exsign is 1
      exsign = abs( rsign(i+2,j,k) * rsign(i+1,j,k) * rsign(i,j,k) )

      if ( exsign < one_half ) return                                    

   else if ( myfront == mpi_proc_null .and. i == iend) then

      ! ∂p/∂ξ
      PressureCurvilinearGradient(1)   =  dc2 * (    one   * q(1,i-2,j,k) &
                                                   - four  * q(1,i-1,j,k) &
                                                   + three * q(1,i  ,j,k) )

      exsign = abs( rsign(i-2,j,k) * rsign(i-1,j,k) * rsign(i,j,k) )

      if ( exsign < one_half ) return                                    

   else

      ! ∂p/∂ξ
      PressureCurvilinearGradient(1) = dc2 * ( q(1,i+1,j,k) - q(1,i-1,j,k) )

      exsign = abs( rsign(i+1,j,k) * rsign(i-1,j,k) )

      if ( exsign < one_half ) return                                    

   end if
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                         
   !                            η - DIRECTION
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                         

   if ( myleft == mpi_proc_null .and. j == jsta ) then

      ! ∂p/∂η
      PressureCurvilinearGradient(2)   =  de2 * (  - one   * q(1,i,j+2,k) &
                                                   + four  * q(1,i,j+1,k) &
                                                   - three * q(1,i,j  ,k) )

      ! if rsign == 1 or -1, then exsign is 1
      exsign = abs( rsign(i,j+2,k) * rsign(i,j+1,k) * rsign(i,j,k) )

      if ( exsign < one_half ) return                                    

   else if ( myright == mpi_proc_null .and. j == jend) then

      ! ∂p/∂η
      PressureCurvilinearGradient(2)   =  de2 * (    one   * q(1,i,j-2,k) &
                                                   - four  * q(1,i,j-1,k) &
                                                   + three * q(1,i,j  ,k) )
      
      ! if rsign == 1 or -1, then exsign is 1
      exsign = abs( rsign(i,j-2,k) * rsign(i,j-1,k) * rsign(i,j,k) )

      if ( exsign < one_half ) return                                    

   else

      ! ∂p/∂η
      PressureCurvilinearGradient(2) = de2 * ( q(2,i,j+1,k) - q(2,i,j-1,k) )

      ! if rsign == 1 or -1, then exsign is 1
      exsign = abs( rsign(i,j+1,k) * rsign(i,j-1,k)  )

      if ( exsign < one_half ) return                                    

   end if


   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                         
   !                            ζ - DIRECTION
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                         

   if ( mydown == mpi_proc_null .and. k == ksta ) then

      ! ∂p/∂ζ
      PressureCurvilinearGradient(3)   =  dz2 * (  - one   * q(1,i,j,k+2) &
                                                   + four  * q(1,i,j,k+1) &
                                                   - three * q(1,i,j,k  ) )

      ! if rsign == 1 or -1, then exsign is 1
      exsign = abs( rsign(i,j,k+2) * rsign(i,j,k+1) * rsign(i,j,k) )

      if ( exsign < one_half ) return                                    

   else if ( myup == mpi_proc_null .and. k == kend) then

      ! ∂p/∂ζ
      PressureCurvilinearGradient(3)   =  dz2 * (    one   * q(1,i,j,k-2) &
                                                   - four  * q(1,i,j,k-1) &
                                                   + three * q(1,i,j,k  ) )

      ! if rsign == 1 or -1, then exsign is 1
      exsign = abs( rsign(i,j,k-2) * rsign(i,j,k-1) * rsign(i,j,k) )

      if ( exsign < one_half ) return                                    

   else

      ! ∂p/∂ζ
      PressureCurvilinearGradient(3) = dz2 * ( q(1,i,j,k+1) - q(1,i,j,k-1) )

      ! if rsign == 1 or -1, then exsign is 1
      exsign = abs( rsign(i,j,k+1) * rsign(i,j,k-1) )

      if ( exsign < one_half ) return                                    

   end if


   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                         
   !                                _                    _                                          
   !                               | ∂ξ/∂x  ∂ξ/∂y  ∂ξ/∂z  |  
   !               MetricsTensor = | ∂η/∂x  ∂η/∂y  ∂η/∂z  |                               
   !                               | ∂ζ/∂x  ∂ζ/∂y  ∂ζ/∂z  |  
   !                               •-                    -•  
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                         

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
       
       PressureGradient(m) =   PressureGradient(m) &
                             + PressureCurvilinearGradient(p) * MetricsTensor(p,m)
    end do
   end do

end subroutine calc_pressure_gradient1