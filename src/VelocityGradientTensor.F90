subroutine VelocityGradientTensor( i , j , k , VelocityGradient )

   ! Input and output variables
   integer, intent(in) :: i,j,k
   real ( kind = rdf ), dimension(3,3), intent(out) :: VelocityGradient

   ! Local variables
   real ( kind = rdf ), dimension(3,3) :: VelocityCurvilinearGradient, MetricsTensor
   real ( kind = rdf ) :: dc2, de2, dz2
   integer :: m,l,p
   integer :: ista , jsta , ksta , iend , jend , kend

   logical :: BlankingFlagFront , BlankingFlagBack , BlankingFlagLeft , BlankingFlagRight

   ! Physical boundaries

   ista = il ; jsta = jl ; ksta = kl 
   iend = iu ; jend = ju ; kend = ku 

   if ( myback  == mpi_proc_null )  ista = il + igp 
   if ( myleft  == mpi_proc_null )  jsta = jl + jgp 
   if ( mydown  == mpi_proc_null )  ksta = kl + kgp 

   if ( myfront == mpi_proc_null )  iend = iu - igp
   if ( myright == mpi_proc_null )  jend = ju - jgp
   if ( myup    == mpi_proc_null )  kend = ku - kgp

   VelocityGradient            = zero
   VelocityCurvilinearGradient = zero
   MetricsTensor               = zero

   dc2 = one_half * dc
   de2 = one_half * de
   dz2 = one_half * dz 

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                         
   !                                      _                    _                                          
   !                                     | ∂u/∂ξ  ∂u/∂η  ∂u/∂ζ  |  
   !       VelocityCurvilinearGradient = | ∂v/∂ξ  ∂v/∂η  ∂v/∂ζ  |                               
   !                                     | ∂w/∂ξ  ∂w/∂η  ∂w/∂ζ  |  
   !                                      •-                   -• 
   !       We calculate the derivatives using a second order centred
   !       difference scheme. If the stencil is bounded, then a second
   !       order biased derivative is applied.
   ! 
   !       It assumes the routine that calls it has ista, jsta, ksta and
   !       iend, jend, kend, defined.
   !
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                         
      
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

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                         
   !                            ξ - DIRECTION
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                         

   if ( ( myback == mpi_proc_null .and. i == ista ) .or. BlankingFlagBack ) then

      ! ∂u/∂ξ
      VelocityCurvilinearGradient(1,1)   =  dc2 * (  - one   * q(2,i+2,j,k) &
                                                     + four  * q(2,i+1,j,k) &
                                                     - three * q(2,i  ,j,k) )

      ! ∂v/∂ξ
      VelocityCurvilinearGradient(2,1)   =  dc2 * (  - one   * q(3,i+2,j,k) &
                                                     + four  * q(3,i+1,j,k) &
                                                     - three * q(3,i  ,j,k) )

      ! ∂w/∂ξ
      VelocityCurvilinearGradient(3,1)   =  dc2 * (  - one   * q(4,i+2,j,k) &
                                                     + four  * q(4,i+1,j,k) &
                                                     - three * q(4,i  ,j,k) )

   else if ( ( myfront == mpi_proc_null .and. i == iend ) .or. BlankingFlagFront ) then

      ! ∂u/∂ξ
      VelocityCurvilinearGradient(1,1)   =  dc2 * (    one   * q(2,i-2,j,k) &
                                                     - four  * q(2,i-1,j,k) &
                                                     + three * q(2,i  ,j,k) )

      ! ∂v/∂ξ
      VelocityCurvilinearGradient(2,1)   =  dc2 * (    one   * q(3,i-2,j,k) &
                                                     - four  * q(3,i-1,j,k) &
                                                     + three * q(3,i  ,j,k) )

      ! ∂w/∂ξ
      VelocityCurvilinearGradient(3,1)   =  dc2 * (    one   * q(4,i-2,j,k) &
                                                     - four  * q(4,i-1,j,k) &
                                                     + three * q(4,i  ,j,k) )

   else

      ! ∂u/∂ξ
      VelocityCurvilinearGradient(1,1) = dc2 * ( q(2,i+1,j,k) - q(2,i-1,j,k) )
      ! ∂v/∂ξ
      VelocityCurvilinearGradient(2,1) = dc2 * ( q(3,i+1,j,k) - q(3,i-1,j,k) )
      ! ∂w/∂ξ
      VelocityCurvilinearGradient(3,1) = dc2 * ( q(4,i+1,j,k) - q(4,i-1,j,k) )

   end if
   
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                         
   !                            η - DIRECTION
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                         

   if ( ( myleft == mpi_proc_null .and. j == jsta ) .or. BlankingFlagRight ) then

      ! ∂u/∂η
      VelocityCurvilinearGradient(1,2)   =  de2 * (  - one   * q(2,i,j+2,k) &
                                                     + four  * q(2,i,j+1,k) &
                                                     - three * q(2,i,j  ,k) )

      ! ∂v/∂η
      VelocityCurvilinearGradient(2,2)   =  de2 * (  - one   * q(3,i,j+2,k) &
                                                     + four  * q(3,i,j+1,k) &
                                                     - three * q(3,i,j  ,k) )

      ! ∂w/∂η
      VelocityCurvilinearGradient(3,2)   =  de2 * (  - one   * q(4,i,j+2,k) &
                                                     + four  * q(4,i,j+1,k) &
                                                     - three * q(4,i,j  ,k) )

   else if ( ( myright == mpi_proc_null .and. j == jend ) .or. BlankingFlagLeft ) then

      ! ∂u/∂η
      VelocityCurvilinearGradient(1,2)   =  de2 * (    one   * q(2,i,j-2,k) &
                                                     - four  * q(2,i,j-1,k) &
                                                     + three * q(2,i,j  ,k) )

      ! ∂v/∂η
      VelocityCurvilinearGradient(2,2)   =  de2 * (    one   * q(3,i,j-2,k) &
                                                     - four  * q(3,i,j-1,k) &
                                                     + three * q(3,i,j  ,k) )

      ! ∂w/∂η
      VelocityCurvilinearGradient(3,2)   =  de2 * (    one   * q(4,i,j-2,k) &
                                                     - four  * q(4,i,j-1,k) &
                                                     + three * q(4,i,j  ,k) )

   else

      ! ∂u/∂η
      VelocityCurvilinearGradient(1,2) = de2 * ( q(2,i,j+1,k) - q(2,i,j-1,k) )
      ! ∂v/∂η
      VelocityCurvilinearGradient(2,2) = de2 * ( q(3,i,j+1,k) - q(3,i,j-1,k) )
      ! ∂w/∂η
      VelocityCurvilinearGradient(3,2) = de2 * ( q(4,i,j+1,k) - q(4,i,j-1,k) )

   end if


   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                         
   !                            ζ - DIRECTION
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                         

   if ( mydown == mpi_proc_null .and. k == ksta ) then

      ! ∂u/∂ζ
      VelocityCurvilinearGradient(1,3)   =  dz2 * (  - one   * q(2,i,j,k+2) &
                                                     + four  * q(2,i,j,k+1) &
                                                     - three * q(2,i,j,k  ) )

      ! ∂v/∂ζ
      VelocityCurvilinearGradient(2,3)   =  dz2 * (  - one   * q(3,i,j,k+2) &
                                                     + four  * q(3,i,j,k+1) &
                                                     - three * q(3,i,j,k  ) )

      ! ∂w/∂ζ
      VelocityCurvilinearGradient(3,3)   =  dz2 * (  - one   * q(4,i,j,k+2) &
                                                     + four  * q(4,i,j,k+1) &
                                                     - three * q(4,i,j,k  ) )

   else if ( myup == mpi_proc_null .and. k == kend) then

      ! ∂u/∂ζ
      VelocityCurvilinearGradient(1,3)   =  dz2 * (    one   * q(2,i,j,k-2) &
                                                     - four  * q(2,i,j,k-1) &
                                                     + three * q(2,i,j,k  ) )

      ! ∂v/∂ζ
      VelocityCurvilinearGradient(2,3)   =  dz2 * (    one   * q(3,i,j,k-2) &
                                                     - four  * q(3,i,j,k-1) &
                                                     + three * q(3,i,j,k  ) )

      ! ∂w/∂ζ
      VelocityCurvilinearGradient(3,3)   =  dz2 * (    one   * q(4,i,j,k-2) &
                                                     - four  * q(4,i,j,k-1) &
                                                     + three * q(4,i,j,k  ) )

   else

      ! ∂u/∂ζ
      VelocityCurvilinearGradient(1,3) = dz2 * ( q(2,i,j,k+1) - q(2,i,j,k-1) )
      ! ∂v/∂ζ
      VelocityCurvilinearGradient(2,3) = dz2 * ( q(3,i,j,k+1) - q(3,i,j,k-1) )
      ! ∂w/∂ζ
      VelocityCurvilinearGradient(3,3) = dz2 * ( q(4,i,j,k+1) - q(4,i,j,k-1) )

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

   ! VelocityGradient(m,l) = ∂um/∂xl = ( ∂um/∂ξp ) * ( ∂ξp/∂xl )

   do m = 1,3
   do l = 1,3
      do p = 1,3
         
         VelocityGradient(m,l) =   VelocityGradient(m,l) &
                                 + VelocityCurvilinearGradient(m,p) * MetricsTensor(p,l)
      end do
   end do
   end do

end subroutine VelocityGradientTensor
