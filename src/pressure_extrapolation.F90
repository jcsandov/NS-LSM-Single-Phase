subroutine pressure_extrapolation( )

  use InterpolationMethods

  implicit none

  ! local variables

  integer :: i,j,k
  integer :: ii, jj, kk
  integer :: isearch, jsearch, ksearch
  integer :: iic, jic, kic
  integer :: is, js, ks
  integer :: ioc, joc, koc
  integer :: iaux, jaux, kaux

   !indexes
   integer  :: i_mysta , i_myend
   integer  :: j_mysta , j_myend
   integer  :: k_mysta , k_myend

   integer  :: ista , iend
   integer  :: jsta , jend
   integer  :: ksta , kend

  real(kind=rdf) :: rAir2Water(3)

  real(kind=rdf) :: PhiAir, PhiWater, pWater , pmaxmod
  real(kind=rdf) :: xIntpNode, yIntpNode, zIntpNode
  real(kind=rdf) :: dWater, dAir, dTotal
  real(kind=rdf) :: pAux, LocalWeight, TotalWeight
  real(kind=rdf) :: dx   
  real(kind=rdf) , parameter :: extp_threshold = 1.0E-05
  real(kind=rdf) , parameter :: pmax = two

  integer :: nWaterNodesExtp


   !Nodos incluyendo el borde
   i_mysta = il + igp
   j_mysta = jl + jgp
   k_mysta = kl + kgp

   i_myend = iu - igp
   j_myend = ju - jgp
   k_myend = ku - kgp
         
   ! processes on the domain boundaries

   if ( myback  == mpi_proc_null )  i_mysta = il + igp + 1
   if ( myleft  == mpi_proc_null )  j_mysta = jl + jgp + 1
   if ( mydown  == mpi_proc_null )  k_mysta = kl + kgp + 1
   
   if ( myfront == mpi_proc_null )  i_myend = iu - igp - 1
   if ( myright == mpi_proc_null )  j_myend = ju - jgp - 1
   if ( myup    == mpi_proc_null )  k_myend = ku - kgp - 1


   ! Physical boundaries of the processor
   
   ista = il ; jsta = jl ; ksta = kl 
   iend = iu ; jend = ju ; kend = ku 

   if ( myback  == mpi_proc_null )  ista = il + igp 
   if ( myleft  == mpi_proc_null )  jsta = jl + jgp 
   if ( mydown  == mpi_proc_null )  ksta = kl + kgp 

   if ( myfront == mpi_proc_null )  iend = iu - igp
   if ( myright == mpi_proc_null )  jend = ju - jgp
   if ( myup    == mpi_proc_null )  kend = ku - kgp


  pmaxmod = zero
   
  ! Loop which considers the boundaries
  do k = k_mysta - 1 , k_myend + 1 !ksta , kend !k_mysta - 1 , k_myend + 1
  do j = j_mysta - 1 , j_myend + 1 !jsta , jend !j_mysta - 1 , j_myend + 1
  do i = i_mysta - 1 , i_myend + 1 !ista , iend !i_mysta - 1 , i_myend + 1
            
    ! ----------------------------------------------------------------------------
    ! 1. Identify if the node i,j,k has to be extrapolated based on rsign value
    ! ----------------------------------------------------------------------------
    ! if rsign(i,j,k) = 0 (air-phase), but one of its neighbors is in the water
    ! phase (rsign = 1), then this if is true. This condition identifies air nodes
    ! next to the free-surface
    ! ----------------------------------------------------------------------------

    ! If the node is closer to the free-surface than the machine-precision, then
    ! I'm on the free surface and the pressure is zero  
    if ( abs( phi(i,j,k) ) < eps_sims ) then

      q(1,i,j,k)   = zero
      rsign(i,j,k) = -one

      cycle

    end if

    ! local lenght scale
    dx = aj(i,j,k)**( -one_third )

    ! If I'm in the water phase or too far from the fs, I just cycle
    !if ( rsign(i,j,k) > one_half .or. abs( phi(i,j,k) ) > 3.5 * dx ) cycle 
    if ( rsign(i,j,k) > one_half .or. abs( phi(i,j,k) ) > 10.0 * dx ) cycle 
               
    ! Check if there are water or near water nodes in the neighbourhood 
    nWaterNodesExtp = 0

    searchloop2 : &
    do kk = max( k-1 , ksta ) , min( k+1 , kend )
    do jj = max( j-1 , jsta ) , min( j+1 , jend )
    do ii = max( i-1 , ista ) , min( i+1 , iend )
      
      if ( rsign (ii,jj,kk) > one_half ) then
            
        nWaterNodesExtp = 1
            
        exit searchloop2
            
      end if

    end do
    end do
    end do searchloop2

    if ( nWaterNodesExtp == 0 ) cycle

    ! variables initialisation
    pAux        = zero
    TotalWeight = zero
    PhiAir      = abs( phi(i,j,k) )

    do kk = max( k-2 , ksta ) , min( k+2 , kend )
    do jj = max( j-2 , jsta ) , min( j+2 , jend )
    do ii = max( i-2 , ista ) , min( i+2 , iend )

      ! If the node is too close to the fs, I don't use
      ! if for interpolation. Also, if the gradient is too large, maybe
      ! I'm taking nodes away from the geometric reinitialisation range,
      ! so I discard those nodes too.
      !if ( ( ii == i .and. jj == j .and. kk == k )                 .or.  &
      !       phi(ii,jj,kk)                       < extp_threshold  .or.  & ! ϕ not deep enough
      !       q(1,ii,jj,kk)                       < extp_threshold  .or.  & ! water -ve pressure
      !       norm2( phi_gradient(1:3,ii,jj,kk) ) > two             .or.  & ! ∇ϕ badly defined
      !       norm2( phi_gradient(1:3,ii,jj,kk) ) < pt1                    ) cycle

      if ( ( ii == i .and. jj == j .and. kk == k )                 .or.  &
             phi(ii,jj,kk)                       < extp_threshold  .or.  & ! ϕ not deep enough
             norm2( phi_gradient(1:3,ii,jj,kk) ) > two             .or.  & ! ∇ϕ badly defined
             norm2( phi_gradient(1:3,ii,jj,kk) ) < pt1                    ) cycle

      PhiWater = phi(ii,jj,kk)

      pWater   = q(1,ii,jj,kk)

      if ( abs(pWater) > abs(pmaxmod) ) pmaxmod = -pWater

      rAir2Water = (/ x(ii,jj,kk) - x(i,j,k) , &
                      y(ii,jj,kk) - y(i,j,k) , &
                      z(ii,jj,kk) - z(i,j,k)    /)

      dTotal = norm2( rAir2Water )

      ! SmoothStepFunction returns always a number between 0 and 1
      dWater = PhiWater  / ( PhiWater + PhiAir  ) * dTotal ! SmoothStepFunction( PhiWater / ( PhiWater + PhiAir ) ) * dTotal 
      dAir   = dTotal - dWater

      ! If dWater is 0, then paux is undetermined
      if ( dWater < eps_sims ) cycle
      
      ! If the locally extrapolated pressure is too big, then it's not considered
      if ( abs( ( dAir /  dWater ) * pWater ) > pmax ) cycle

      LocalWeight = abs(  dot_product ( rAir2Water  / norm2( rAir2Water ) , &
                          phi_gradient(1:3,ii,jj,kk) / norm2( phi_gradient(1:3,ii,jj,kk) ) ) )!**two

      pAux        = pAux - LocalWeight * ( dAir /  dWater ) * pWater 
      TotalWeight = TotalWeight + LocalWeight

      !if ( i==32 .and. j==100 .and. k==73 ) then
      !  
      !  print * , ' '
      !  print * , 'ii,jj,kk    = ' , ii,jj,kk
      !  print * , 'PhiWater    = ' , PhiWater
      !  print * , 'pWater      = ' , pWater
      !  print * , 'rAir2Water  = ' , rAir2Water
      !  print * , 'dTotal      = ' , dTotal
      !  print * , 'dWater      = ' , dWater
      !  print * , 'dAir        = ' , dAir
      !  print * , 'LocalWeight = ' , LocalWeight
      !  print * , 'TotalWeight = ' , TotalWeight
      !  print * , 'pAux        = ' , LocalWeight * ( dAir /  dWater ) * pWater
      !  print * , 'pAux accum  = ' , pAux
      !  print * , ' '
      !  
      !end if 

    end do ! ii = max( i-2 , ista ) , min( i+2 , iend )
    end do ! jj = max( j-2 , jsta ) , min( j+2 , jend )
    end do ! kk = max( k-2 , ksta ) , min( k+2 , kend )

    if ( abs( TotalWeight ) > eps_sims ) then 
      
      ! DEBUG:
      if ( abs(pAux / TotalWeight) < abs(pmaxmod) ) then
        q(1,i,j,k)   = pAux / TotalWeight
      else 
        q(1,i,j,k) = pmaxmod
      end if

      rsign(i,j,k) = -one
      
      !if ( i==32 .and. j==100 .and. k==73 ) then
      !  print * , ' '
      !  print * , 'q(1,i,j,k)            = ' , q(1,i,j,k)
      !  print * , 'h(i,j,k)              = ' , h(i,j,k)
      !  print * , 'q(1,i,j,k) - h(i,j,k) = ' , q(1,i,j,k) - h(i,j,k) 
      !  print * , ' '
      !end if
      
      !q(1,i,j,k) = two * q(1,i,j,k-1) - q(1,i,j,k-2)  
      
      
    else

      q(1,i,j,k)   = two * q(1,i,j,k-1) - q(1,i,j,k-2)
      rsign(i,j,k) = -one

    end if

  end do ! i = i_mysta - 1 , i_myend + 1
  end do ! j = j_mysta - 1 , j_myend + 1
  end do ! k = k_mysta - 1 , k_myend + 1


end subroutine pressure_extrapolation